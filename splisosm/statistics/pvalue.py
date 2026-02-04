"""P-value computation for HSIC statistics.

Implements multiple approximation methods for the weighted chi-squared
mixture distribution: T → Σᵢⱼ λᵢμⱼzᵢⱼ²

Methods:
- Liu's method: Noncentral chi-squared approximation (recommended, matches official SPLISOSM)
- Gamma approximation: Simpler moment-matching approach
- Permutation testing: Exact but slow
"""

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy import stats


@dataclass
class PValueResult:
    """Result of p-value computation.

    Attributes
    ----------
    pvalue : float
        P-value
    statistic : float
        Test statistic
    method : str
        Method used for p-value computation
    mean : float
        Expected value under null
    variance : float
        Variance under null
    df : float
        Degrees of freedom (for Liu's method)
    delta : float
        Noncentrality parameter (for Liu's method)
    """
    pvalue: float
    statistic: float
    method: str = "liu"
    mean: float = 0.0
    variance: float = 0.0
    df: Optional[float] = None
    delta: Optional[float] = None


def liu_sf(
    t: float,
    eigenvalues: NDArray,
    dofs: Optional[NDArray] = None,
    deltas: Optional[NDArray] = None
) -> Tuple[float, float, float]:
    """Liu's approximation for linear combination of chi-squared variables.

    Approximates: X = Σ λᵢ χ²(hᵢ, δᵢ) using χ²(l, δ)

    Based on:
    - Liu et al. (2009) "A new chi-square approximation to the distribution
      of non-negative definite quadratic forms in non-central normal variables"
    - Implementation follows official SPLISOSM likelihood.py

    Parameters
    ----------
    t : float
        Test statistic value
    eigenvalues : ndarray
        Weights λᵢ for the chi-squared mixture
    dofs : ndarray, optional
        Degrees of freedom hᵢ for each component (default: all 1)
    deltas : ndarray, optional
        Noncentrality parameters δᵢ (default: all 0, central chi-squared)

    Returns
    -------
    pvalue : float
        Survival probability P(X > t)
    df : float
        Approximating degrees of freedom
    delta : float
        Approximating noncentrality parameter
    """
    eigenvalues = np.asarray(eigenvalues).flatten()

    # Filter out near-zero eigenvalues
    eigenvalues = eigenvalues[np.abs(eigenvalues) > 1e-10]

    if len(eigenvalues) == 0:
        return 1.0, 1.0, 0.0

    n_terms = len(eigenvalues)

    # Default: central chi-squared with df=1
    if dofs is None:
        dofs = np.ones(n_terms)
    else:
        dofs = np.asarray(dofs).flatten()

    if deltas is None:
        deltas = np.zeros(n_terms)
    else:
        deltas = np.asarray(deltas).flatten()

    # Compute cumulants c_1, c_2, c_3, c_4
    # c_r = Σᵢ λᵢʳ (hᵢ + r·δᵢ)
    c1 = np.sum(eigenvalues * (dofs + deltas))
    c2 = np.sum(eigenvalues**2 * (dofs + 2 * deltas))
    c3 = np.sum(eigenvalues**3 * (dofs + 3 * deltas))
    c4 = np.sum(eigenvalues**4 * (dofs + 4 * deltas))

    # Standardized cumulants
    # s1 = c3 / c2^{3/2} (skewness-like)
    # s2 = c4 / c2^2 (kurtosis-like)
    if c2 <= 0:
        return 1.0, 1.0, 0.0

    s1 = c3 / (c2 ** 1.5)
    s2 = c4 / (c2 ** 2)

    # Liu's approximation parameters
    # Choose between two approximations based on s1² vs s2
    if s1 ** 2 > s2:
        # Case 1: Use skewness-based approximation
        a = 1.0 / (s1 - np.sqrt(s1**2 - s2))
        delta_approx = s1 * a**3 - a**2
        df_approx = a**2 - 2 * delta_approx
    else:
        # Case 2: Simpler approximation when s1² ≤ s2
        a = 1.0 / s1
        delta_approx = 0.0
        df_approx = 1.0 / (s1**2)

    # Ensure valid parameters
    df_approx = max(df_approx, 0.5)
    delta_approx = max(delta_approx, 0.0)

    # Mean and std of the approximating distribution
    mu_x = c1
    sigma_x = np.sqrt(2 * c2)

    # Mean and std of χ²(l, δ)
    mu_chi = df_approx + delta_approx
    sigma_chi = np.sqrt(2 * (df_approx + 2 * delta_approx))

    # Standardize and transform test statistic
    # t* = (t - μ_X) / σ_X
    # t_final = t* × σ_χ + μ_χ
    if sigma_x > 0:
        t_star = (t - mu_x) / sigma_x
        t_final = t_star * sigma_chi + mu_chi
    else:
        t_final = t

    # Compute p-value using noncentral chi-squared survival function
    if delta_approx < 1e-10:
        # Use central chi-squared for numerical stability
        pvalue = stats.chi2.sf(t_final, df_approx)
    else:
        pvalue = stats.ncx2.sf(t_final, df_approx, delta_approx)

    # Clamp to valid range
    pvalue = np.clip(pvalue, 1e-300, 1.0)

    return float(pvalue), float(df_approx), float(delta_approx)


def liu_pvalue(
    statistic: float,
    eigenvalues: NDArray,
    n: int,
    clamp_min: float = 1e-300,
    clamp_max: float = 1.0
) -> PValueResult:
    """Compute p-value using Liu's method for HSIC.

    This is the recommended method, matching the official SPLISOSM implementation.

    Parameters
    ----------
    statistic : float
        HSIC test statistic (scaled by n)
    eigenvalues : ndarray
        Product eigenvalues λᵢμⱼ from centered kernels
    n : int
        Sample size
    clamp_min : float
        Minimum p-value
    clamp_max : float
        Maximum p-value

    Returns
    -------
    result : PValueResult
    """
    # Scale statistic for asymptotic distribution
    scaled_stat = n * statistic

    pvalue, df, delta = liu_sf(scaled_stat, eigenvalues)
    pvalue = np.clip(pvalue, clamp_min, clamp_max)

    # Compute moments for diagnostics
    eigenvalues = np.asarray(eigenvalues).flatten()
    eigenvalues = eigenvalues[np.abs(eigenvalues) > 1e-10]
    mean = np.sum(eigenvalues) if len(eigenvalues) > 0 else 0.0
    variance = 2 * np.sum(eigenvalues**2) if len(eigenvalues) > 0 else 0.0

    return PValueResult(
        pvalue=pvalue,
        statistic=statistic,
        method="liu",
        mean=mean / n,
        variance=variance / (n**2),
        df=df,
        delta=delta
    )


def compute_composite_eigenvalues(
    eigenvalues_x: NDArray,
    eigenvalues_y: NDArray,
    threshold: float = 1e-10
) -> NDArray:
    """Compute composite eigenvalues for HSIC null distribution.

    The null distribution of HSIC involves the product of eigenvalues
    from both kernel matrices: λ_xy = {λᵢμⱼ : i,j}

    Parameters
    ----------
    eigenvalues_x : ndarray
        Eigenvalues of centered feature kernel
    eigenvalues_y : ndarray
        Eigenvalues of centered spatial kernel
    threshold : float
        Filter eigenvalues below this threshold

    Returns
    -------
    composite : ndarray
        Flattened array of eigenvalue products
    """
    # Filter eigenvalues
    lambda_x = eigenvalues_x[np.abs(eigenvalues_x) > threshold]
    mu_y = eigenvalues_y[np.abs(eigenvalues_y) > threshold]

    if len(lambda_x) == 0 or len(mu_y) == 0:
        return np.array([])

    # Compute outer product and flatten
    composite = np.outer(lambda_x, mu_y).flatten()

    return composite


def hsic_pvalue(
    statistic: float,
    eigenvalues_x: NDArray,
    eigenvalues_y: NDArray,
    n: int,
    method: str = "liu",
    clamp_min: float = 1e-300,
    clamp_max: float = 1.0
) -> PValueResult:
    """Compute p-value for HSIC test statistic.

    Parameters
    ----------
    statistic : float
        HSIC test statistic
    eigenvalues_x : ndarray
        Eigenvalues of centered feature kernel
    eigenvalues_y : ndarray
        Eigenvalues of centered spatial kernel
    n : int
        Sample size
    method : str
        Method: "liu" (recommended) or "gamma"
    clamp_min : float
        Minimum p-value
    clamp_max : float
        Maximum p-value

    Returns
    -------
    result : PValueResult
    """
    if method == "liu":
        # Compute composite eigenvalues for Liu's method
        composite = compute_composite_eigenvalues(eigenvalues_x, eigenvalues_y)
        return liu_pvalue(statistic, composite, n, clamp_min, clamp_max)

    elif method == "gamma":
        return gamma_approximation_pvalue(
            statistic, eigenvalues_x, eigenvalues_y, n, clamp_min, clamp_max
        )

    else:
        raise ValueError(f"Unknown method: {method}. Use 'liu' or 'gamma'.")


def gamma_approximation_pvalue(
    statistic: float,
    eigenvalues_x: NDArray,
    eigenvalues_y: NDArray,
    n: int,
    clamp_min: float = 1e-300,
    clamp_max: float = 1.0
) -> PValueResult:
    """Compute p-value using gamma approximation (legacy method).

    The null distribution of n * HSIC follows a weighted chi-squared mixture:
        n * T_HSIC → Σᵢⱼ λᵢμⱼzᵢⱼ²

    This is approximated using a gamma distribution via moment matching.

    Parameters
    ----------
    statistic : float
        Observed HSIC statistic (unscaled)
    eigenvalues_x : ndarray
        Eigenvalues of centered feature kernel
    eigenvalues_y : ndarray
        Eigenvalues of centered spatial kernel
    n : int
        Sample size (number of spots)
    clamp_min : float
        Minimum p-value (for numerical stability)
    clamp_max : float
        Maximum p-value

    Returns
    -------
    result : PValueResult
        P-value and diagnostic information
    """
    # Filter eigenvalues
    lambda_x = eigenvalues_x[np.abs(eigenvalues_x) > 1e-10]
    mu_y = eigenvalues_y[np.abs(eigenvalues_y) > 1e-10]

    if len(lambda_x) == 0 or len(mu_y) == 0:
        return PValueResult(
            pvalue=1.0,
            statistic=statistic,
            method="gamma",
            mean=0.0,
            variance=0.0
        )

    # Compute moments of null distribution for n * T
    sum_lambda = np.sum(lambda_x)
    sum_mu = np.sum(mu_y)
    sum_lambda_sq = np.sum(lambda_x ** 2)
    sum_mu_sq = np.sum(mu_y ** 2)

    E_nT = sum_lambda * sum_mu / n
    Var_nT = 2 * sum_lambda_sq * sum_mu_sq / (n ** 2)

    if Var_nT <= 0 or E_nT <= 0:
        return PValueResult(
            pvalue=1.0,
            statistic=statistic,
            method="gamma",
            mean=E_nT,
            variance=Var_nT
        )

    # Gamma parameters via method of moments
    beta = Var_nT / E_nT
    alpha = E_nT ** 2 / Var_nT

    # Scale statistic by n for asymptotic distribution
    scaled_statistic = n * statistic

    # Compute p-value
    pvalue = stats.gamma.sf(scaled_statistic, a=alpha, scale=beta)
    pvalue = np.clip(pvalue, clamp_min, clamp_max)

    return PValueResult(
        pvalue=pvalue,
        statistic=statistic,
        method="gamma",
        mean=E_nT / n,
        variance=Var_nT / (n ** 2),
        df=alpha,  # Store gamma shape as df
        delta=beta  # Store gamma scale as delta
    )


def compute_pvalue_from_kernels(
    statistic: float,
    K_X: NDArray,
    L_Y: NDArray,
    method: str = "liu",
    eigenvalue_threshold: float = 1e-10
) -> PValueResult:
    """Compute p-value from kernel matrices.

    Convenience function that extracts eigenvalues from
    centered kernel matrices.

    Parameters
    ----------
    statistic : float
        Observed HSIC statistic
    K_X : ndarray
        Feature kernel matrix (will be centered if not already)
    L_Y : ndarray
        Spatial kernel matrix (will be centered if not already)
    method : str
        Method: "liu" (recommended) or "gamma"
    eigenvalue_threshold : float
        Filter eigenvalues below this threshold

    Returns
    -------
    result : PValueResult
    """
    n = K_X.shape[0]

    # Center kernels
    from .hsic import center_kernel
    K_X_centered = center_kernel(K_X)
    L_Y_centered = center_kernel(L_Y)

    # Extract eigenvalues
    eigenvalues_x = np.linalg.eigvalsh(K_X_centered)
    eigenvalues_y = np.linalg.eigvalsh(L_Y_centered)

    return hsic_pvalue(
        statistic, eigenvalues_x, eigenvalues_y, n, method=method
    )


def permutation_pvalue(
    K_X: NDArray,
    L_Y: NDArray,
    n_permutations: int = 1000,
    seed: int = 42
) -> Tuple[float, NDArray]:
    """Compute p-value via permutation testing.

    This is slower but provides exact p-values without distributional assumptions.

    Parameters
    ----------
    K_X : ndarray
        Feature kernel matrix
    L_Y : ndarray
        Spatial kernel matrix
    n_permutations : int
        Number of permutations
    seed : int
        Random seed

    Returns
    -------
    pvalue : float
        Permutation p-value
    null_statistics : ndarray
        HSIC statistics under permutations
    """
    from .hsic import compute_hsic_statistic, center_kernel

    rng = np.random.default_rng(seed)
    n = K_X.shape[0]

    # Center kernels once
    K_X_centered = center_kernel(K_X)
    L_Y_centered = center_kernel(L_Y)

    # Observed statistic
    T_obs = compute_hsic_statistic(K_X_centered, L_Y_centered, centered=True)

    # Permutation null distribution
    null_statistics = np.zeros(n_permutations)
    for i in range(n_permutations):
        perm = rng.permutation(n)
        K_X_perm = K_X_centered[perm][:, perm]
        T_perm = compute_hsic_statistic(K_X_perm, L_Y_centered, centered=True)
        null_statistics[i] = T_perm

    # P-value: fraction of permutations >= observed
    pvalue = (np.sum(null_statistics >= T_obs) + 1) / (n_permutations + 1)

    return pvalue, null_statistics

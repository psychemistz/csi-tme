"""Conditional HSIC for controlling spatial confounding.

When correlating isoform usage X with RBP expression Z,
spatial confounding can create spurious associations.
This module implements conditional independence testing.
"""

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import RBFInterpolator

from .hsic import compute_hsic_statistic, center_kernel, get_kernel_eigenvalues
from .pvalue import hsic_pvalue, PValueResult


@dataclass
class ConditionalHSICResult:
    """Result of conditional HSIC test.

    Attributes
    ----------
    statistic : float
        Conditional HSIC statistic
    pvalue : float
        P-value
    residuals_x : ndarray
        Residuals of X after spatial regression
    residuals_z : ndarray
        Residuals of Z after spatial regression
    """
    statistic: float
    pvalue: float
    residuals_x: Optional[NDArray] = None
    residuals_z: Optional[NDArray] = None


def spatial_residualize(
    values: NDArray,
    coords: NDArray,
    method: str = 'rbf',
    rbf_smoothing: float = 1.0
) -> NDArray:
    """Remove spatial trend from values.

    Learns spatial smoother f: Y → values and returns residuals.

    Parameters
    ----------
    values : ndarray
        Values to residualize, shape (n,) or (n, d)
    coords : ndarray
        Spatial coordinates, shape (n, 2)
    method : str
        'rbf' for RBF interpolation, 'polynomial' for polynomial
    rbf_smoothing : float
        Smoothing parameter for RBF (higher = smoother)

    Returns
    -------
    residuals : ndarray
        Residuals after removing spatial trend
    """
    if values.ndim == 1:
        values = values.reshape(-1, 1)
        squeeze = True
    else:
        squeeze = False

    # Normalize coordinates to zero mean, unit variance so that the
    # smoothing parameter has a consistent effect regardless of the
    # coordinate scale (pixel coords vs array indices, etc.)
    coord_mean = coords.mean(axis=0)
    coord_std = coords.std(axis=0)
    coord_std = np.where(coord_std > 0, coord_std, 1.0)
    coords_norm = (coords - coord_mean) / coord_std

    if method == 'rbf':
        # RBF interpolation as spatial smoother
        interpolator = RBFInterpolator(
            coords_norm, values,
            smoothing=rbf_smoothing,
            kernel='thin_plate_spline'
        )
        predicted = interpolator(coords_norm)

    elif method == 'polynomial':
        # Polynomial regression (degree 2)
        x, y = coords_norm[:, 0], coords_norm[:, 1]
        # Design matrix: [1, x, y, x^2, xy, y^2]
        X_design = np.column_stack([
            np.ones_like(x),
            x, y,
            x ** 2, x * y, y ** 2
        ])
        # Solve least squares for each column of values
        coeffs, _, _, _ = np.linalg.lstsq(X_design, values, rcond=None)
        predicted = X_design @ coeffs

    else:
        raise ValueError(f"Unknown method: {method}")

    residuals = values - predicted

    if squeeze:
        residuals = residuals.ravel()

    return residuals


def conditional_hsic(
    X: NDArray,
    Z: NDArray,
    coords: NDArray,
    method: str = 'rbf',
    rbf_smoothing: float = 1.0
) -> ConditionalHSICResult:
    """Test independence of X and Z conditional on spatial location.

    This controls for spatial confounding when testing whether
    isoform usage (X) correlates with RBP expression (Z).

    Algorithm:
    1. Learn spatial smoothers: f_X: Y → X, f_Z: Y → Z
    2. Compute residuals: ε_X = X - f_X(Y), ε_Z = Z - f_Z(Y)
    3. Test: HSIC(ε_X, ε_Z)

    Parameters
    ----------
    X : ndarray
        First variable (e.g., isoform proportions), shape (n,) or (n, d)
    Z : ndarray
        Second variable (e.g., RBP expression), shape (n,) or (n, k)
    coords : ndarray
        Spatial coordinates, shape (n, 2)
    method : str
        Spatial smoothing method ('rbf' or 'polynomial')
    rbf_smoothing : float
        Smoothing parameter for RBF method

    Returns
    -------
    result : ConditionalHSICResult
    """
    n = coords.shape[0]

    # Residualize both variables
    eps_X = spatial_residualize(X, coords, method, rbf_smoothing)
    eps_Z = spatial_residualize(Z, coords, method, rbf_smoothing)

    # Build kernels on residuals
    if eps_X.ndim == 1:
        eps_X = eps_X.reshape(-1, 1)
    if eps_Z.ndim == 1:
        eps_Z = eps_Z.reshape(-1, 1)

    K_X = eps_X @ eps_X.T
    K_Z = eps_Z @ eps_Z.T

    # Center kernels
    K_X_centered = center_kernel(K_X)
    K_Z_centered = center_kernel(K_Z)

    # Compute HSIC
    statistic = compute_hsic_statistic(K_X_centered, K_Z_centered, centered=True)

    # Get eigenvalues for p-value
    eigenvalues_x = get_kernel_eigenvalues(K_X_centered)
    eigenvalues_z = get_kernel_eigenvalues(K_Z_centered)

    # Compute p-value using Liu's method (default)
    pvalue_result = hsic_pvalue(
        statistic, eigenvalues_x, eigenvalues_z, n, method='liu'
    )

    return ConditionalHSICResult(
        statistic=statistic,
        pvalue=pvalue_result.pvalue,
        residuals_x=eps_X.squeeze(),
        residuals_z=eps_Z.squeeze()
    )


def partial_correlation_spatial(
    X: NDArray,
    Z: NDArray,
    coords: NDArray,
    method: str = 'rbf'
) -> Tuple[float, float]:
    """Compute partial correlation controlling for space.

    Simpler alternative to conditional HSIC for scalar variables.

    Parameters
    ----------
    X : ndarray
        First variable, shape (n,)
    Z : ndarray
        Second variable, shape (n,)
    coords : ndarray
        Spatial coordinates, shape (n, 2)
    method : str
        Spatial smoothing method

    Returns
    -------
    correlation : float
        Partial correlation coefficient
    pvalue : float
        P-value (Fisher's z-test)
    """
    from scipy import stats as scipy_stats

    n = len(X)

    # Residualize
    eps_X = spatial_residualize(X, coords, method)
    eps_Z = spatial_residualize(Z, coords, method)

    # Pearson correlation on residuals
    correlation, pvalue = scipy_stats.pearsonr(eps_X, eps_Z)

    return correlation, pvalue

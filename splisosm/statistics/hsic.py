"""HSIC (Hilbert-Schmidt Independence Criterion) test statistics.

Implements three variants matching the official SPLISOSM:
- HSIC-GC: Tests total gene expression spatial variation
- HSIC-IR: Tests isoform ratio spatial variation (independent of expression)
- HSIC-IC: Tests combined expression + composition spatial variation

Uses unbiased (n-1)² normalization following the official implementation.

Supports GPU acceleration via CuPy when use_gpu=True.
"""

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..kernels.compositional import clr_transform
from ..utils.backend import (
    get_backend,
    to_numpy,
    gpu_center_kernel,
    gpu_trace_product,
    gpu_kernel_compute,
    gpu_eigvalsh,
    CUPY_AVAILABLE,
)


@dataclass
class HSICResult:
    """Result of an HSIC test.

    Attributes
    ----------
    statistic : float
        HSIC test statistic value
    pvalue : float
        P-value
    eigenvalues_x : ndarray
        Eigenvalues of feature kernel (for diagnostics)
    eigenvalues_y : ndarray
        Eigenvalues of spatial kernel (for diagnostics)
    """
    statistic: float
    pvalue: float
    eigenvalues_x: Optional[NDArray] = None
    eigenvalues_y: Optional[NDArray] = None


def center_kernel(K: NDArray, use_gpu: bool = False) -> NDArray:
    """Center kernel matrix: HKH.

    Parameters
    ----------
    K : ndarray
        Kernel matrix of shape (n, n)
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    K_centered : ndarray
        Centered kernel matrix
    """
    if use_gpu and CUPY_AVAILABLE:
        return gpu_center_kernel(K, use_gpu=True)

    n = K.shape[0]
    row_means = K.mean(axis=1, keepdims=True)
    col_means = K.mean(axis=0, keepdims=True)
    grand_mean = K.mean()

    return K - row_means - col_means + grand_mean


def compute_hsic_statistic(
    K_X: NDArray,
    L_Y: NDArray,
    centered: bool = False,
    unbiased: bool = True,
    use_gpu: bool = False
) -> float:
    """Compute HSIC test statistic.

    For unbiased=True (default, matches official SPLISOSM):
        T_HSIC = tr(K_X H L_Y H) / (n-1)²

    For unbiased=False (legacy):
        T_HSIC = tr(K_X H L_Y H) / n

    Parameters
    ----------
    K_X : ndarray
        Feature kernel matrix (n, n)
    L_Y : ndarray
        Spatial kernel matrix (n, n)
    centered : bool
        If True, kernels are already centered
    unbiased : bool
        If True, use (n-1)² normalization (matches official SPLISOSM)
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    T : float
        HSIC statistic
    """
    n = K_X.shape[0]

    if not centered:
        K_X = center_kernel(K_X, use_gpu=use_gpu)
        L_Y = center_kernel(L_Y, use_gpu=use_gpu)

    # T = tr(K_X @ L_Y) / normalization
    # Use trace(A @ B) = sum(A * B.T) for efficiency
    if use_gpu and CUPY_AVAILABLE:
        trace_val = gpu_trace_product(K_X, L_Y, use_gpu=True)
    else:
        trace_val = np.sum(K_X * L_Y)

    if unbiased:
        # Unbiased estimator using (n-1)² - matches official SPLISOSM
        T = trace_val / ((n - 1) ** 2)
    else:
        # Biased V-statistic using n
        T = trace_val / n

    return float(T)


def compute_hsic_linear(
    X: NDArray,
    Y: NDArray,
    L_sp: NDArray,
    center: bool = True
) -> Tuple[float, NDArray, NDArray]:
    """Compute linear HSIC statistic efficiently.

    This matches the official SPLISOSM linear_hsic_test formulation:
        hsic_scaled = ||Y^T @ X||_F² / (n-1)²

    Parameters
    ----------
    X : ndarray
        Feature matrix of shape (n, p) - e.g., centered isoform ratios
    Y : ndarray
        Spatial kernel features of shape (n, k) or kernel (n, n)
    L_sp : ndarray
        Centered spatial kernel of shape (n, n)
    center : bool
        If True, center X before computation

    Returns
    -------
    T : float
        HSIC statistic
    eigenvalues_x : ndarray
        Eigenvalues of X^T X
    eigenvalues_y : ndarray
        Eigenvalues of centered spatial kernel
    """
    n = X.shape[0]

    if center:
        X = X - X.mean(axis=0, keepdims=True)

    # Compute statistic as ||X^T @ L_sp @ X||_F normalized
    # This is equivalent to tr(K_X @ H @ L_sp @ H) where K_X = X @ X^T
    XtL = X.T @ L_sp
    hsic_trace = np.sum(XtL * X.T)  # tr(X^T L X) = sum((X^T L) * X^T)

    T = hsic_trace / ((n - 1) ** 2)

    # Eigenvalues for p-value computation
    eigenvalues_x = np.linalg.eigvalsh(X.T @ X)
    eigenvalues_y = np.linalg.eigvalsh(L_sp)

    return T, eigenvalues_x, eigenvalues_y


def compute_hsic_gc(
    expression: NDArray,
    spatial_kernel: NDArray,
    centered: bool = False,
    unbiased: bool = True,
    use_gpu: bool = False
) -> Tuple[float, NDArray]:
    """Compute HSIC-GC for gene expression spatial variation.

    Tests whether total gene expression varies spatially.
    K_g = g @ g^T (linear kernel on expression)

    Parameters
    ----------
    expression : ndarray
        Gene expression values of shape (n_spots,)
    spatial_kernel : ndarray
        Spatial kernel matrix of shape (n_spots, n_spots)
    centered : bool
        If True, spatial kernel is already centered
    unbiased : bool
        If True, use (n-1)² normalization
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    T : float
        HSIC-GC statistic
    K_g : ndarray
        Gene expression kernel (for p-value computation)
    """
    if use_gpu and CUPY_AVAILABLE:
        xp = get_backend(use_gpu=True)
        expression = xp.asarray(expression)
        g = expression - expression.mean()
        K_g = xp.outer(g, g)
        K_g_np = to_numpy(K_g)
    else:
        # Center expression
        g = expression - expression.mean()
        # Expression kernel: K_g = g @ g^T
        K_g = np.outer(g, g)
        K_g_np = K_g

    T = compute_hsic_statistic(K_g_np, spatial_kernel, centered=centered,
                                unbiased=unbiased, use_gpu=use_gpu)

    return T, K_g_np


def compute_hsic_ir(
    isoform_proportions: NDArray,
    spatial_kernel: NDArray,
    pseudocount: float = 1e-10,
    centered: bool = False,
    unbiased: bool = True,
    transform: str = "none",
    use_gpu: bool = False
) -> Tuple[float, NDArray]:
    """Compute HSIC-IR for isoform ratio spatial variation.

    Tests whether isoform ratios vary spatially, independent of
    total gene expression. This is the primary statistic for
    detecting microenvironment-driven isoform switching.

    Parameters
    ----------
    isoform_proportions : ndarray
        Isoform proportions of shape (n_spots, n_isoforms)
    spatial_kernel : ndarray
        Spatial kernel matrix
    pseudocount : float
        Added to proportions before log-ratio transform
    centered : bool
        If True, spatial kernel is already centered
    unbiased : bool
        If True, use (n-1)² normalization
    transform : str
        Compositional transform: "clr", "ilr", "alr", or "none"
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    T : float
        HSIC-IR statistic
    K_X : ndarray
        Compositional kernel (for p-value computation)
    """
    from ..kernels.compositional import clr_transform, ilr_transform, alr_transform

    # Apply compositional transform (CPU - transforms use numpy)
    if transform == "clr":
        X = clr_transform(isoform_proportions, pseudocount)
    elif transform == "ilr":
        X = ilr_transform(isoform_proportions, pseudocount)
    elif transform == "alr":
        X = alr_transform(isoform_proportions, pseudocount)
    elif transform == "none":
        # Simple centering without log-ratio
        X = isoform_proportions - isoform_proportions.mean(axis=0, keepdims=True)
    else:
        raise ValueError(f"Unknown transform: {transform}")

    # Compositional kernel (GPU-accelerated if requested)
    if use_gpu and CUPY_AVAILABLE:
        K_X = gpu_kernel_compute(X, use_gpu=True)
    else:
        K_X = X @ X.T

    T = compute_hsic_statistic(K_X, spatial_kernel, centered=centered,
                                unbiased=unbiased, use_gpu=use_gpu)

    return T, K_X


def compute_hsic_ic(
    isoform_counts: NDArray,
    spatial_kernel: NDArray,
    centered: bool = False,
    unbiased: bool = True,
    use_gpu: bool = False
) -> Tuple[float, NDArray]:
    """Compute HSIC-IC for isoform counts spatial variation.

    Tests whether raw isoform counts vary spatially.
    This is the combined measure that captures both expression
    and composition changes.

    Parameters
    ----------
    isoform_counts : ndarray
        Isoform counts of shape (n_spots, n_isoforms)
    spatial_kernel : ndarray
        Spatial kernel matrix
    centered : bool
        If True, spatial kernel is already centered
    unbiased : bool
        If True, use (n-1)² normalization
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    T : float
        HSIC-IC statistic
    K_X : ndarray
        Isoform count kernel (for p-value computation)
    """
    # Center counts column-wise
    X = isoform_counts - isoform_counts.mean(axis=0, keepdims=True)

    # Count kernel (GPU-accelerated if requested)
    if use_gpu and CUPY_AVAILABLE:
        K_X = gpu_kernel_compute(X, use_gpu=True)
    else:
        K_X = X @ X.T

    T = compute_hsic_statistic(K_X, spatial_kernel, centered=centered,
                                unbiased=unbiased, use_gpu=use_gpu)

    return T, K_X


def compute_hsic_combined(
    isoform_proportions: NDArray,
    expression: NDArray,
    spatial_kernel: NDArray,
    alpha: float = 0.5,
    pseudocount: float = 1e-10,
    centered: bool = False,
    unbiased: bool = True,
    use_gpu: bool = False
) -> Tuple[float, NDArray]:
    """Compute combined expression + composition HSIC.

    Tests for spatial variation in both total expression and isoform
    composition using a combined kernel.

    K_XG = K_X ⊙ K_G + α * (K_X + K_G)

    Parameters
    ----------
    isoform_proportions : ndarray
        Isoform proportions of shape (n_spots, n_isoforms)
    expression : ndarray
        Total gene expression of shape (n_spots,)
    spatial_kernel : ndarray
        Spatial kernel matrix
    alpha : float
        Weight for additive component
    pseudocount : float
        Added to proportions before CLR transform
    centered : bool
        If True, spatial kernel is already centered
    unbiased : bool
        If True, use (n-1)² normalization
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    T : float
        Combined HSIC statistic
    K_XG : ndarray
        Combined kernel (for p-value computation)
    """
    # CLR transform (CPU)
    X_clr = clr_transform(isoform_proportions, pseudocount)

    if use_gpu and CUPY_AVAILABLE:
        xp = get_backend(use_gpu=True)
        X_clr_gpu = xp.asarray(X_clr)
        K_X = X_clr_gpu @ X_clr_gpu.T

        expression_gpu = xp.asarray(expression)
        g = expression_gpu - expression_gpu.mean()
        g = g / (g.std() + 1e-10)
        K_G = xp.outer(g, g)

        # Combined kernel
        K_XG = K_X * K_G + alpha * (K_X + K_G)
        K_XG = to_numpy(K_XG)
    else:
        K_X = X_clr @ X_clr.T

        # Expression kernel
        g = expression - expression.mean()
        g = g / (g.std() + 1e-10)
        K_G = np.outer(g, g)

        # Combined kernel
        K_XG = K_X * K_G + alpha * (K_X + K_G)

    T = compute_hsic_statistic(K_XG, spatial_kernel, centered=centered,
                                unbiased=unbiased, use_gpu=use_gpu)

    return T, K_XG


def get_kernel_eigenvalues(K: NDArray, threshold: float = 1e-10,
                           use_gpu: bool = False) -> NDArray:
    """Extract eigenvalues from kernel matrix.

    Parameters
    ----------
    K : ndarray
        Kernel matrix (should be symmetric)
    threshold : float
        Filter eigenvalues with |λ| < threshold
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    eigenvalues : ndarray
        Filtered eigenvalues
    """
    if use_gpu and CUPY_AVAILABLE:
        eigenvalues = gpu_eigvalsh(K, use_gpu=True)
    else:
        eigenvalues = np.linalg.eigvalsh(K)

    eigenvalues = eigenvalues[np.abs(eigenvalues) > threshold]
    return eigenvalues

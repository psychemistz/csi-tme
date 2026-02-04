"""Compositional data kernels and transforms.

Implements log-ratio transforms for compositional data (e.g., isoform proportions):
- CLR (Centered Log-Ratio): Maps to D-dimensional space, rows sum to zero
- ILR (Isometric Log-Ratio): Maps to (D-1)-dimensional space, preserves distances
- ALR (Additive Log-Ratio): Maps to (D-1)-dimensional space, uses reference component

These transforms respect the Aitchison geometry of simplex-constrained data.
"""

from typing import Optional

import numpy as np
from numpy.typing import NDArray

from .base import Kernel


def clr_transform(proportions: NDArray, pseudocount: float = 1e-10) -> NDArray:
    """Apply Centered Log-Ratio (CLR) transform to compositional data.

    CLR respects the Aitchison geometry of simplex-constrained data.
    For proportions p = (p₁, ..., p_q):
        clr(p)_j = log(p_j / geometric_mean(p))

    Parameters
    ----------
    proportions : ndarray
        Compositional data of shape (n_samples, n_components).
        Each row should sum to 1 (or close to it).
    pseudocount : float
        Small value added before log transform to handle zeros.

    Returns
    -------
    clr_data : ndarray
        CLR-transformed data of shape (n_samples, n_components).
        Rows sum to zero.
    """
    # Add pseudocount for numerical stability
    p = proportions + pseudocount

    # Normalize rows (in case they don't sum to 1)
    p = p / p.sum(axis=1, keepdims=True)

    # Geometric mean per row
    log_p = np.log(p)
    log_gm = log_p.mean(axis=1, keepdims=True)

    # CLR transform
    clr_data = log_p - log_gm

    return clr_data


def ilr_transform(proportions: NDArray, pseudocount: float = 1e-10) -> NDArray:
    """Apply Isometric Log-Ratio (ILR) transform to compositional data.

    ILR maps from D-dimensional simplex to (D-1)-dimensional Euclidean space.
    Uses the Helmert sub-composition approach.

    Parameters
    ----------
    proportions : ndarray
        Compositional data of shape (n_samples, n_components)
    pseudocount : float
        Small value added before log transform

    Returns
    -------
    ilr_data : ndarray
        ILR-transformed data of shape (n_samples, n_components-1)
    """
    n, d = proportions.shape
    if d < 2:
        raise ValueError("ILR requires at least 2 components")

    # CLR first
    clr_data = clr_transform(proportions, pseudocount)

    # Helmert contrast matrix (orthonormal basis for simplex)
    # Creates d-1 orthonormal contrasts
    V = np.zeros((d, d - 1))
    for i in range(d - 1):
        V[:i + 1, i] = 1.0 / (i + 1)
        V[i + 1, i] = -1.0
        V[:, i] *= np.sqrt((i + 1) / (i + 2))

    ilr_data = clr_data @ V

    return ilr_data


def alr_transform(
    proportions: NDArray,
    pseudocount: float = 1e-10,
    reference: int = -1
) -> NDArray:
    """Apply Additive Log-Ratio (ALR) transform to compositional data.

    ALR maps from D-dimensional simplex to (D-1)-dimensional space
    using one component as reference:
        alr(p)_j = log(p_j / p_ref)

    Parameters
    ----------
    proportions : ndarray
        Compositional data of shape (n_samples, n_components)
    pseudocount : float
        Small value added before log transform
    reference : int
        Index of reference component (default: last component)

    Returns
    -------
    alr_data : ndarray
        ALR-transformed data of shape (n_samples, n_components-1)
    """
    n, d = proportions.shape
    if d < 2:
        raise ValueError("ALR requires at least 2 components")

    # Add pseudocount
    p = proportions + pseudocount
    p = p / p.sum(axis=1, keepdims=True)

    # Reference component
    p_ref = p[:, reference:reference+1]

    # Compute ALR (exclude reference column)
    indices = [i for i in range(d) if i != (reference % d)]
    alr_data = np.log(p[:, indices] / p_ref)

    return alr_data


def radial_transform(proportions: NDArray, pseudocount: float = 1e-10) -> NDArray:
    """Apply radial (unit-norm) transform to compositional data.

    Simply normalizes each row to unit norm without log transformation.
    This is the simplest transform that handles the simplex constraint.

    Parameters
    ----------
    proportions : ndarray
        Compositional data of shape (n_samples, n_components)
    pseudocount : float
        Small value added for numerical stability

    Returns
    -------
    radial_data : ndarray
        Unit-normalized data of shape (n_samples, n_components)
    """
    p = proportions + pseudocount
    norms = np.linalg.norm(p, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)  # Avoid division by zero
    return p / norms


class CompositionalKernel(Kernel):
    """Kernel for compositional (simplex-constrained) data.

    Uses CLR transform followed by linear kernel:
        K_X = clr(X) @ clr(X)^T

    This is appropriate for isoform ratio analysis where
    the proportions must sum to 1.

    Parameters
    ----------
    pseudocount : float
        Added to proportions before log transform
    use_ilr : bool
        If True, use ILR instead of CLR (reduces dimensionality)
    """

    def __init__(self, pseudocount: float = 1e-10, use_ilr: bool = False):
        super().__init__()
        self.pseudocount = pseudocount
        self.use_ilr = use_ilr
        self._transformed_data: Optional[NDArray] = None

    def compute(self, proportions: NDArray, **kwargs) -> NDArray:
        """Compute compositional kernel.

        Parameters
        ----------
        proportions : ndarray
            Isoform proportions of shape (n_spots, n_isoforms).
            Each row should sum to 1.

        Returns
        -------
        K : ndarray
            Kernel matrix of shape (n_spots, n_spots)
        """
        if self.use_ilr:
            X = ilr_transform(proportions, self.pseudocount)
        else:
            X = clr_transform(proportions, self.pseudocount)

        self._transformed_data = X
        self._kernel_matrix = X @ X.T
        self._centered = False

        return self._kernel_matrix

    @property
    def transformed_data(self) -> Optional[NDArray]:
        """Return the CLR/ILR transformed data."""
        return self._transformed_data


class ExpressionCompositionalKernel(Kernel):
    """Combined kernel for expression + composition (HSIC-IC).

    Combines total gene expression with isoform composition:
        K_XG = K_X ⊙ K_G + α * (K_X + K_G)

    where K_X is compositional kernel and K_G is expression kernel.

    Parameters
    ----------
    alpha : float
        Weighting for additive term (default 0.5)
    pseudocount : float
        Added to proportions before log transform
    """

    def __init__(self, alpha: float = 0.5, pseudocount: float = 1e-10):
        super().__init__()
        self.alpha = alpha
        self.pseudocount = pseudocount

    def compute(self, proportions: NDArray, expression: NDArray,
                **kwargs) -> NDArray:
        """Compute combined expression-composition kernel.

        Parameters
        ----------
        proportions : ndarray
            Isoform proportions of shape (n_spots, n_isoforms)
        expression : ndarray
            Total gene expression of shape (n_spots,)

        Returns
        -------
        K : ndarray
            Combined kernel matrix
        """
        # Compositional kernel
        X_clr = clr_transform(proportions, self.pseudocount)
        K_X = X_clr @ X_clr.T

        # Expression kernel (outer product of standardized expression)
        g = (expression - expression.mean()) / (expression.std() + 1e-10)
        K_G = np.outer(g, g)

        # Combined kernel
        K = K_X * K_G + self.alpha * (K_X + K_G)

        self._kernel_matrix = K
        self._centered = False

        return K

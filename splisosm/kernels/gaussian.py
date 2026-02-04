"""Gaussian (RBF) kernel implementation."""

from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.spatial.distance import cdist

from .base import Kernel


class GaussianKernel(Kernel):
    """Gaussian (Radial Basis Function) kernel.

    K(x, y) = exp(-||x - y||² / (2σ²))

    Parameters
    ----------
    sigma : float or None
        Bandwidth parameter. If None, uses median heuristic.
    """

    def __init__(self, sigma: Optional[float] = None):
        super().__init__()
        self.sigma = sigma
        self._computed_sigma: Optional[float] = None

    def compute(self, X: NDArray, **kwargs) -> NDArray:
        """Compute Gaussian kernel matrix.

        Parameters
        ----------
        X : ndarray
            Data matrix of shape (n_samples, n_features)

        Returns
        -------
        K : ndarray
            Kernel matrix of shape (n_samples, n_samples)
        """
        # Compute pairwise squared distances
        dists_sq = cdist(X, X, metric='sqeuclidean')

        # Determine bandwidth
        if self.sigma is not None:
            sigma = self.sigma
        else:
            # Median heuristic: sigma = median of pairwise distances
            # Use only upper triangle to avoid zeros on diagonal
            upper_tri = dists_sq[np.triu_indices_from(dists_sq, k=1)]
            median_dist = np.median(np.sqrt(upper_tri))
            sigma = median_dist if median_dist > 0 else 1.0

        self._computed_sigma = sigma

        # Compute kernel
        K = np.exp(-dists_sq / (2 * sigma ** 2))

        self._kernel_matrix = K
        self._centered = False
        return K

    @property
    def bandwidth(self) -> Optional[float]:
        """Return the computed or specified bandwidth."""
        return self._computed_sigma if self._computed_sigma else self.sigma

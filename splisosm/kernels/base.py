"""Abstract base class for kernel implementations."""

from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy import sparse


class Kernel(ABC):
    """Abstract base class for spatial and feature kernels.

    All kernels should implement:
    - compute(): Build the kernel matrix from input data
    - center(): Apply centering transformation H @ K @ H
    """

    def __init__(self):
        self._kernel_matrix: Optional[NDArray] = None
        self._centered: bool = False

    @abstractmethod
    def compute(self, X: NDArray, **kwargs) -> NDArray:
        """Compute the kernel matrix from input data.

        Parameters
        ----------
        X : ndarray
            Input data, shape depends on kernel type
        **kwargs
            Additional kernel-specific parameters

        Returns
        -------
        K : ndarray
            Kernel matrix of shape (n, n)
        """
        pass

    @property
    def matrix(self) -> Optional[NDArray]:
        """Return the computed kernel matrix."""
        return self._kernel_matrix

    @staticmethod
    def centering_matrix(n: int, sparse_format: bool = False) -> NDArray:
        """Compute the centering matrix H = I - (1/n) * 11^T.

        Parameters
        ----------
        n : int
            Matrix dimension
        sparse_format : bool
            If True, return sparse matrix

        Returns
        -------
        H : ndarray or sparse matrix
            Centering matrix of shape (n, n)
        """
        if sparse_format:
            I = sparse.eye(n, format='csr')
            ones = sparse.csr_matrix(np.ones((n, n)) / n)
            return I - ones
        else:
            return np.eye(n) - np.ones((n, n)) / n

    def center(self, K: Optional[NDArray] = None) -> NDArray:
        """Apply centering: H @ K @ H.

        Centers the kernel matrix to have zero mean in both rows and columns.
        This is required for HSIC computation.

        Parameters
        ----------
        K : ndarray, optional
            Kernel matrix to center. If None, uses self._kernel_matrix

        Returns
        -------
        K_centered : ndarray
            Centered kernel matrix
        """
        if K is None:
            K = self._kernel_matrix
        if K is None:
            raise ValueError("No kernel matrix available. Call compute() first.")

        n = K.shape[0]
        H = self.centering_matrix(n)

        # HKH = H @ K @ H (more efficient formulation)
        # = K - K @ 1 @ 1^T/n - 1 @ 1^T @ K/n + 1 @ 1^T @ K @ 1 @ 1^T/n^2
        row_means = K.mean(axis=1, keepdims=True)
        col_means = K.mean(axis=0, keepdims=True)
        grand_mean = K.mean()

        K_centered = K - row_means - col_means + grand_mean
        self._centered = True

        return K_centered

    def get_eigenvalues(self, K: Optional[NDArray] = None,
                        threshold: float = 1e-10) -> NDArray:
        """Compute eigenvalues of kernel matrix.

        Parameters
        ----------
        K : ndarray, optional
            Kernel matrix. If None, uses self._kernel_matrix
        threshold : float
            Filter out eigenvalues with |λ| < threshold

        Returns
        -------
        eigenvalues : ndarray
            Sorted eigenvalues (ascending)
        """
        if K is None:
            K = self._kernel_matrix
        if K is None:
            raise ValueError("No kernel matrix available. Call compute() first.")

        # Use eigvalsh for symmetric matrices (more stable)
        eigenvalues = np.linalg.eigvalsh(K)

        # Filter small eigenvalues
        eigenvalues = eigenvalues[np.abs(eigenvalues) > threshold]

        return eigenvalues


class LinearKernel(Kernel):
    """Linear kernel: K(X, Y) = X @ Y^T."""

    def compute(self, X: NDArray, **kwargs) -> NDArray:
        """Compute linear kernel.

        Parameters
        ----------
        X : ndarray
            Data matrix of shape (n_samples, n_features)

        Returns
        -------
        K : ndarray
            Kernel matrix of shape (n_samples, n_samples)
        """
        self._kernel_matrix = X @ X.T
        self._centered = False
        return self._kernel_matrix

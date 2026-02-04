"""ICAR (Intrinsic Conditional Autoregressive) spatial kernel.

Implements the Gaussian Markov Random Field kernel:
    L = (D - ρW)^{-1}

where W is the adjacency matrix, D is the degree matrix,
and ρ is the autocorrelation parameter.

Also provides low-rank GFT (Graph Fourier Transform) approximation
for efficient computation on large datasets.

Supports GPU acceleration via CuPy when use_gpu=True.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu, eigsh
from scipy.spatial import KDTree

from .base import Kernel
from ..utils.backend import (
    get_backend,
    to_numpy,
    gpu_eigh,
    CUPY_AVAILABLE,
)


class ICARKernel(Kernel):
    """ICAR spatial kernel using k-NN graph structure.

    The ICAR kernel respects tissue topology through a graph-based
    spatial correlation structure. Low eigenvalues capture smooth
    large-scale patterns; high eigenvalues capture local fluctuations.

    Parameters
    ----------
    k_neighbors : int
        Number of nearest neighbors for graph construction.
        Typical values: 20-50 for neural tissue, 10-15 for tumor core.
    rho : float
        Autocorrelation strength, must be in (0, 1).
        Higher values = stronger spatial smoothing.
        Typical values: 0.95-0.99
    regularization : float
        Small value added to diagonal for numerical stability.
    """

    def __init__(self, k_neighbors: int = 20, rho: float = 0.99,
                 regularization: float = 1e-6):
        super().__init__()
        if not 0 < rho < 1:
            raise ValueError(f"rho must be in (0, 1), got {rho}")
        self.k_neighbors = k_neighbors
        self.rho = rho
        self.regularization = regularization
        self._adjacency: Optional[csr_matrix] = None
        self._precision: Optional[csr_matrix] = None

    def build_adjacency(self, coords: NDArray) -> csr_matrix:
        """Build symmetrized k-NN adjacency matrix.

        Parameters
        ----------
        coords : ndarray
            Spatial coordinates of shape (n_spots, 2)

        Returns
        -------
        W : csr_matrix
            Symmetric adjacency matrix
        """
        n = coords.shape[0]
        k = min(self.k_neighbors, n - 1)

        # Build KDTree for efficient neighbor search
        tree = KDTree(coords)

        # Query k+1 neighbors (includes self)
        distances, indices = tree.query(coords, k=k + 1)

        # Build adjacency matrix (excluding self-connections)
        row_indices = []
        col_indices = []

        for i in range(n):
            for j_idx in range(1, k + 1):  # Skip self (index 0)
                j = indices[i, j_idx]
                row_indices.append(i)
                col_indices.append(j)

        # Create sparse matrix with ones
        data = np.ones(len(row_indices))
        W = csr_matrix((data, (row_indices, col_indices)), shape=(n, n))

        # Symmetrize: W = (W + W^T) / 2, then threshold to binary
        W = W + W.T
        W.data = np.ones_like(W.data)  # Binary adjacency

        self._adjacency = W
        return W

    def build_precision(self, W: Optional[csr_matrix] = None) -> csr_matrix:
        """Build ICAR precision matrix Q = D - ρW.

        Parameters
        ----------
        W : csr_matrix, optional
            Adjacency matrix. If None, uses stored adjacency.

        Returns
        -------
        Q : csr_matrix
            Precision matrix
        """
        if W is None:
            W = self._adjacency
        if W is None:
            raise ValueError("No adjacency matrix. Call build_adjacency first.")

        n = W.shape[0]

        # Degree matrix (diagonal)
        degrees = np.array(W.sum(axis=1)).flatten()
        D = sparse.diags(degrees, format='csr')

        # Precision matrix: Q = D - ρW + εI (regularization for stability)
        Q = D - self.rho * W
        Q = Q + sparse.eye(n, format='csr') * self.regularization

        self._precision = Q
        return Q

    def compute(self, coords: NDArray, **kwargs) -> NDArray:
        """Compute ICAR kernel matrix L = Q^{-1}.

        Parameters
        ----------
        coords : ndarray
            Spatial coordinates of shape (n_spots, 2)

        Returns
        -------
        L : ndarray
            ICAR kernel matrix (covariance matrix)
        """
        n = coords.shape[0]

        # Build adjacency and precision
        W = self.build_adjacency(coords)
        Q = self.build_precision(W)

        # Compute inverse via sparse LU factorization
        # For moderate n, we compute full inverse
        # For large n (>5000), consider using solve instead
        if n <= 5000:
            # Sparse LU factorization
            lu = splu(Q.tocsc())
            # Compute inverse column by column
            L = np.zeros((n, n))
            I = np.eye(n)
            for i in range(n):
                L[:, i] = lu.solve(I[:, i])
        else:
            # For very large matrices, use sparse solve
            # This returns dense matrix anyway, so we compute columnwise
            lu = splu(Q.tocsc())
            L = np.zeros((n, n))
            I = np.eye(n)
            for i in range(n):
                L[:, i] = lu.solve(I[:, i])

        # Ensure symmetry (numerical errors can break this)
        L = (L + L.T) / 2

        self._kernel_matrix = L
        self._centered = False
        return L

    def get_graph_fourier_eigenvalues(self, use_gpu: bool = False) -> Tuple[NDArray, NDArray]:
        """Compute eigendecomposition for Graph Fourier analysis.

        Parameters
        ----------
        use_gpu : bool
            If True, use GPU acceleration via CuPy

        Returns
        -------
        eigenvalues : ndarray
            Eigenvalues sorted ascending (low = smooth patterns)
        eigenvectors : ndarray
            Corresponding eigenvectors (columns)
        """
        if self._kernel_matrix is None:
            raise ValueError("Kernel not computed. Call compute() first.")

        if use_gpu and CUPY_AVAILABLE:
            eigenvalues, eigenvectors = gpu_eigh(self._kernel_matrix, use_gpu=True)
        else:
            eigenvalues, eigenvectors = np.linalg.eigh(self._kernel_matrix)

        # Sort ascending
        idx = np.argsort(eigenvalues)
        return eigenvalues[idx], eigenvectors[:, idx]

    def apply_boundary_weighting(self, boundary_scores: NDArray) -> NDArray:
        """Modify kernel for boundary-weighted analysis.

        Parameters
        ----------
        boundary_scores : ndarray
            Score for each spot indicating proximity to boundary.
            Shape (n_spots,), values typically in [0.1, 1].

        Returns
        -------
        L_boundary : ndarray
            Boundary-weighted kernel matrix
        """
        if self._kernel_matrix is None:
            raise ValueError("Kernel not computed. Call compute() first.")

        L = self._kernel_matrix
        n = L.shape[0]

        # Build boundary weight matrix: B_ij = sqrt(b_i * b_j)
        B = np.sqrt(np.outer(boundary_scores, boundary_scores))

        # Element-wise multiplication: L_boundary = L ⊙ B
        L_boundary = L * B

        return L_boundary

    def compute_low_rank(self, coords: NDArray, n_components: int = 50,
                         **kwargs) -> Tuple[NDArray, NDArray, NDArray]:
        """Compute low-rank GFT approximation of ICAR kernel.

        Uses the Graph Fourier Transform to approximate L ≈ U_k @ Λ_k @ U_k.T
        where U_k contains the top-k eigenvectors (by eigenvalue magnitude).

        CAUTION (Theorem 1): Low-rank approximations sacrifice statistical power.
        High-frequency eigenvectors capture biologically meaningful local patterns
        (e.g., tumor boundaries). Use full-rank ICAR when possible.

        Parameters
        ----------
        coords : ndarray
            Spatial coordinates of shape (n_spots, 2)
        n_components : int
            Number of eigenvectors to retain (default: 50)

        Returns
        -------
        L_approx : ndarray
            Low-rank approximation of kernel matrix
        eigenvalues : ndarray
            Retained eigenvalues
        eigenvectors : ndarray
            Retained eigenvectors
        """
        n = coords.shape[0]
        n_components = min(n_components, n - 1)

        # Build adjacency and precision
        W = self.build_adjacency(coords)
        Q = self.build_precision(W)

        # Compute top-k eigenvalues/vectors of precision matrix
        # For ICAR, we want largest eigenvalues of L = Q^{-1}
        # which are smallest eigenvalues of Q
        # Use shift-invert mode: solve (Q - sigma*I)^{-1} @ x = theta @ x
        try:
            eigenvalues_q, eigenvectors = eigsh(
                Q.tocsc(), k=n_components, which='SM',  # Smallest magnitude
                return_eigenvectors=True
            )
        except Exception:
            # Fallback to dense computation if sparse fails
            Q_dense = Q.toarray()
            eigenvalues_q, eigenvectors_full = np.linalg.eigh(Q_dense)
            idx = np.argsort(np.abs(eigenvalues_q))[:n_components]
            eigenvalues_q = eigenvalues_q[idx]
            eigenvectors = eigenvectors_full[:, idx]

        # Eigenvalues of L = 1/eigenvalues of Q
        eigenvalues_l = 1.0 / eigenvalues_q

        # Low-rank approximation: L ≈ U @ diag(λ) @ U.T
        L_approx = eigenvectors @ np.diag(eigenvalues_l) @ eigenvectors.T

        # Ensure symmetry
        L_approx = (L_approx + L_approx.T) / 2

        self._kernel_matrix = L_approx
        self._centered = False

        return L_approx, eigenvalues_l, eigenvectors


class GFTKernel(Kernel):
    """Graph Fourier Transform kernel with band-pass filtering.

    Provides frequency-selective spatial analysis by filtering the graph
    spectrum. Useful for isolating specific spatial scales:
    - Low frequencies: smooth, large-scale patterns
    - High frequencies: local fluctuations, boundaries

    Parameters
    ----------
    k_neighbors : int
        Number of nearest neighbors for graph construction
    rho : float
        ICAR autocorrelation strength
    freq_filter : str
        Frequency filter type: 'low', 'high', 'band', or 'all'
    low_cutoff : float
        Lower frequency cutoff (fraction of total, 0-1)
    high_cutoff : float
        Upper frequency cutoff (fraction of total, 0-1)
    """

    def __init__(
        self,
        k_neighbors: int = 20,
        rho: float = 0.99,
        freq_filter: str = 'all',
        low_cutoff: float = 0.0,
        high_cutoff: float = 1.0,
        regularization: float = 1e-6
    ):
        super().__init__()
        self.k_neighbors = k_neighbors
        self.rho = rho
        self.freq_filter = freq_filter
        self.low_cutoff = low_cutoff
        self.high_cutoff = high_cutoff
        self.regularization = regularization
        self._eigenvalues: Optional[NDArray] = None
        self._eigenvectors: Optional[NDArray] = None
        self._filter_mask: Optional[NDArray] = None

    def compute(self, coords: NDArray, use_gpu: bool = False, **kwargs) -> NDArray:
        """Compute GFT-filtered kernel matrix.

        Parameters
        ----------
        coords : ndarray
            Spatial coordinates of shape (n_spots, 2)
        use_gpu : bool
            If True, use GPU acceleration via CuPy

        Returns
        -------
        L_filtered : ndarray
            Frequency-filtered kernel matrix
        """
        # Build ICAR kernel
        icar = ICARKernel(
            k_neighbors=self.k_neighbors,
            rho=self.rho,
            regularization=self.regularization
        )
        L = icar.compute(coords)

        # Eigendecomposition (GPU accelerated if requested)
        if use_gpu and CUPY_AVAILABLE:
            eigenvalues, eigenvectors = gpu_eigh(L, use_gpu=True)
        else:
            eigenvalues, eigenvectors = np.linalg.eigh(L)

        # Sort by eigenvalue magnitude (ascending = low frequency first)
        idx = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        n = len(eigenvalues)

        # Create frequency filter mask
        freq_idx = np.arange(n) / n  # Normalized frequency index [0, 1]

        if self.freq_filter == 'low':
            mask = freq_idx <= self.high_cutoff
        elif self.freq_filter == 'high':
            mask = freq_idx >= self.low_cutoff
        elif self.freq_filter == 'band':
            mask = (freq_idx >= self.low_cutoff) & (freq_idx <= self.high_cutoff)
        else:  # 'all'
            mask = np.ones(n, dtype=bool)

        # Apply filter
        eigenvalues_filtered = eigenvalues.copy()
        eigenvalues_filtered[~mask] = 0

        # Reconstruct filtered kernel
        L_filtered = eigenvectors @ np.diag(eigenvalues_filtered) @ eigenvectors.T

        # Ensure symmetry
        L_filtered = (L_filtered + L_filtered.T) / 2

        self._kernel_matrix = L_filtered
        self._eigenvalues = eigenvalues
        self._eigenvectors = eigenvectors
        self._filter_mask = mask
        self._centered = False

        return L_filtered

    @property
    def eigenvalues(self) -> Optional[NDArray]:
        """Return eigenvalues of the full kernel (before filtering)."""
        return self._eigenvalues

    @property
    def eigenvectors(self) -> Optional[NDArray]:
        """Return eigenvectors of the full kernel."""
        return self._eigenvectors

    @property
    def filter_mask(self) -> Optional[NDArray]:
        """Return boolean mask indicating retained frequencies."""
        return self._filter_mask

    def get_frequency_spectrum(self) -> Optional[NDArray]:
        """Return the frequency spectrum (eigenvalue distribution)."""
        return self._eigenvalues


def compute_gft_features(
    coords: NDArray,
    data: NDArray,
    k_neighbors: int = 20,
    rho: float = 0.99,
    n_components: Optional[int] = None,
    use_gpu: bool = False
) -> Tuple[NDArray, NDArray]:
    """Compute Graph Fourier Transform features.

    Projects spatial data onto graph frequency basis for analysis
    of spatial patterns at different scales.

    Parameters
    ----------
    coords : ndarray
        Spatial coordinates of shape (n_spots, 2)
    data : ndarray
        Data matrix of shape (n_spots, n_features)
    k_neighbors : int
        Number of nearest neighbors for graph
    rho : float
        ICAR autocorrelation strength
    n_components : int, optional
        Number of frequency components to return (default: all)
    use_gpu : bool
        If True, use GPU acceleration via CuPy

    Returns
    -------
    gft_features : ndarray
        GFT coefficients of shape (n_components, n_features)
    frequencies : ndarray
        Corresponding graph frequencies (eigenvalues)
    """
    # Build ICAR kernel
    icar = ICARKernel(k_neighbors=k_neighbors, rho=rho)
    L = icar.compute(coords)

    # Eigendecomposition (GPU accelerated if requested)
    if use_gpu and CUPY_AVAILABLE:
        eigenvalues, eigenvectors = gpu_eigh(L, use_gpu=True)
    else:
        eigenvalues, eigenvectors = np.linalg.eigh(L)

    # Sort ascending (low frequency first)
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Ensure data is 2D
    if data.ndim == 1:
        data = data.reshape(-1, 1)

    # GFT: project data onto eigenvector basis
    gft_coeffs = eigenvectors.T @ data  # Shape: (n_spots, n_features)

    if n_components is not None:
        n_components = min(n_components, len(eigenvalues))
        gft_coeffs = gft_coeffs[:n_components]
        eigenvalues = eigenvalues[:n_components]

    return gft_coeffs, eigenvalues


def inverse_gft(
    gft_coeffs: NDArray,
    eigenvectors: NDArray,
    n_components: Optional[int] = None
) -> NDArray:
    """Compute inverse Graph Fourier Transform.

    Reconstructs spatial data from GFT coefficients.

    Parameters
    ----------
    gft_coeffs : ndarray
        GFT coefficients of shape (n_components, n_features)
    eigenvectors : ndarray
        Graph eigenvectors of shape (n_spots, n_spots)
    n_components : int, optional
        Number of components to use for reconstruction

    Returns
    -------
    data_reconstructed : ndarray
        Reconstructed spatial data
    """
    if n_components is not None:
        gft_coeffs = gft_coeffs[:n_components]
        eigenvectors = eigenvectors[:, :n_components]

    return eigenvectors @ gft_coeffs

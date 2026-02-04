"""Sparse matrix utilities for efficient large-scale computation."""

from typing import Tuple, Union

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import splu, spsolve

# Type alias for sparse matrices
SparseMatrix = Union[csr_matrix, csc_matrix]


def sparse_centering_trace(
    K: NDArray,
    L: SparseMatrix
) -> float:
    """Compute tr(H @ K @ H @ L) efficiently for sparse L.

    Uses the identity:
        tr(HKH @ L) = tr(K @ H @ L @ H)
                    = sum_ij K_ij * (HLH)_ij

    Parameters
    ----------
    K : ndarray
        Dense kernel matrix
    L : sparse matrix
        Sparse kernel matrix

    Returns
    -------
    trace : float
    """
    n = K.shape[0]

    # Center K: HKH
    row_means_K = K.mean(axis=1, keepdims=True)
    col_means_K = K.mean(axis=0, keepdims=True)
    grand_mean_K = K.mean()
    K_centered = K - row_means_K - col_means_K + grand_mean_K

    # Center L: HLH (need to convert to dense for centering)
    L_dense = L.toarray() if sparse.issparse(L) else L
    row_means_L = L_dense.mean(axis=1, keepdims=True)
    col_means_L = L_dense.mean(axis=0, keepdims=True)
    grand_mean_L = L_dense.mean()
    L_centered = L_dense - row_means_L - col_means_L + grand_mean_L

    # Trace = sum of element-wise product
    trace = np.sum(K_centered * L_centered)

    return trace


def sparse_inverse_column(
    Q: SparseMatrix,
    column_idx: int,
    lu_factorization=None
) -> NDArray:
    """Compute a single column of Q^{-1}.

    Parameters
    ----------
    Q : sparse matrix
        Matrix to invert
    column_idx : int
        Column index to compute
    lu_factorization : SuperLU, optional
        Pre-computed LU factorization

    Returns
    -------
    column : ndarray
    """
    n = Q.shape[0]
    e = np.zeros(n)
    e[column_idx] = 1.0

    if lu_factorization is not None:
        return lu_factorization.solve(e)
    else:
        return spsolve(Q.tocsc(), e)


def sparse_inverse_diagonal(
    Q: SparseMatrix,
    lu_factorization=None
) -> NDArray:
    """Compute diagonal of Q^{-1}.

    More efficient than computing full inverse when only
    diagonal is needed (e.g., for variance estimates).

    Parameters
    ----------
    Q : sparse matrix
        Matrix to invert
    lu_factorization : SuperLU, optional
        Pre-computed LU factorization

    Returns
    -------
    diagonal : ndarray
    """
    n = Q.shape[0]

    if lu_factorization is None:
        lu_factorization = splu(Q.tocsc())

    diagonal = np.zeros(n)
    for i in range(n):
        e = np.zeros(n)
        e[i] = 1.0
        diagonal[i] = lu_factorization.solve(e)[i]

    return diagonal


def efficient_knn_adjacency(
    coords: NDArray,
    k: int,
    max_chunk_size: int = 10000
) -> csr_matrix:
    """Build k-NN adjacency matrix efficiently for large datasets.

    Uses chunked processing to manage memory.

    Parameters
    ----------
    coords : ndarray
        Coordinates (n, 2)
    k : int
        Number of neighbors
    max_chunk_size : int
        Maximum number of points per chunk

    Returns
    -------
    W : csr_matrix
        Symmetric adjacency matrix
    """
    from scipy.spatial import KDTree

    n = coords.shape[0]
    k = min(k, n - 1)

    tree = KDTree(coords)

    row_indices = []
    col_indices = []

    # Process in chunks
    n_chunks = (n + max_chunk_size - 1) // max_chunk_size

    for chunk_idx in range(n_chunks):
        start = chunk_idx * max_chunk_size
        end = min(start + max_chunk_size, n)
        chunk_coords = coords[start:end]

        # Query neighbors
        _, indices = tree.query(chunk_coords, k=k + 1)

        for i, local_idx in enumerate(range(start, end)):
            for j_idx in range(1, k + 1):  # Skip self
                j = indices[i, j_idx]
                row_indices.append(local_idx)
                col_indices.append(j)

    # Create sparse matrix
    data = np.ones(len(row_indices))
    W = csr_matrix((data, (row_indices, col_indices)), shape=(n, n))

    # Symmetrize
    W = W + W.T
    W.data = np.ones_like(W.data)

    return W


def compute_eigenvalues_sparse(
    Q: SparseMatrix,
    k: int = 100,
    which: str = 'SM'
) -> NDArray:
    """Compute eigenvalues of sparse matrix.

    Parameters
    ----------
    Q : sparse matrix
        Matrix for eigendecomposition
    k : int
        Number of eigenvalues to compute
    which : str
        Which eigenvalues: 'SM' (smallest magnitude),
        'LM' (largest magnitude), etc.

    Returns
    -------
    eigenvalues : ndarray
    """
    from scipy.sparse.linalg import eigsh

    # Ensure symmetric
    Q_sym = (Q + Q.T) / 2

    eigenvalues, _ = eigsh(Q_sym, k=k, which=which)

    return eigenvalues

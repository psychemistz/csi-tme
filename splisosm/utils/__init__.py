"""Utility functions for SPLISOSM."""

from .sparse import (
    sparse_centering_trace,
    sparse_inverse_column,
    sparse_inverse_diagonal,
    efficient_knn_adjacency,
    compute_eigenvalues_sparse,
)
from .backend import (
    is_gpu_available,
    get_backend,
    get_array_module,
    to_numpy,
    to_backend,
    ensure_same_backend,
    BackendContext,
    gpu_eigvalsh,
    gpu_eigh,
    gpu_matmul,
    gpu_kernel_compute,
    gpu_trace_product,
    gpu_center_kernel,
    get_gpu_memory_info,
    clear_gpu_memory,
    CUPY_AVAILABLE,
)

__all__ = [
    # Sparse utilities
    "sparse_centering_trace",
    "sparse_inverse_column",
    "sparse_inverse_diagonal",
    "efficient_knn_adjacency",
    "compute_eigenvalues_sparse",
    # GPU backend
    "is_gpu_available",
    "get_backend",
    "get_array_module",
    "to_numpy",
    "to_backend",
    "ensure_same_backend",
    "BackendContext",
    "gpu_eigvalsh",
    "gpu_eigh",
    "gpu_matmul",
    "gpu_kernel_compute",
    "gpu_trace_product",
    "gpu_center_kernel",
    "get_gpu_memory_info",
    "clear_gpu_memory",
    "CUPY_AVAILABLE",
]

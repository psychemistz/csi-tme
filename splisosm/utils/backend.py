"""Backend abstraction for NumPy/CuPy GPU acceleration.

Provides seamless switching between CPU (NumPy) and GPU (CuPy) backends.
CuPy is optional - falls back to NumPy if not available.

Usage:
    from splisosm.utils.backend import get_backend, to_numpy, to_backend

    xp = get_backend(use_gpu=True)  # Returns cupy if available, else numpy
    arr = xp.array([1, 2, 3])
    result = to_numpy(arr)  # Always returns numpy array
"""

from typing import Any, Optional, Union
import numpy as np
from numpy.typing import NDArray

# Try to import CuPy
try:
    import cupy as cp
    import cupyx.scipy.sparse as cp_sparse
    import cupyx.scipy.linalg as cp_linalg
    CUPY_AVAILABLE = True
except ImportError:
    cp = None
    cp_sparse = None
    cp_linalg = None
    CUPY_AVAILABLE = False

# Type alias for array module
ArrayModule = Any  # numpy or cupy module


def is_gpu_available() -> bool:
    """Check if GPU acceleration is available via CuPy."""
    return CUPY_AVAILABLE


def get_backend(use_gpu: bool = False) -> ArrayModule:
    """Get the appropriate array backend.

    Parameters
    ----------
    use_gpu : bool
        If True, return CuPy if available, else NumPy.
        If False, always return NumPy.

    Returns
    -------
    xp : module
        Either numpy or cupy module
    """
    if use_gpu and CUPY_AVAILABLE:
        return cp
    return np


def get_array_module(arr: Any) -> ArrayModule:
    """Get the array module for an array (numpy or cupy).

    Parameters
    ----------
    arr : array-like
        Input array

    Returns
    -------
    xp : module
        The module that owns the array
    """
    if CUPY_AVAILABLE:
        return cp.get_array_module(arr)
    return np


def to_numpy(arr: Any) -> NDArray:
    """Convert array to numpy (CPU) array.

    Parameters
    ----------
    arr : array-like
        Input array (numpy or cupy)

    Returns
    -------
    result : ndarray
        NumPy array on CPU
    """
    if CUPY_AVAILABLE and isinstance(arr, cp.ndarray):
        return cp.asnumpy(arr)
    return np.asarray(arr)


def to_backend(arr: NDArray, use_gpu: bool = False) -> Any:
    """Convert numpy array to specified backend.

    Parameters
    ----------
    arr : ndarray
        NumPy array
    use_gpu : bool
        If True, convert to CuPy array

    Returns
    -------
    result : array
        Array on specified backend
    """
    if use_gpu and CUPY_AVAILABLE:
        return cp.asarray(arr)
    return np.asarray(arr)


def ensure_same_backend(*arrays) -> tuple:
    """Ensure all arrays are on the same backend.

    If any array is on GPU, move all to GPU.
    Otherwise, keep all on CPU.

    Parameters
    ----------
    *arrays : array-like
        Input arrays

    Returns
    -------
    tuple
        Arrays all on the same backend
    """
    if not CUPY_AVAILABLE:
        return tuple(np.asarray(a) for a in arrays)

    # Check if any array is on GPU
    has_gpu = any(isinstance(a, cp.ndarray) for a in arrays)

    if has_gpu:
        return tuple(cp.asarray(a) for a in arrays)
    return tuple(np.asarray(a) for a in arrays)


class BackendContext:
    """Context manager for backend operations.

    Provides convenient access to backend-specific functions.

    Parameters
    ----------
    use_gpu : bool
        Whether to use GPU acceleration

    Examples
    --------
    >>> with BackendContext(use_gpu=True) as ctx:
    ...     arr = ctx.xp.array([1, 2, 3])
    ...     result = ctx.xp.linalg.norm(arr)
    """

    def __init__(self, use_gpu: bool = False):
        self.use_gpu = use_gpu and CUPY_AVAILABLE
        self.xp = get_backend(use_gpu)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Synchronize GPU if needed
        if self.use_gpu and CUPY_AVAILABLE:
            cp.cuda.Stream.null.synchronize()
        return False

    @property
    def linalg(self):
        """Get linalg module for current backend."""
        if self.use_gpu:
            return cp.linalg
        return np.linalg

    @property
    def scipy_linalg(self):
        """Get scipy.linalg equivalent for current backend."""
        if self.use_gpu:
            return cp_linalg
        from scipy import linalg
        return linalg

    def to_numpy(self, arr: Any) -> NDArray:
        """Convert array to numpy."""
        return to_numpy(arr)

    def asarray(self, arr: Any) -> Any:
        """Convert to current backend's array type."""
        return self.xp.asarray(arr)


def gpu_eigvalsh(K: Any, use_gpu: bool = False) -> NDArray:
    """Compute eigenvalues of symmetric matrix with optional GPU.

    Parameters
    ----------
    K : array-like
        Symmetric matrix
    use_gpu : bool
        Whether to use GPU

    Returns
    -------
    eigenvalues : ndarray
        Eigenvalues (always returned as numpy array)
    """
    if use_gpu and CUPY_AVAILABLE:
        K_gpu = cp.asarray(K)
        eigenvalues = cp.linalg.eigvalsh(K_gpu)
        return cp.asnumpy(eigenvalues)
    return np.linalg.eigvalsh(K)


def gpu_eigh(K: Any, use_gpu: bool = False) -> tuple:
    """Compute eigendecomposition of symmetric matrix with optional GPU.

    Parameters
    ----------
    K : array-like
        Symmetric matrix
    use_gpu : bool
        Whether to use GPU

    Returns
    -------
    eigenvalues : ndarray
        Eigenvalues (numpy array)
    eigenvectors : ndarray
        Eigenvectors (numpy array)
    """
    if use_gpu and CUPY_AVAILABLE:
        K_gpu = cp.asarray(K)
        eigenvalues, eigenvectors = cp.linalg.eigh(K_gpu)
        return cp.asnumpy(eigenvalues), cp.asnumpy(eigenvectors)
    return np.linalg.eigh(K)


def gpu_matmul(A: Any, B: Any, use_gpu: bool = False) -> NDArray:
    """Matrix multiplication with optional GPU.

    Parameters
    ----------
    A, B : array-like
        Input matrices
    use_gpu : bool
        Whether to use GPU

    Returns
    -------
    result : ndarray
        Matrix product (numpy array)
    """
    if use_gpu and CUPY_AVAILABLE:
        A_gpu = cp.asarray(A)
        B_gpu = cp.asarray(B)
        result = A_gpu @ B_gpu
        return cp.asnumpy(result)
    return np.asarray(A) @ np.asarray(B)


def gpu_kernel_compute(X: Any, use_gpu: bool = False) -> NDArray:
    """Compute linear kernel X @ X.T with optional GPU.

    Parameters
    ----------
    X : array-like
        Feature matrix of shape (n, p)
    use_gpu : bool
        Whether to use GPU

    Returns
    -------
    K : ndarray
        Kernel matrix of shape (n, n)
    """
    if use_gpu and CUPY_AVAILABLE:
        X_gpu = cp.asarray(X)
        K = X_gpu @ X_gpu.T
        return cp.asnumpy(K)
    X = np.asarray(X)
    return X @ X.T


def gpu_trace_product(A: Any, B: Any, use_gpu: bool = False) -> float:
    """Compute tr(A @ B) efficiently as sum(A * B.T).

    Parameters
    ----------
    A, B : array-like
        Input matrices
    use_gpu : bool
        Whether to use GPU

    Returns
    -------
    trace : float
        Trace of the product
    """
    if use_gpu and CUPY_AVAILABLE:
        A_gpu = cp.asarray(A)
        B_gpu = cp.asarray(B)
        # tr(A @ B) = sum(A * B.T) for symmetric B
        # For general case: sum(A * B.T)
        result = cp.sum(A_gpu * B_gpu)
        return float(cp.asnumpy(result))
    return float(np.sum(np.asarray(A) * np.asarray(B)))


def gpu_center_kernel(K: Any, use_gpu: bool = False) -> NDArray:
    """Center kernel matrix: HKH with optional GPU.

    Parameters
    ----------
    K : array-like
        Kernel matrix of shape (n, n)
    use_gpu : bool
        Whether to use GPU

    Returns
    -------
    K_centered : ndarray
        Centered kernel matrix
    """
    if use_gpu and CUPY_AVAILABLE:
        K_gpu = cp.asarray(K)
        n = K_gpu.shape[0]
        row_means = K_gpu.mean(axis=1, keepdims=True)
        col_means = K_gpu.mean(axis=0, keepdims=True)
        grand_mean = K_gpu.mean()
        K_centered = K_gpu - row_means - col_means + grand_mean
        return cp.asnumpy(K_centered)

    K = np.asarray(K)
    row_means = K.mean(axis=1, keepdims=True)
    col_means = K.mean(axis=0, keepdims=True)
    grand_mean = K.mean()
    return K - row_means - col_means + grand_mean


def get_gpu_memory_info() -> Optional[dict]:
    """Get GPU memory information if available.

    Returns
    -------
    info : dict or None
        Dictionary with 'free', 'total', 'used' in bytes, or None if no GPU
    """
    if not CUPY_AVAILABLE:
        return None

    try:
        mempool = cp.get_default_memory_pool()
        free, total = cp.cuda.runtime.memGetInfo()
        return {
            'free': free,
            'total': total,
            'used': total - free,
            'pool_used': mempool.used_bytes(),
            'pool_total': mempool.total_bytes(),
        }
    except Exception:
        return None


def clear_gpu_memory():
    """Clear GPU memory pool if CuPy is available."""
    if CUPY_AVAILABLE:
        mempool = cp.get_default_memory_pool()
        mempool.free_all_blocks()

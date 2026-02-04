"""Boundary-weighted kernel modifications for tumor boundary analysis."""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import KDTree

from .base import Kernel


def compute_boundary_scores(
    coords: NDArray,
    region_labels: NDArray,
    tumor_label: int = 1,
    stroma_label: int = 0,
    smoothing_k: int = 10,
    interior_weight: float = 0.1
) -> NDArray:
    """Compute boundary proximity scores for each spot.

    Spots near the tumor-stroma interface get high scores (close to 1),
    spots in the interior get low scores (interior_weight).

    Parameters
    ----------
    coords : ndarray
        Spatial coordinates of shape (n_spots, 2)
    region_labels : ndarray
        Region label for each spot (e.g., 0=stroma, 1=tumor)
    tumor_label : int
        Label value for tumor regions
    stroma_label : int
        Label value for stroma regions
    smoothing_k : int
        Number of neighbors for smoothing boundary detection
    interior_weight : float
        Minimum weight for interior spots

    Returns
    -------
    boundary_scores : ndarray
        Boundary proximity scores in [interior_weight, 1]
    """
    n = len(coords)
    tree = KDTree(coords)

    # Find neighbors for each spot
    _, indices = tree.query(coords, k=smoothing_k + 1)

    # Compute fraction of neighbors with different label
    boundary_scores = np.zeros(n)

    for i in range(n):
        neighbors = indices[i, 1:]  # Exclude self
        neighbor_labels = region_labels[neighbors]
        my_label = region_labels[i]

        # Fraction of neighbors with different label
        diff_fraction = np.mean(neighbor_labels != my_label)
        boundary_scores[i] = diff_fraction

    # Normalize to [interior_weight, 1]
    if boundary_scores.max() > 0:
        boundary_scores = (boundary_scores / boundary_scores.max() *
                          (1 - interior_weight) + interior_weight)
    else:
        boundary_scores = np.ones(n) * interior_weight

    return boundary_scores


def compute_distance_to_interface(
    coords: NDArray,
    region_labels: NDArray,
    target_label: int = 1
) -> NDArray:
    """Compute distance from each spot to nearest interface.

    Parameters
    ----------
    coords : ndarray
        Spatial coordinates of shape (n_spots, 2)
    region_labels : ndarray
        Region label for each spot
    target_label : int
        Label to measure distance from

    Returns
    -------
    distances : ndarray
        Distance to nearest spot with different label
    """
    target_mask = region_labels == target_label
    other_mask = ~target_mask

    if not np.any(other_mask):
        return np.zeros(len(coords))

    # Build tree of other-region spots
    other_coords = coords[other_mask]
    tree = KDTree(other_coords)

    # Query distance to nearest other-region spot
    distances, _ = tree.query(coords, k=1)

    return distances


class BoundaryWeightedKernel(Kernel):
    """Kernel weighted by boundary proximity.

    Modifies a base kernel by element-wise multiplication with
    boundary weight matrix: L_boundary = L ⊙ B

    Parameters
    ----------
    base_kernel : Kernel
        Base kernel to modify
    interior_weight : float
        Minimum weight for interior spots
    """

    def __init__(self, base_kernel: Kernel, interior_weight: float = 0.1):
        super().__init__()
        self.base_kernel = base_kernel
        self.interior_weight = interior_weight
        self._boundary_scores: Optional[NDArray] = None

    def compute(self, coords: NDArray, region_labels: NDArray,
                **kwargs) -> NDArray:
        """Compute boundary-weighted kernel.

        Parameters
        ----------
        coords : ndarray
            Spatial coordinates
        region_labels : ndarray
            Region labels for boundary detection

        Returns
        -------
        K : ndarray
            Boundary-weighted kernel matrix
        """
        # Compute base kernel
        L = self.base_kernel.compute(coords, **kwargs)

        # Compute boundary scores
        boundary_scores = compute_boundary_scores(
            coords, region_labels,
            interior_weight=self.interior_weight
        )
        self._boundary_scores = boundary_scores

        # Build boundary weight matrix: B_ij = sqrt(b_i * b_j)
        B = np.sqrt(np.outer(boundary_scores, boundary_scores))

        # Apply weighting
        K = L * B

        self._kernel_matrix = K
        self._centered = False
        return K

    @property
    def boundary_scores(self) -> Optional[NDArray]:
        """Return computed boundary scores."""
        return self._boundary_scores

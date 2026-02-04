"""Kernel implementations for SPLISOSM spatial analysis."""

from .base import Kernel, LinearKernel
from .icar import (
    ICARKernel,
    GFTKernel,
    compute_gft_features,
    inverse_gft,
)
from .gaussian import GaussianKernel
from .compositional import (
    CompositionalKernel,
    ExpressionCompositionalKernel,
    clr_transform,
    ilr_transform,
    alr_transform,
    radial_transform,
)
from .boundary import (
    BoundaryWeightedKernel,
    compute_boundary_scores,
    compute_distance_to_interface,
)

__all__ = [
    "Kernel",
    "LinearKernel",
    "ICARKernel",
    "GFTKernel",
    "compute_gft_features",
    "inverse_gft",
    "GaussianKernel",
    "CompositionalKernel",
    "ExpressionCompositionalKernel",
    "BoundaryWeightedKernel",
    "clr_transform",
    "ilr_transform",
    "alr_transform",
    "radial_transform",
    "compute_boundary_scores",
    "compute_distance_to_interface",
]

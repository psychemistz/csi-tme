"""Statistical tests for SPLISOSM spatial isoform analysis."""

from .hsic import (
    HSICResult,
    center_kernel,
    compute_hsic_statistic,
    compute_hsic_gc,
    compute_hsic_ir,
    compute_hsic_ic,
    compute_hsic_combined,
    get_kernel_eigenvalues,
)
from .pvalue import (
    PValueResult,
    liu_sf,
    liu_pvalue,
    hsic_pvalue,
    compute_composite_eigenvalues,
    gamma_approximation_pvalue,
    compute_pvalue_from_kernels,
    permutation_pvalue,
)
from .conditional import (
    ConditionalHSICResult,
    spatial_residualize,
    conditional_hsic,
    partial_correlation_spatial,
)

__all__ = [
    # HSIC
    "HSICResult",
    "center_kernel",
    "compute_hsic_statistic",
    "compute_hsic_gc",
    "compute_hsic_ir",
    "compute_hsic_ic",
    "compute_hsic_combined",
    "get_kernel_eigenvalues",
    # P-values
    "PValueResult",
    "liu_sf",
    "liu_pvalue",
    "hsic_pvalue",
    "compute_composite_eigenvalues",
    "gamma_approximation_pvalue",
    "compute_pvalue_from_kernels",
    "permutation_pvalue",
    # Conditional
    "ConditionalHSICResult",
    "spatial_residualize",
    "conditional_hsic",
    "partial_correlation_spatial",
]

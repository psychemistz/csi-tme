"""Validation module for SPLISOSM results.

Provides integration with single-cell validation methods like SCALPEL
for validating spatially variable isoform patterns.
"""

from .scalpel import (
    load_scalpel_results,
    validate_svp_genes,
    compute_isoform_concordance,
)

__all__ = [
    "load_scalpel_results",
    "validate_svp_genes",
    "compute_isoform_concordance",
]

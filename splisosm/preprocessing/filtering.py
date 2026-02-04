"""Quality filtering for spots and genes."""

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import NDArray


@dataclass
class FilterResult:
    """Result of filtering operation.

    Attributes
    ----------
    mask : ndarray
        Boolean mask of items passing filter
    n_passed : int
        Number of items passing
    n_failed : int
        Number of items failing
    """
    mask: NDArray
    n_passed: int
    n_failed: int

    @property
    def pass_rate(self) -> float:
        """Fraction of items passing filter."""
        return self.n_passed / (self.n_passed + self.n_failed)


def filter_spots_by_counts(
    counts: NDArray,
    min_counts: int = 0,
    min_genes: int = 0
) -> FilterResult:
    """Filter spots by total counts and gene detection.

    Parameters
    ----------
    counts : ndarray
        Count matrix of shape (n_spots, n_genes)
    min_counts : int
        Minimum total UMI counts per spot
    min_genes : int
        Minimum number of detected genes per spot

    Returns
    -------
    result : FilterResult
    """
    total_counts = counts.sum(axis=1)
    n_genes = (counts > 0).sum(axis=1)

    mask = (total_counts >= min_counts) & (n_genes >= min_genes)

    return FilterResult(
        mask=mask,
        n_passed=mask.sum(),
        n_failed=(~mask).sum()
    )


def filter_genes_by_expression(
    counts: NDArray,
    min_spots: int = 0,
    min_counts: int = 0
) -> FilterResult:
    """Filter genes by expression across spots.

    Parameters
    ----------
    counts : ndarray
        Count matrix of shape (n_spots, n_genes)
    min_spots : int
        Minimum number of spots with non-zero expression
    min_counts : int
        Minimum total counts across all spots

    Returns
    -------
    result : FilterResult
    """
    n_spots_expressed = (counts > 0).sum(axis=0)
    total_counts = counts.sum(axis=0)

    mask = (n_spots_expressed >= min_spots) & (total_counts >= min_counts)

    return FilterResult(
        mask=mask,
        n_passed=mask.sum(),
        n_failed=(~mask).sum()
    )


def filter_by_isoform_coverage(
    isoform_counts: NDArray,
    gene_totals: NDArray,
    min_spots_with_expression: int = 50,
    min_fraction_detected: float = 0.1
) -> Tuple[bool, dict]:
    """Check if a gene has sufficient isoform coverage for analysis.

    Parameters
    ----------
    isoform_counts : ndarray
        Isoform counts of shape (n_spots, n_isoforms)
    gene_totals : ndarray
        Total gene counts per spot
    min_spots_with_expression : int
        Minimum number of spots with gene expression
    min_fraction_detected : float
        Minimum fraction of isoforms that must be detected

    Returns
    -------
    passes : bool
        Whether gene passes filters
    diagnostics : dict
        Filtering diagnostics
    """
    n_spots, n_isoforms = isoform_counts.shape

    # Spots with any expression
    spots_expressed = gene_totals > 0
    n_spots_expressed = spots_expressed.sum()

    # Isoforms with any counts
    isoform_totals = isoform_counts.sum(axis=0)
    isoforms_detected = isoform_totals > 0
    n_isoforms_detected = isoforms_detected.sum()
    fraction_detected = n_isoforms_detected / n_isoforms

    diagnostics = {
        'n_spots_expressed': int(n_spots_expressed),
        'n_isoforms_detected': int(n_isoforms_detected),
        'fraction_detected': fraction_detected,
    }

    passes = (n_spots_expressed >= min_spots_with_expression and
              fraction_detected >= min_fraction_detected)

    return passes, diagnostics


def identify_variable_isoforms(
    ratios: NDArray,
    min_variance: float = 0.01,
    min_mean: float = 0.01
) -> NDArray:
    """Identify isoforms with sufficient variability for analysis.

    Filters out isoforms that are nearly constant or rarely expressed.

    Parameters
    ----------
    ratios : ndarray
        Isoform ratios of shape (n_spots, n_isoforms)
    min_variance : float
        Minimum variance across spots
    min_mean : float
        Minimum mean ratio

    Returns
    -------
    mask : ndarray
        Boolean mask of variable isoforms
    """
    # Compute statistics ignoring NaN
    means = np.nanmean(ratios, axis=0)
    variances = np.nanvar(ratios, axis=0)

    mask = (variances >= min_variance) & (means >= min_mean)

    return mask

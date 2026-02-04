"""Missing ratio imputation for compositional data.

Implements Theorem 2 from MATH.md: Mean imputation preserves test validity
because the centering matrix H projects out the contribution of imputed values.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray


def compute_isoform_ratios(
    isoform_counts: NDArray,
    gene_totals: Optional[NDArray] = None,
    min_total: float = 0.0
) -> Tuple[NDArray, NDArray]:
    """Compute isoform ratios from counts.

    Parameters
    ----------
    isoform_counts : ndarray
        Isoform counts of shape (n_spots, n_isoforms)
    gene_totals : ndarray, optional
        Total gene counts per spot. If None, computed from isoform_counts.
    min_total : float
        Minimum total count to consider valid (spots below get NaN)

    Returns
    -------
    ratios : ndarray
        Isoform ratios of shape (n_spots, n_isoforms).
        Rows with zero totals have NaN values.
    valid_mask : ndarray
        Boolean mask of valid spots (total > min_total)
    """
    if gene_totals is None:
        gene_totals = isoform_counts.sum(axis=1)

    valid_mask = gene_totals > min_total

    # Initialize with NaN
    ratios = np.full_like(isoform_counts, np.nan, dtype=float)

    # Compute ratios for valid spots
    ratios[valid_mask] = (isoform_counts[valid_mask].T /
                          gene_totals[valid_mask]).T

    return ratios, valid_mask


def mean_impute_ratios(
    ratios: NDArray,
    valid_mask: Optional[NDArray] = None
) -> NDArray:
    """Impute missing ratios with global mean.

    Per Theorem 2: Mean imputation preserves test validity because
    the centering matrix H projects out the contribution of
    constant (imputed) values to the test statistic.

    Parameters
    ----------
    ratios : ndarray
        Isoform ratios of shape (n_spots, n_isoforms).
        Missing values indicated by NaN.
    valid_mask : ndarray, optional
        Boolean mask of valid spots. If None, inferred from NaN.

    Returns
    -------
    imputed_ratios : ndarray
        Ratios with NaN replaced by column means
    """
    if valid_mask is None:
        valid_mask = ~np.isnan(ratios[:, 0])

    imputed = ratios.copy()

    # Compute mean from valid spots
    global_mean = np.nanmean(ratios[valid_mask], axis=0)

    # Impute missing spots
    imputed[~valid_mask] = global_mean

    return imputed


def impute_with_counts(
    isoform_counts: NDArray,
    gene_totals: Optional[NDArray] = None
) -> NDArray:
    """Compute imputed ratios directly from counts.

    For spots with zero counts, uses global ratio computed from
    total isoform counts across all valid spots.

    Parameters
    ----------
    isoform_counts : ndarray
        Isoform counts of shape (n_spots, n_isoforms)
    gene_totals : ndarray, optional
        Total gene counts per spot

    Returns
    -------
    imputed_ratios : ndarray
        Isoform ratios with valid imputation
    """
    if gene_totals is None:
        gene_totals = isoform_counts.sum(axis=1)

    valid_mask = gene_totals > 0

    # Global mean ratio from valid spots
    total_isoform_counts = isoform_counts[valid_mask].sum(axis=0)
    global_mean = total_isoform_counts / total_isoform_counts.sum()

    # Compute ratios
    ratios = np.zeros_like(isoform_counts, dtype=float)

    # Valid spots: actual ratios
    ratios[valid_mask] = (isoform_counts[valid_mask].T /
                          gene_totals[valid_mask]).T

    # Invalid spots: global mean
    ratios[~valid_mask] = global_mean

    return ratios


def weighted_mean_impute(
    ratios: NDArray,
    weights: NDArray,
    valid_mask: Optional[NDArray] = None
) -> NDArray:
    """Impute missing ratios with weighted mean.

    Uses weights (e.g., total counts) to compute more representative
    global mean from valid spots.

    Parameters
    ----------
    ratios : ndarray
        Isoform ratios of shape (n_spots, n_isoforms)
    weights : ndarray
        Weights for each spot (e.g., total gene counts)
    valid_mask : ndarray, optional
        Boolean mask of valid spots

    Returns
    -------
    imputed_ratios : ndarray
    """
    if valid_mask is None:
        valid_mask = ~np.isnan(ratios[:, 0])

    imputed = ratios.copy()

    # Weighted mean from valid spots
    valid_weights = weights[valid_mask]
    valid_ratios = ratios[valid_mask]
    weighted_mean = np.average(valid_ratios, weights=valid_weights, axis=0)

    # Impute
    imputed[~valid_mask] = weighted_mean

    return imputed

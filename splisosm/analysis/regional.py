"""Region-stratified analysis for tumor microenvironment zones."""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from .single_sample import SplisomAnalyzer, GeneResult


@dataclass
class RegionalResult:
    """Result comparing SPLISOSM statistics across regions.

    Attributes
    ----------
    gene : str
        Gene name
    global_result : GeneResult
        Result from full tissue analysis
    tumor_result : GeneResult
        Result from tumor region only
    stroma_result : GeneResult
        Result from stroma region only
    boundary_result : GeneResult
        Result from boundary zone
    regional_enrichment : dict
        Enrichment scores for each region
    """
    gene: str
    global_result: Optional[GeneResult]
    tumor_result: Optional[GeneResult]
    stroma_result: Optional[GeneResult]
    boundary_result: Optional[GeneResult]
    regional_enrichment: Optional[dict] = None


def identify_boundary_spots(
    region_labels: NDArray,
    coords: NDArray,
    boundary_distance: float = 2.0
) -> NDArray:
    """Identify spots near region boundaries.

    Parameters
    ----------
    region_labels : ndarray
        Region label for each spot
    coords : ndarray
        Spatial coordinates
    boundary_distance : float
        Maximum distance from boundary to be considered boundary zone

    Returns
    -------
    boundary_mask : ndarray
        Boolean mask of boundary spots
    """
    from scipy.spatial import KDTree

    unique_labels = np.unique(region_labels)
    if len(unique_labels) < 2:
        return np.zeros(len(coords), dtype=bool)

    boundary_mask = np.zeros(len(coords), dtype=bool)

    # For each label, find spots near the interface with other labels
    for label in unique_labels:
        label_mask = region_labels == label
        other_mask = ~label_mask

        if not np.any(other_mask):
            continue

        # Build tree of other-region spots
        tree = KDTree(coords[other_mask])

        # Query distances from this region to other regions
        distances, _ = tree.query(coords[label_mask], k=1)

        # Mark spots within boundary distance
        near_boundary = distances <= boundary_distance
        boundary_mask[label_mask] = near_boundary

    return boundary_mask


class RegionalAnalyzer:
    """Analyzer for region-stratified SPLISOSM analysis.

    Parameters
    ----------
    base_analyzer : SplisomAnalyzer
        Base analyzer with kernel parameters
    min_spots_per_region : int
        Minimum spots required for region analysis
    """

    def __init__(
        self,
        base_analyzer: Optional[SplisomAnalyzer] = None,
        min_spots_per_region: int = 30
    ):
        if base_analyzer is None:
            base_analyzer = SplisomAnalyzer()
        self.base_analyzer = base_analyzer
        self.min_spots_per_region = min_spots_per_region

    def analyze_by_region(
        self,
        isoform_counts: NDArray,
        coords: NDArray,
        region_labels: NDArray,
        gene_name: str,
        tumor_label: int = 1,
        stroma_label: int = 0
    ) -> RegionalResult:
        """Analyze a gene in different tissue regions.

        Parameters
        ----------
        isoform_counts : ndarray
            Isoform counts (n_spots, n_isoforms)
        coords : ndarray
            Spatial coordinates
        region_labels : ndarray
            Region label per spot
        gene_name : str
            Gene name
        tumor_label : int
            Label value for tumor region
        stroma_label : int
            Label value for stroma region

        Returns
        -------
        result : RegionalResult
        """
        # Define region masks
        tumor_mask = region_labels == tumor_label
        stroma_mask = region_labels == stroma_label
        boundary_mask = identify_boundary_spots(region_labels, coords)

        results = {'global': None, 'tumor': None, 'stroma': None, 'boundary': None}

        # Global analysis
        self.base_analyzer.build_spatial_kernel(coords)
        results['global'] = self.base_analyzer.analyze_gene(
            isoform_counts, gene_name
        )

        # Tumor region
        if tumor_mask.sum() >= self.min_spots_per_region:
            self.base_analyzer.build_spatial_kernel(coords[tumor_mask])
            results['tumor'] = self.base_analyzer.analyze_gene(
                isoform_counts[tumor_mask], gene_name
            )

        # Stroma region
        if stroma_mask.sum() >= self.min_spots_per_region:
            self.base_analyzer.build_spatial_kernel(coords[stroma_mask])
            results['stroma'] = self.base_analyzer.analyze_gene(
                isoform_counts[stroma_mask], gene_name
            )

        # Boundary region
        if boundary_mask.sum() >= self.min_spots_per_region:
            self.base_analyzer.build_spatial_kernel(coords[boundary_mask])
            results['boundary'] = self.base_analyzer.analyze_gene(
                isoform_counts[boundary_mask], gene_name
            )

        # Compute enrichment scores
        enrichment = self._compute_enrichment(results)

        return RegionalResult(
            gene=gene_name,
            global_result=results['global'],
            tumor_result=results['tumor'],
            stroma_result=results['stroma'],
            boundary_result=results['boundary'],
            regional_enrichment=enrichment
        )

    def _compute_enrichment(self, results: dict) -> dict:
        """Compute regional enrichment scores.

        Higher enrichment = stronger spatial isoform signal in that region
        compared to global.
        """
        enrichment = {}

        global_stat = (results['global'].statistic_ir
                       if results['global'] else 0)

        for region in ['tumor', 'stroma', 'boundary']:
            if results[region] and global_stat > 0:
                region_stat = results[region].statistic_ir
                enrichment[region] = region_stat / global_stat
            else:
                enrichment[region] = None

        return enrichment


def compare_tumor_stroma_ratios(
    isoform_counts: NDArray,
    region_labels: NDArray,
    tumor_label: int = 1,
    stroma_label: int = 0
) -> dict:
    """Compare isoform ratios between tumor and stroma.

    Parameters
    ----------
    isoform_counts : ndarray
        Isoform counts (n_spots, n_isoforms)
    region_labels : ndarray
        Region labels
    tumor_label : int
        Tumor region label
    stroma_label : int
        Stroma region label

    Returns
    -------
    comparison : dict
        Contains mean ratios and statistical comparison
    """
    from scipy import stats as scipy_stats

    from ..preprocessing.imputation import impute_with_counts

    tumor_mask = region_labels == tumor_label
    stroma_mask = region_labels == stroma_label

    # Compute ratios
    ratios = impute_with_counts(isoform_counts)

    tumor_ratios = ratios[tumor_mask]
    stroma_ratios = ratios[stroma_mask]

    comparison = {
        'tumor_mean': tumor_ratios.mean(axis=0),
        'stroma_mean': stroma_ratios.mean(axis=0),
        'tumor_n': tumor_mask.sum(),
        'stroma_n': stroma_mask.sum(),
    }

    # Statistical test for each isoform
    n_isoforms = ratios.shape[1]
    pvalues = []

    for i in range(n_isoforms):
        stat, pval = scipy_stats.mannwhitneyu(
            tumor_ratios[:, i], stroma_ratios[:, i],
            alternative='two-sided'
        )
        pvalues.append(pval)

    comparison['pvalues'] = np.array(pvalues)
    comparison['effect_size'] = comparison['tumor_mean'] - comparison['stroma_mean']

    return comparison

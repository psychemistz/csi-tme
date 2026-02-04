"""Analysis pipelines for SPLISOSM."""

from .single_sample import (
    GeneResult,
    SampleResult,
    SplisomAnalyzer,
    build_isoform_mapping,
)
from .meta import (
    MetaAnalysisResult,
    fisher_combine,
    stouffer_combine,
    hierarchical_fdr,
    meta_analyze_results,
    results_to_dataframe,
)
from .regional import (
    RegionalResult,
    RegionalAnalyzer,
    identify_boundary_spots,
    compare_tumor_stroma_ratios,
)

__all__ = [
    # Single sample
    "GeneResult",
    "SampleResult",
    "SplisomAnalyzer",
    "build_isoform_mapping",
    # Meta-analysis
    "MetaAnalysisResult",
    "fisher_combine",
    "stouffer_combine",
    "hierarchical_fdr",
    "meta_analyze_results",
    "results_to_dataframe",
    # Regional
    "RegionalResult",
    "RegionalAnalyzer",
    "identify_boundary_spots",
    "compare_tumor_stroma_ratios",
]

"""SPLISOSM: Spatial Isoform Statistical Omnibus Model.

A framework for detecting spatially variable isoform patterns in
spatial transcriptomics data, with application to cancer cytokine
and secreted protein analysis in the tumor microenvironment.

Main components:
- api: High-level SplisosmNP class for AnnData integration
- kernels: Spatial (ICAR) and compositional kernels
- statistics: HSIC-based test statistics and p-value computation
- analysis: Single-sample and multi-sample analysis pipelines
- io: Data loading for Visium and AnnData formats
- preprocessing: Imputation and filtering utilities
- gene_lists: Curated cytokine/RBP gene lists
"""

__version__ = "0.1.0"

# High-level API (recommended entry point)
from .api import SplisosmNP

# Core analysis
from .analysis.single_sample import (
    SplisomAnalyzer,
    GeneResult,
    SampleResult,
    build_isoform_mapping,
)

# Kernels
from .kernels import (
    ICARKernel,
    CompositionalKernel,
    GaussianKernel,
    BoundaryWeightedKernel,
    clr_transform,
)

# Statistics
from .statistics import (
    compute_hsic_gc,
    compute_hsic_ir,
    compute_hsic_ic,
    gamma_approximation_pvalue,
    hsic_pvalue,
)

# I/O
from .io import (
    load_visium_tsv,
    extract_isoform_counts,
    find_sample_files,
)

# Meta-analysis
from .analysis.meta import (
    meta_analyze_results,
    fisher_combine,
    stouffer_combine,
)

# Gene lists
from .gene_lists import (
    CYTOKINE_TARGETS,
    RBP_REGULATORS,
    HYPOXIA_MARKERS,
    STROMAL_MARKERS,
    IMMUNE_MARKERS,
    get_all_cytokines,
    get_category_genes,
    get_rbp_regulators,
    filter_genes_by_category,
)

# GPU utilities
from .utils.backend import (
    is_gpu_available,
    get_gpu_memory_info,
    clear_gpu_memory,
    CUPY_AVAILABLE,
)

__all__ = [
    # Version
    "__version__",
    # High-level API
    "SplisosmNP",
    # Core
    "SplisomAnalyzer",
    "GeneResult",
    "SampleResult",
    "build_isoform_mapping",
    # Kernels
    "ICARKernel",
    "CompositionalKernel",
    "GaussianKernel",
    "BoundaryWeightedKernel",
    "clr_transform",
    # Statistics
    "compute_hsic_gc",
    "compute_hsic_ir",
    "compute_hsic_ic",
    "gamma_approximation_pvalue",
    "hsic_pvalue",
    # I/O
    "load_visium_tsv",
    "extract_isoform_counts",
    "find_sample_files",
    # Meta-analysis
    "meta_analyze_results",
    "fisher_combine",
    "stouffer_combine",
    # Gene lists
    "CYTOKINE_TARGETS",
    "RBP_REGULATORS",
    "HYPOXIA_MARKERS",
    "STROMAL_MARKERS",
    "IMMUNE_MARKERS",
    "get_all_cytokines",
    "get_category_genes",
    "get_rbp_regulators",
    "filter_genes_by_category",
    # GPU utilities
    "is_gpu_available",
    "get_gpu_memory_info",
    "clear_gpu_memory",
    "CUPY_AVAILABLE",
]

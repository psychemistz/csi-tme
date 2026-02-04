"""Single-sample SPLISOSM analysis pipeline.

Main entry point for analyzing spatial isoform variation in a single sample.
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from ..kernels.icar import ICARKernel
from ..kernels.compositional import clr_transform
from ..statistics.hsic import (
    compute_hsic_gc,
    compute_hsic_ir,
    compute_hsic_ic,
    compute_hsic_combined,
    center_kernel,
    get_kernel_eigenvalues,
)
from ..statistics.pvalue import hsic_pvalue, PValueResult
from ..preprocessing.imputation import impute_with_counts
from ..preprocessing.filtering import filter_by_isoform_coverage


@dataclass
class GeneResult:
    """Result for a single gene's spatial isoform analysis.

    Attributes
    ----------
    gene : str
        Gene name
    statistic_ir : float
        HSIC-IR statistic (isoform ratio variation)
    pvalue_ir : float
        P-value for HSIC-IR
    statistic_gc : float
        HSIC-GC statistic (gene expression variation)
    pvalue_gc : float
        P-value for HSIC-GC
    statistic_ic : float
        HSIC-IC statistic (combined)
    pvalue_ic : float
        P-value for HSIC-IC
    n_spots_expressed : int
        Number of spots with gene expression
    n_isoforms : int
        Number of isoforms
    """
    gene: str
    statistic_ir: float
    pvalue_ir: float
    statistic_gc: float
    pvalue_gc: float
    statistic_ic: float
    pvalue_ic: float
    n_spots_expressed: int
    n_isoforms: int


@dataclass
class SampleResult:
    """Result for a full sample analysis.

    Attributes
    ----------
    sample_id : str
        Sample identifier
    n_spots : int
        Number of spatial spots
    n_genes_tested : int
        Number of genes tested
    gene_results : list
        List of GeneResult objects
    spatial_kernel_eigenvalues : ndarray
        Eigenvalues of centered spatial kernel (for diagnostics)
    """
    sample_id: str
    n_spots: int
    n_genes_tested: int
    gene_results: list = field(default_factory=list)
    spatial_kernel_eigenvalues: Optional[NDArray] = None

    def to_dataframe(self):
        """Convert results to pandas DataFrame."""
        import pandas as pd

        records = []
        for r in self.gene_results:
            records.append({
                'gene': r.gene,
                'statistic_ir': r.statistic_ir,
                'pvalue_ir': r.pvalue_ir,
                'statistic_gc': r.statistic_gc,
                'pvalue_gc': r.pvalue_gc,
                'statistic_ic': r.statistic_ic,
                'pvalue_ic': r.pvalue_ic,
                'n_spots_expressed': r.n_spots_expressed,
                'n_isoforms': r.n_isoforms,
            })

        df = pd.DataFrame(records)
        return df

    def get_significant_genes(
        self,
        test: str = 'ir',
        fdr_threshold: float = 0.05
    ) -> list:
        """Get genes with significant spatial isoform variation.

        Parameters
        ----------
        test : str
            Which test to use: 'ir', 'gc', or 'ic'
        fdr_threshold : float
            FDR threshold for significance

        Returns
        -------
        genes : list
            Significant gene names
        """
        from scipy import stats

        pvalue_key = f'pvalue_{test}'
        pvalues = [getattr(r, pvalue_key) for r in self.gene_results]

        # BH correction
        n = len(pvalues)
        sorted_idx = np.argsort(pvalues)
        sorted_pvalues = np.array(pvalues)[sorted_idx]

        thresholds = fdr_threshold * np.arange(1, n + 1) / n
        significant_mask = sorted_pvalues <= thresholds

        # Find largest k where p_k <= k * fdr / n
        if not np.any(significant_mask):
            return []

        k = np.max(np.where(significant_mask)[0]) + 1
        significant_idx = sorted_idx[:k]

        genes = [self.gene_results[i].gene for i in significant_idx]
        return genes


class SplisomAnalyzer:
    """Main analyzer class for SPLISOSM spatial isoform analysis.

    Parameters
    ----------
    k_neighbors : int
        Number of neighbors for ICAR kernel
    rho : float
        Autocorrelation parameter for ICAR
    min_spots_expressed : int
        Minimum spots with gene expression to test
    pseudocount : float
        Pseudocount for CLR transform
    use_gpu : bool
        If True, use GPU acceleration via CuPy (requires cupy package)
    """

    def __init__(
        self,
        k_neighbors: int = 20,
        rho: float = 0.99,
        min_spots_expressed: int = 50,
        pseudocount: float = 1e-10,
        use_gpu: bool = False
    ):
        self.k_neighbors = k_neighbors
        self.rho = rho
        self.min_spots_expressed = min_spots_expressed
        self.pseudocount = pseudocount
        self.use_gpu = use_gpu

        self._spatial_kernel: Optional[NDArray] = None
        self._spatial_eigenvalues: Optional[NDArray] = None

        # Check GPU availability if requested
        if use_gpu:
            from ..utils.backend import is_gpu_available
            if not is_gpu_available():
                import warnings
                warnings.warn(
                    "GPU requested but CuPy not available. Falling back to CPU.",
                    RuntimeWarning
                )
                self.use_gpu = False

    def build_spatial_kernel(self, coords: NDArray) -> NDArray:
        """Build and cache spatial ICAR kernel.

        Parameters
        ----------
        coords : ndarray
            Spatial coordinates (n_spots, 2)

        Returns
        -------
        L : ndarray
            Spatial kernel matrix
        """
        kernel = ICARKernel(
            k_neighbors=self.k_neighbors,
            rho=self.rho
        )
        L = kernel.compute(coords)
        self._spatial_kernel = L

        # Cache eigenvalues of centered kernel (GPU accelerated if enabled)
        L_centered = center_kernel(L, use_gpu=self.use_gpu)
        self._spatial_eigenvalues = get_kernel_eigenvalues(L_centered, use_gpu=self.use_gpu)

        return L

    def analyze_gene(
        self,
        isoform_counts: NDArray,
        gene_name: str
    ) -> Optional[GeneResult]:
        """Analyze spatial isoform variation for a single gene.

        Parameters
        ----------
        isoform_counts : ndarray
            Isoform counts (n_spots, n_isoforms)
        gene_name : str
            Gene name

        Returns
        -------
        result : GeneResult or None
            None if gene doesn't pass filters
        """
        if self._spatial_kernel is None:
            raise ValueError("Call build_spatial_kernel() first")

        n_spots, n_isoforms = isoform_counts.shape

        # Compute gene totals
        gene_totals = isoform_counts.sum(axis=1)
        n_expressed = (gene_totals > 0).sum()

        # Filter check
        if n_expressed < self.min_spots_expressed:
            return None

        # Skip single-isoform genes
        if n_isoforms < 2:
            return None

        # Impute ratios for spots with zero expression
        ratios = impute_with_counts(isoform_counts, gene_totals)

        # Center spatial kernel (cached computation via eigenvalues)
        L_centered = center_kernel(self._spatial_kernel, use_gpu=self.use_gpu)

        # HSIC-IR: Isoform ratio variation (GPU accelerated if enabled)
        T_ir, K_X = compute_hsic_ir(
            ratios, self._spatial_kernel,
            pseudocount=self.pseudocount,
            use_gpu=self.use_gpu
        )
        K_X_centered = center_kernel(K_X, use_gpu=self.use_gpu)
        eigenvalues_x = get_kernel_eigenvalues(K_X_centered, use_gpu=self.use_gpu)

        pval_ir = hsic_pvalue(
            T_ir, eigenvalues_x, self._spatial_eigenvalues, n_spots, method='liu'
        )

        # HSIC-GC: Gene expression variation (GPU accelerated if enabled)
        T_gc, K_g = compute_hsic_gc(gene_totals, self._spatial_kernel,
                                     use_gpu=self.use_gpu)
        K_g_centered = center_kernel(K_g, use_gpu=self.use_gpu)
        eigenvalues_g = get_kernel_eigenvalues(K_g_centered, use_gpu=self.use_gpu)

        pval_gc = hsic_pvalue(
            T_gc, eigenvalues_g, self._spatial_eigenvalues, n_spots, method='liu'
        )

        # HSIC-IC: Combined variation (GPU accelerated if enabled)
        T_ic, K_XG = compute_hsic_combined(
            ratios, gene_totals, self._spatial_kernel,
            pseudocount=self.pseudocount,
            use_gpu=self.use_gpu
        )
        K_XG_centered = center_kernel(K_XG, use_gpu=self.use_gpu)
        eigenvalues_xg = get_kernel_eigenvalues(K_XG_centered, use_gpu=self.use_gpu)

        pval_ic = hsic_pvalue(
            T_ic, eigenvalues_xg, self._spatial_eigenvalues, n_spots, method='liu'
        )

        return GeneResult(
            gene=gene_name,
            statistic_ir=T_ir,
            pvalue_ir=pval_ir.pvalue,
            statistic_gc=T_gc,
            pvalue_gc=pval_gc.pvalue,
            statistic_ic=T_ic,
            pvalue_ic=pval_ic.pvalue,
            n_spots_expressed=int(n_expressed),
            n_isoforms=n_isoforms
        )

    def analyze_sample(
        self,
        counts: NDArray,
        coords: NDArray,
        genes: list,
        isoform_mapping: dict,
        sample_id: str = "sample",
        show_progress: bool = True
    ) -> SampleResult:
        """Run full SPLISOSM analysis on a sample.

        Parameters
        ----------
        counts : ndarray
            Full count matrix (n_spots, n_genes)
        coords : ndarray
            Spatial coordinates (n_spots, 2)
        genes : list
            Gene/isoform names
        isoform_mapping : dict
            Maps gene name -> list of column indices for its isoforms
        sample_id : str
            Sample identifier
        show_progress : bool
            Show progress bar

        Returns
        -------
        result : SampleResult
        """
        n_spots = counts.shape[0]

        # Build spatial kernel
        self.build_spatial_kernel(coords)

        # Analyze each gene
        gene_results = []
        gene_iter = isoform_mapping.items()
        if show_progress:
            gene_iter = tqdm(gene_iter, desc=f"Analyzing {sample_id}")

        for gene_name, isoform_indices in gene_iter:
            isoform_counts = counts[:, isoform_indices]
            result = self.analyze_gene(isoform_counts, gene_name)
            if result is not None:
                gene_results.append(result)

        return SampleResult(
            sample_id=sample_id,
            n_spots=n_spots,
            n_genes_tested=len(gene_results),
            gene_results=gene_results,
            spatial_kernel_eigenvalues=self._spatial_eigenvalues
        )


def build_isoform_mapping(genes: list, separator: str = '-') -> dict:
    """Build mapping from base gene to isoform indices.

    Parameters
    ----------
    genes : list
        Gene/isoform names (e.g., ["VEGF-001", "VEGF-002", "IL6-001"])
    separator : str
        Character separating gene name from isoform suffix

    Returns
    -------
    mapping : dict
        Maps gene name -> list of column indices
    """
    mapping = {}
    for i, gene in enumerate(genes):
        parts = gene.rsplit(separator, 1)
        base_name = parts[0]

        if base_name not in mapping:
            mapping[base_name] = []
        mapping[base_name].append(i)

    # Filter to genes with multiple isoforms
    mapping = {k: v for k, v in mapping.items() if len(v) > 1}

    return mapping

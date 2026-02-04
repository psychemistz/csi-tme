"""High-level API for SPLISOSM analysis.

This module provides the SplisosmNP class, which matches the API documented
in RESEARCH.md for easy integration with AnnData workflows.
"""

from dataclasses import dataclass
from typing import Optional, List, Dict, Any, Union

import numpy as np
from numpy.typing import NDArray

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

try:
    import anndata as ad
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False

from .analysis.single_sample import SplisomAnalyzer, build_isoform_mapping, SampleResult
from .statistics.conditional import conditional_hsic, partial_correlation_spatial
from .statistics.pvalue import hsic_pvalue
from .statistics.hsic import center_kernel, get_kernel_eigenvalues
from .kernels.icar import ICARKernel
from .io.anndata import to_splisosm_format


def _check_pandas():
    """Check if pandas is available."""
    if not PANDAS_AVAILABLE:
        raise ImportError(
            "pandas is required for DataFrame output. "
            "Install with: pip install pandas"
        )


def _check_anndata():
    """Check if AnnData is available."""
    if not ANNDATA_AVAILABLE:
        raise ImportError(
            "anndata is required for SplisosmNP. "
            "Install with: pip install anndata"
        )


def benjamini_hochberg(pvalues: NDArray, alpha: float = 0.05) -> NDArray:
    """Apply Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    pvalues : ndarray
        Array of p-values
    alpha : float
        FDR threshold (not used in adjustment, just for reference)

    Returns
    -------
    padj : ndarray
        Adjusted p-values
    """
    n = len(pvalues)
    if n == 0:
        return np.array([])

    # Sort p-values and get order
    sorted_idx = np.argsort(pvalues)
    sorted_pvalues = pvalues[sorted_idx]

    # Compute adjusted p-values
    adjusted = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if i == n - 1:
            adjusted[i] = sorted_pvalues[i]
        else:
            adjusted[i] = min(adjusted[i + 1], sorted_pvalues[i] * n / (i + 1))

    # Restore original order
    padj = np.zeros(n)
    padj[sorted_idx] = adjusted

    # Clamp to [0, 1]
    padj = np.clip(padj, 0, 1)

    return padj


class SplisosmNP:
    """High-level SPLISOSM API for AnnData integration.

    This class provides a user-friendly interface matching the API documented
    in RESEARCH.md. It wraps the lower-level SplisomAnalyzer for easy
    integration with scanpy/AnnData workflows.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing spatial transcriptomics data.
        Should have isoform-level counts (genes named like "GENE-001", "GENE-002")
        or gene-level counts for gene family analysis.
    spatial_key : str, optional
        Key in adata.obsm containing spatial coordinates (default: 'spatial')
    k_neighbors : int, optional
        Number of neighbors for ICAR spatial kernel (default: 20)
    rho : float, optional
        Autocorrelation parameter for ICAR (default: 0.99)
    min_spots_expressed : int, optional
        Minimum spots with expression to test a gene (default: 50)
    use_gpu : bool, optional
        If True, use GPU acceleration via CuPy (default: False)
    isoform_separator : str, optional
        Character separating gene name from isoform suffix (default: '-')

    Examples
    --------
    >>> import scanpy as sc
    >>> from splisosm import SplisosmNP
    >>>
    >>> adata = sc.read_h5ad('spatial_data.h5ad')
    >>> model = SplisosmNP(adata, spatial_key='spatial')
    >>> results = model.test_spatial_variability(test_type='HSIC-IR')
    >>> sig_genes = results[results['padj_ir'] < 0.05]['gene'].tolist()
    """

    def __init__(
        self,
        adata: "ad.AnnData",
        spatial_key: str = 'spatial',
        k_neighbors: int = 20,
        rho: float = 0.99,
        min_spots_expressed: int = 50,
        use_gpu: bool = False,
        isoform_separator: str = '-'
    ):
        _check_anndata()
        _check_pandas()

        self.adata = adata
        self.spatial_key = spatial_key
        self.k_neighbors = k_neighbors
        self.rho = rho
        self.min_spots_expressed = min_spots_expressed
        self.use_gpu = use_gpu
        self.isoform_separator = isoform_separator

        # Extract data from AnnData
        self._counts, self._coords, self._genes = to_splisosm_format(
            adata, spatial_key=spatial_key
        )

        # Build isoform mapping
        self._isoform_mapping = build_isoform_mapping(
            self._genes, separator=isoform_separator
        )

        # Initialize analyzer
        self._analyzer = SplisomAnalyzer(
            k_neighbors=k_neighbors,
            rho=rho,
            min_spots_expressed=min_spots_expressed,
            use_gpu=use_gpu
        )

        # Cache for spatial kernel
        self._spatial_kernel = None
        self._spatial_kernel_eigenvalues = None

    @property
    def n_spots(self) -> int:
        """Number of spatial spots."""
        return self._counts.shape[0]

    @property
    def n_genes(self) -> int:
        """Number of genes/isoforms."""
        return self._counts.shape[1]

    @property
    def n_testable_genes(self) -> int:
        """Number of genes with multiple isoforms (testable)."""
        return len(self._isoform_mapping)

    @property
    def genes_with_isoforms(self) -> List[str]:
        """List of genes with multiple isoforms."""
        return list(self._isoform_mapping.keys())

    def _build_spatial_kernel(self) -> None:
        """Build and cache the spatial kernel."""
        if self._spatial_kernel is None:
            kernel = ICARKernel(k_neighbors=self.k_neighbors, rho=self.rho)
            self._spatial_kernel = kernel.compute(self._coords)

            # Cache eigenvalues for p-value computation
            L_centered = center_kernel(self._spatial_kernel, use_gpu=self.use_gpu)
            self._spatial_kernel_eigenvalues = get_kernel_eigenvalues(
                L_centered, use_gpu=self.use_gpu
            )

    def test_spatial_variability(
        self,
        test_type: str = 'HSIC-IR',
        genes: Optional[List[str]] = None,
        kernel_type: str = 'icar',
        rho: Optional[float] = None,
        k: Optional[int] = None,
        fdr_threshold: float = 0.05,
        show_progress: bool = True
    ) -> "pd.DataFrame":
        """Test for spatially variable isoform patterns.

        This is the primary analysis function that tests whether isoform
        usage varies spatially across the tissue.

        Parameters
        ----------
        test_type : str
            Type of test to run:
            - 'HSIC-IR': Isoform ratio variation (primary, default)
            - 'HSIC-GC': Gene expression variation
            - 'HSIC-IC': Combined expression + composition
            - 'all': Run all three tests
        genes : list, optional
            Specific genes to test. If None, test all genes with isoforms.
        kernel_type : str
            Spatial kernel type: 'icar' (default) or 'gaussian'
        rho : float, optional
            Override ICAR autocorrelation parameter
        k : int, optional
            Override number of neighbors
        fdr_threshold : float
            FDR threshold for significance (used for adjusted p-values)
        show_progress : bool
            Show progress bar

        Returns
        -------
        results : DataFrame
            DataFrame with columns:
            - gene: Gene name
            - statistic_ir, pvalue_ir, padj_ir: HSIC-IR results
            - statistic_gc, pvalue_gc, padj_gc: HSIC-GC results
            - statistic_ic, pvalue_ic, padj_ic: HSIC-IC results
            - n_spots_expressed: Number of spots with expression
            - n_isoforms: Number of isoforms
        """
        _check_pandas()

        # Override parameters if provided
        if rho is not None:
            self._analyzer.rho = rho
            self._spatial_kernel = None  # Force rebuild
        if k is not None:
            self._analyzer.k_neighbors = k
            self._spatial_kernel = None  # Force rebuild

        # Filter genes if specified
        if genes is not None:
            mapping = {g: v for g, v in self._isoform_mapping.items() if g in genes}
        else:
            mapping = self._isoform_mapping

        if len(mapping) == 0:
            return pd.DataFrame()

        # Run analysis
        result = self._analyzer.analyze_sample(
            counts=self._counts,
            coords=self._coords,
            genes=self._genes,
            isoform_mapping=mapping,
            sample_id="sample",
            show_progress=show_progress
        )

        # Convert to DataFrame
        df = result.to_dataframe()

        if len(df) == 0:
            return df

        # Apply BH correction
        df['padj_ir'] = benjamini_hochberg(df['pvalue_ir'].values)
        df['padj_gc'] = benjamini_hochberg(df['pvalue_gc'].values)
        df['padj_ic'] = benjamini_hochberg(df['pvalue_ic'].values)

        # Sort by primary test p-value
        if test_type == 'HSIC-IR' or test_type == 'all':
            df = df.sort_values('pvalue_ir')
        elif test_type == 'HSIC-GC':
            df = df.sort_values('pvalue_gc')
        elif test_type == 'HSIC-IC':
            df = df.sort_values('pvalue_ic')

        return df.reset_index(drop=True)

    def test_differential_usage(
        self,
        genes: List[str],
        covariates: Union[List[str], NDArray],
        conditional: bool = True,
        method: str = 'hsic',
        rbf_smoothing: float = 1.0,
        fdr_threshold: float = 0.05
    ) -> "pd.DataFrame":
        """Test differential isoform usage with RBP or other covariates.

        Tests whether isoform usage correlates with covariate expression
        (e.g., RNA-binding proteins) while optionally controlling for
        spatial confounding.

        Parameters
        ----------
        genes : list
            Genes to test (must have multiple isoforms)
        covariates : list or ndarray
            Either:
            - List of covariate gene names (looked up in adata)
            - ndarray of shape (n_spots,) or (n_spots, n_covariates)
        conditional : bool
            If True, use conditional HSIC to control for spatial confounding
        method : str
            Test method: 'hsic' (default) or 'correlation'
        rbf_smoothing : float
            Smoothing parameter for spatial residualization (conditional HSIC)
        fdr_threshold : float
            FDR threshold for adjusted p-values

        Returns
        -------
        results : DataFrame
            DataFrame with columns:
            - gene: Gene name
            - covariate: Covariate name or index
            - statistic: Test statistic
            - pvalue: Raw p-value
            - padj: BH-adjusted p-value
        """
        _check_pandas()

        self._build_spatial_kernel()

        results = []

        # Get covariate matrix
        if isinstance(covariates, list):
            covariate_names = covariates
            covariate_matrix = self._get_covariate_matrix(covariates)
        else:
            covariate_matrix = np.atleast_2d(covariates)
            if covariate_matrix.shape[0] == self.n_spots:
                covariate_matrix = covariate_matrix
            else:
                covariate_matrix = covariate_matrix.T
            covariate_names = [f'cov_{i}' for i in range(covariate_matrix.shape[1])]

        # Test each gene-covariate pair
        for gene in genes:
            if gene not in self._isoform_mapping:
                continue

            # Get isoform proportions for this gene
            isoform_indices = self._isoform_mapping[gene]
            isoform_counts = self._counts[:, isoform_indices]
            gene_totals = isoform_counts.sum(axis=1)

            # Skip genes with insufficient expression
            n_expressed = (gene_totals > 0).sum()
            if n_expressed < self.min_spots_expressed:
                continue

            # Compute isoform proportions with mean imputation
            from .preprocessing.imputation import impute_with_counts
            proportions = impute_with_counts(isoform_counts, gene_totals)

            # Test against each covariate
            for i, cov_name in enumerate(covariate_names):
                cov_values = covariate_matrix[:, i]

                if conditional:
                    # Conditional HSIC controlling for spatial location
                    result = conditional_hsic(
                        proportions, cov_values, self._coords,
                        method='rbf', rbf_smoothing=rbf_smoothing
                    )
                    stat = result.statistic
                    pval = result.pvalue
                else:
                    # Standard correlation/HSIC
                    if method == 'correlation':
                        # Average isoform proportion change
                        from scipy import stats as scipy_stats
                        # Use first principal component of proportions
                        prop_centered = proportions - proportions.mean(axis=0)
                        if prop_centered.shape[1] > 1:
                            u, s, vh = np.linalg.svd(prop_centered, full_matrices=False)
                            prop_pc1 = u[:, 0] * s[0]
                        else:
                            prop_pc1 = prop_centered[:, 0]
                        stat, pval = scipy_stats.pearsonr(prop_pc1, cov_values)
                    else:
                        # HSIC between proportions and covariate
                        from .statistics.hsic import compute_hsic_ir
                        # Build covariate kernel
                        cov_centered = cov_values - cov_values.mean()
                        K_cov = np.outer(cov_centered, cov_centered)

                        # Build isoform kernel
                        X_centered = proportions - proportions.mean(axis=0)
                        K_iso = X_centered @ X_centered.T

                        # Compute HSIC
                        from .statistics.hsic import compute_hsic_statistic
                        stat = compute_hsic_statistic(K_iso, K_cov)

                        # Get p-value
                        K_iso_centered = center_kernel(K_iso)
                        K_cov_centered = center_kernel(K_cov)
                        eig_iso = get_kernel_eigenvalues(K_iso_centered)
                        eig_cov = get_kernel_eigenvalues(K_cov_centered)
                        pval_result = hsic_pvalue(
                            stat, eig_iso, eig_cov, self.n_spots, method='liu'
                        )
                        pval = pval_result.pvalue

                results.append({
                    'gene': gene,
                    'covariate': cov_name,
                    'statistic': stat,
                    'pvalue': pval
                })

        if len(results) == 0:
            return pd.DataFrame(columns=['gene', 'covariate', 'statistic', 'pvalue', 'padj'])

        df = pd.DataFrame(results)
        df['padj'] = benjamini_hochberg(df['pvalue'].values)
        df = df.sort_values('pvalue')

        return df.reset_index(drop=True)

    def _get_covariate_matrix(self, covariate_names: List[str]) -> NDArray:
        """Extract covariate expression matrix from adata."""
        matrices = []

        for name in covariate_names:
            if name in self._genes:
                idx = self._genes.index(name)
                matrices.append(self._counts[:, idx])
            elif name in self.adata.var_names:
                idx = list(self.adata.var_names).index(name)
                if hasattr(self.adata.X, 'toarray'):
                    matrices.append(self.adata.X[:, idx].toarray().flatten())
                else:
                    matrices.append(np.array(self.adata.X[:, idx]).flatten())
            else:
                raise ValueError(f"Covariate '{name}' not found in data")

        return np.column_stack(matrices)

    def get_significant_genes(
        self,
        results: "pd.DataFrame",
        test: str = 'ir',
        fdr_threshold: float = 0.05
    ) -> List[str]:
        """Get genes with significant spatial isoform variation.

        Parameters
        ----------
        results : DataFrame
            Results from test_spatial_variability()
        test : str
            Which test to use: 'ir', 'gc', or 'ic'
        fdr_threshold : float
            FDR threshold for significance

        Returns
        -------
        genes : list
            Significant gene names
        """
        padj_col = f'padj_{test}'
        if padj_col not in results.columns:
            raise ValueError(f"Column '{padj_col}' not found. Run test_spatial_variability() first.")

        sig_mask = results[padj_col] < fdr_threshold
        return results.loc[sig_mask, 'gene'].tolist()

    def filter_by_category(
        self,
        results: "pd.DataFrame",
        category: str
    ) -> "pd.DataFrame":
        """Filter results to specific gene category.

        Parameters
        ----------
        results : DataFrame
            Results from test_spatial_variability()
        category : str
            Category from gene_lists (e.g., 'interleukins', 'growth_factors')

        Returns
        -------
        filtered : DataFrame
            Filtered results
        """
        from .gene_lists import CYTOKINE_TARGETS

        if category not in CYTOKINE_TARGETS:
            raise ValueError(f"Unknown category: {category}. "
                           f"Available: {list(CYTOKINE_TARGETS.keys())}")

        genes_in_category = CYTOKINE_TARGETS[category]
        mask = results['gene'].isin(genes_in_category)
        return results[mask].reset_index(drop=True)

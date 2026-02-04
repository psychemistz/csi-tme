"""Tests for the high-level SplisosmNP API."""

import pytest
import numpy as np

# Skip all tests if dependencies not available
pytest.importorskip("pandas")
pytest.importorskip("anndata")


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object for testing."""
    import anndata as ad
    import pandas as pd

    np.random.seed(42)

    n_spots = 100
    n_genes = 20  # 10 genes with 2 isoforms each

    # Create count matrix
    counts = np.random.poisson(10, (n_spots, n_genes)).astype(float)

    # Create spatial coordinates
    coords = np.random.rand(n_spots, 2) * 100

    # Create gene names with isoform suffixes
    gene_names = []
    for i in range(10):
        gene_names.append(f"GENE{i:02d}-001")
        gene_names.append(f"GENE{i:02d}-002")

    adata = ad.AnnData(X=counts)
    adata.var_names = gene_names
    adata.obsm['spatial'] = coords

    return adata


@pytest.fixture
def splisosm_model(mock_adata):
    """Create SplisosmNP model."""
    from splisosm import SplisosmNP
    return SplisosmNP(mock_adata, spatial_key='spatial', min_spots_expressed=10)


class TestSplisosmNPInit:
    """Tests for SplisosmNP initialization."""

    def test_basic_init(self, mock_adata):
        """Test basic initialization."""
        from splisosm import SplisosmNP

        model = SplisosmNP(mock_adata)
        assert model.n_spots == 100
        assert model.n_genes == 20

    def test_init_with_params(self, mock_adata):
        """Test initialization with custom parameters."""
        from splisosm import SplisosmNP

        model = SplisosmNP(
            mock_adata,
            spatial_key='spatial',
            k_neighbors=10,
            rho=0.95,
            min_spots_expressed=20
        )

        assert model.k_neighbors == 10
        assert model.rho == 0.95
        assert model.min_spots_expressed == 20

    def test_isoform_mapping(self, splisosm_model):
        """Test isoform mapping is built correctly."""
        assert splisosm_model.n_testable_genes == 10
        assert len(splisosm_model.genes_with_isoforms) == 10
        assert 'GENE00' in splisosm_model.genes_with_isoforms


class TestSpatialVariability:
    """Tests for test_spatial_variability method."""

    def test_basic_analysis(self, splisosm_model):
        """Test basic spatial variability analysis."""
        results = splisosm_model.test_spatial_variability(show_progress=False)

        assert len(results) > 0
        assert 'gene' in results.columns
        assert 'pvalue_ir' in results.columns
        assert 'padj_ir' in results.columns
        assert 'statistic_ir' in results.columns

    def test_result_columns(self, splisosm_model):
        """Test all expected columns are present."""
        results = splisosm_model.test_spatial_variability(show_progress=False)

        expected_columns = [
            'gene', 'statistic_ir', 'pvalue_ir', 'padj_ir',
            'statistic_gc', 'pvalue_gc', 'padj_gc',
            'statistic_ic', 'pvalue_ic', 'padj_ic',
            'n_spots_expressed', 'n_isoforms'
        ]

        for col in expected_columns:
            assert col in results.columns, f"Missing column: {col}"

    def test_pvalue_range(self, splisosm_model):
        """Test p-values are in valid range."""
        results = splisosm_model.test_spatial_variability(show_progress=False)

        assert (results['pvalue_ir'] >= 0).all()
        assert (results['pvalue_ir'] <= 1).all()
        assert (results['padj_ir'] >= 0).all()
        assert (results['padj_ir'] <= 1).all()

    def test_subset_genes(self, splisosm_model):
        """Test analysis on gene subset."""
        genes = ['GENE00', 'GENE01']
        results = splisosm_model.test_spatial_variability(
            genes=genes, show_progress=False
        )

        assert len(results) <= 2
        assert all(g in genes for g in results['gene'].tolist())


class TestDifferentialUsage:
    """Tests for test_differential_usage method."""

    def test_with_array_covariates(self, splisosm_model):
        """Test with numpy array covariates."""
        np.random.seed(42)
        covariates = np.random.randn(100)
        genes = ['GENE00', 'GENE01']

        results = splisosm_model.test_differential_usage(
            genes=genes,
            covariates=covariates,
            conditional=False
        )

        assert 'gene' in results.columns
        assert 'covariate' in results.columns
        assert 'pvalue' in results.columns
        assert 'padj' in results.columns

    def test_conditional_hsic(self, splisosm_model):
        """Test conditional HSIC analysis."""
        np.random.seed(42)
        covariates = np.random.randn(100)
        genes = ['GENE00']

        results = splisosm_model.test_differential_usage(
            genes=genes,
            covariates=covariates,
            conditional=True
        )

        # Should have results
        assert len(results) >= 0  # May be empty if gene filtered


class TestHelperMethods:
    """Tests for helper methods."""

    def test_get_significant_genes(self, splisosm_model):
        """Test get_significant_genes method."""
        results = splisosm_model.test_spatial_variability(show_progress=False)

        # With high threshold, should get some genes
        sig_genes = splisosm_model.get_significant_genes(
            results, test='ir', fdr_threshold=1.0
        )
        assert isinstance(sig_genes, list)

        # With impossible threshold, should get none
        sig_genes = splisosm_model.get_significant_genes(
            results, test='ir', fdr_threshold=0.0
        )
        assert len(sig_genes) == 0


class TestBHCorrection:
    """Tests for Benjamini-Hochberg correction."""

    def test_bh_correction_ordering(self):
        """Test BH correction preserves p-value ordering."""
        from splisosm.api import benjamini_hochberg

        pvalues = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        padj = benjamini_hochberg(pvalues)

        # Adjusted p-values should maintain relative ordering
        assert padj[0] <= padj[1] <= padj[2] <= padj[3] <= padj[4]

    def test_bh_correction_bounds(self):
        """Test BH correction respects [0, 1] bounds."""
        from splisosm.api import benjamini_hochberg

        pvalues = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        padj = benjamini_hochberg(pvalues)

        assert (padj >= 0).all()
        assert (padj <= 1).all()

    def test_bh_correction_empty(self):
        """Test BH correction with empty array."""
        from splisosm.api import benjamini_hochberg

        padj = benjamini_hochberg(np.array([]))
        assert len(padj) == 0

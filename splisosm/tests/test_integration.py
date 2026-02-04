"""Integration tests for full SPLISOSM pipeline."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from splisosm import (
    SplisomAnalyzer,
    build_isoform_mapping,
    GeneResult,
    SampleResult,
)
from splisosm.io import load_visium_tsv, extract_isoform_counts
from splisosm.preprocessing import impute_with_counts
from splisosm.analysis.meta import fisher_combine, stouffer_combine, meta_analyze_results


class TestSplisomAnalyzer:
    @pytest.fixture
    def simulated_data(self):
        """Generate simulated spatial transcriptomics data."""
        np.random.seed(42)

        # 15x15 grid
        n_spots = 225
        x = np.arange(15)
        y = np.arange(15)
        xx, yy = np.meshgrid(x, y)
        coords = np.column_stack([xx.ravel(), yy.ravel()])

        # 3 genes with 2-3 isoforms each
        genes = ['GENE1-001', 'GENE1-002',
                 'GENE2-001', 'GENE2-002', 'GENE2-003',
                 'GENE3-001', 'GENE3-002']

        # Generate counts
        n_genes = len(genes)
        counts = np.random.poisson(10, (n_spots, n_genes)).astype(float)

        # Make GENE1 have spatially varying isoform ratios
        # (left half has more isoform 1, right half has more isoform 2)
        left_mask = coords[:, 0] < 7
        counts[left_mask, 0] *= 3   # GENE1-001 high on left
        counts[~left_mask, 1] *= 3  # GENE1-002 high on right

        return counts, coords, genes

    def test_analyzer_initialization(self):
        analyzer = SplisomAnalyzer(
            k_neighbors=10,
            rho=0.95,
            min_spots_expressed=30
        )
        assert analyzer.k_neighbors == 10
        assert analyzer.rho == 0.95

    def test_build_spatial_kernel(self, simulated_data):
        counts, coords, genes = simulated_data

        analyzer = SplisomAnalyzer(k_neighbors=10)
        L = analyzer.build_spatial_kernel(coords)

        n = len(coords)
        assert L.shape == (n, n)
        assert analyzer._spatial_kernel is not None
        assert analyzer._spatial_eigenvalues is not None

    def test_analyze_gene(self, simulated_data):
        counts, coords, genes = simulated_data

        analyzer = SplisomAnalyzer(k_neighbors=10, min_spots_expressed=20)
        analyzer.build_spatial_kernel(coords)

        # Get GENE1 isoform counts
        isoform_counts = counts[:, :2]  # GENE1-001, GENE1-002

        result = analyzer.analyze_gene(isoform_counts, 'GENE1')

        assert isinstance(result, GeneResult)
        assert result.gene == 'GENE1'
        assert 0 <= result.pvalue_ir <= 1
        assert 0 <= result.pvalue_gc <= 1
        assert 0 <= result.pvalue_ic <= 1
        assert result.n_isoforms == 2

    def test_spatially_varying_detected(self, simulated_data):
        """Gene with spatial isoform pattern should have higher statistic."""
        counts, coords, genes = simulated_data

        analyzer = SplisomAnalyzer(k_neighbors=10, min_spots_expressed=20)
        analyzer.build_spatial_kernel(coords)

        # GENE1 has spatial pattern
        result_gene1 = analyzer.analyze_gene(counts[:, :2], 'GENE1')

        # GENE2 and GENE3 are random
        result_gene2 = analyzer.analyze_gene(counts[:, 2:5], 'GENE2')

        # GENE1 should have higher HSIC-IR statistic (stronger spatial signal)
        assert result_gene1.statistic_ir > result_gene2.statistic_ir

    def test_analyze_sample(self, simulated_data):
        counts, coords, genes = simulated_data

        isoform_mapping = build_isoform_mapping(genes, separator='-')

        analyzer = SplisomAnalyzer(k_neighbors=10, min_spots_expressed=20)
        result = analyzer.analyze_sample(
            counts, coords, genes, isoform_mapping,
            sample_id='test_sample',
            show_progress=False
        )

        assert isinstance(result, SampleResult)
        assert result.sample_id == 'test_sample'
        assert result.n_spots == 225
        assert result.n_genes_tested > 0
        assert len(result.gene_results) > 0

    def test_to_dataframe(self, simulated_data):
        counts, coords, genes = simulated_data

        isoform_mapping = build_isoform_mapping(genes, separator='-')

        analyzer = SplisomAnalyzer(k_neighbors=10, min_spots_expressed=20)
        result = analyzer.analyze_sample(
            counts, coords, genes, isoform_mapping,
            sample_id='test',
            show_progress=False
        )

        df = result.to_dataframe()

        assert 'gene' in df.columns
        assert 'pvalue_ir' in df.columns
        assert 'statistic_ir' in df.columns
        assert len(df) == result.n_genes_tested


class TestBuildIsoformMapping:
    def test_simple_mapping(self):
        genes = ['VEGF-001', 'VEGF-002', 'IL6-001', 'IL6-002', 'IL6-003']
        mapping = build_isoform_mapping(genes, separator='-')

        assert 'VEGF' in mapping
        assert 'IL6' in mapping
        assert mapping['VEGF'] == [0, 1]
        assert mapping['IL6'] == [2, 3, 4]

    def test_single_isoform_excluded(self):
        genes = ['VEGF-001', 'VEGF-002', 'SOLO-001', 'IL6-001', 'IL6-002']
        mapping = build_isoform_mapping(genes, separator='-')

        # SOLO should be excluded (only one isoform)
        assert 'SOLO' not in mapping
        assert 'VEGF' in mapping
        assert 'IL6' in mapping


class TestMetaAnalysis:
    def test_fisher_combine(self):
        pvalues = np.array([0.01, 0.05, 0.1, 0.5])
        combined, chi2 = fisher_combine(pvalues)

        assert 0 <= combined <= 1
        assert chi2 > 0

    def test_stouffer_combine(self):
        pvalues = np.array([0.01, 0.05, 0.1, 0.5])
        combined, z = stouffer_combine(pvalues)

        assert 0 <= combined <= 1

    def test_fisher_extreme_pvalues(self):
        # Very small p-values
        pvalues = np.array([1e-10, 1e-15, 1e-20])
        combined, _ = fisher_combine(pvalues)

        assert combined < 1e-20

    def test_stouffer_weighted(self):
        pvalues = np.array([0.01, 0.5])
        weights = np.array([100, 10])  # First sample much larger

        combined_weighted, _ = stouffer_combine(pvalues, weights)
        combined_unweighted, _ = stouffer_combine(pvalues)

        # Both should produce valid p-values
        assert 0 <= combined_weighted <= 1
        assert 0 <= combined_unweighted <= 1


class TestPreprocessing:
    def test_impute_with_counts(self):
        # Some spots have zero counts
        isoform_counts = np.array([
            [10, 20, 30],
            [0, 0, 0],    # Zero total
            [5, 10, 15],
            [0, 0, 0],    # Zero total
            [8, 12, 20],
        ])

        ratios = impute_with_counts(isoform_counts)

        # All rows should sum to 1
        assert_allclose(ratios.sum(axis=1), 1.0)

        # No NaN values
        assert not np.any(np.isnan(ratios))

        # Zero-total spots should have imputed values
        assert ratios[1, 0] > 0
        assert ratios[3, 0] > 0


class TestVisiumIO:
    def test_parse_coordinate(self):
        from splisosm.io.visium import parse_coordinate

        row, col = parse_coordinate('5x10')
        assert row == 5
        assert col == 10

        row, col = parse_coordinate('1x1')
        assert row == 1
        assert col == 1

    def test_parse_invalid_coordinate(self):
        from splisosm.io.visium import parse_coordinate

        with pytest.raises(ValueError):
            parse_coordinate('invalid')

    def test_get_study_info(self):
        from splisosm.io.visium import get_study_info

        info = get_study_info('/data/BRCA_2021_Wu/sample123_counts.tsv')

        assert info['study'] == 'BRCA_2021_Wu'
        assert info['sample_id'] == 'sample123'
        assert info['cancer_type'] == 'BRCA'


class TestEndToEnd:
    """Full pipeline test with synthetic data."""

    def test_complete_pipeline(self):
        np.random.seed(123)

        # Generate two "samples"
        samples = []
        for sample_idx in range(2):
            # 10x10 grid
            n_spots = 100
            x = np.arange(10)
            y = np.arange(10)
            xx, yy = np.meshgrid(x, y)
            coords = np.column_stack([xx.ravel(), yy.ravel()])

            # Gene with spatial isoform pattern
            genes = ['TESTGENE-001', 'TESTGENE-002']
            counts = np.random.poisson(20, (n_spots, 2)).astype(float)

            # Add spatial pattern
            left = coords[:, 0] < 5
            counts[left, 0] *= 2
            counts[~left, 1] *= 2

            samples.append({
                'counts': counts,
                'coords': coords,
                'genes': genes,
            })

        # Analyze each sample
        sample_results = []
        for i, sample in enumerate(samples):
            isoform_mapping = build_isoform_mapping(sample['genes'], '-')

            analyzer = SplisomAnalyzer(k_neighbors=8, min_spots_expressed=30)
            result = analyzer.analyze_sample(
                sample['counts'],
                sample['coords'],
                sample['genes'],
                isoform_mapping,
                sample_id=f'sample_{i}',
                show_progress=False
            )
            sample_results.append(result)

        # Meta-analysis
        meta_results = meta_analyze_results(sample_results, test='ir')

        assert 'TESTGENE' in meta_results
        assert meta_results['TESTGENE'].n_samples == 2
        # Just verify meta-analysis ran and produced valid p-value
        assert 0 <= meta_results['TESTGENE'].combined_pvalue <= 1

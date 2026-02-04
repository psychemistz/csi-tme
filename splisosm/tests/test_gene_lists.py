"""Tests for gene_lists module."""

import pytest


class TestCytokineTargets:
    """Tests for CYTOKINE_TARGETS."""

    def test_cytokine_categories_exist(self):
        """Test all expected categories exist."""
        from splisosm.gene_lists import CYTOKINE_TARGETS

        expected_categories = [
            'interleukins',
            'growth_factors',
            'chemokines',
            'death_ligands',
            'cytokine_inhibitors',
            'proteases',
            'convertases',
        ]

        for category in expected_categories:
            assert category in CYTOKINE_TARGETS, f"Missing category: {category}"

    def test_cytokine_lists_not_empty(self):
        """Test each category has genes."""
        from splisosm.gene_lists import CYTOKINE_TARGETS

        for category, genes in CYTOKINE_TARGETS.items():
            assert len(genes) > 0, f"Empty gene list for {category}"

    def test_key_cytokines_present(self):
        """Test key cytokines are in the lists."""
        from splisosm.gene_lists import CYTOKINE_TARGETS

        # Check specific important genes
        assert 'IL33' in CYTOKINE_TARGETS['interleukins']
        assert 'VEGFA' in CYTOKINE_TARGETS['growth_factors']
        assert 'CXCL12' in CYTOKINE_TARGETS['chemokines']
        assert 'TNFSF10' in CYTOKINE_TARGETS['death_ligands']


class TestRBPRegulators:
    """Tests for RBP_REGULATORS."""

    def test_rbp_list_not_empty(self):
        """Test RBP list has genes."""
        from splisosm.gene_lists import RBP_REGULATORS

        assert len(RBP_REGULATORS) > 0

    def test_key_rbps_present(self):
        """Test key RBPs are in the list."""
        from splisosm.gene_lists import RBP_REGULATORS

        key_rbps = ['RBFOX1', 'CELF1', 'QKI', 'PTBP1', 'SRSF1']
        for rbp in key_rbps:
            assert rbp in RBP_REGULATORS, f"Missing RBP: {rbp}"


class TestHelperFunctions:
    """Tests for helper functions."""

    def test_get_all_cytokines(self):
        """Test get_all_cytokines returns flat list."""
        from splisosm.gene_lists import get_all_cytokines, CYTOKINE_TARGETS

        all_cytokines = get_all_cytokines()

        # Should be a flat list
        assert isinstance(all_cytokines, list)
        assert all(isinstance(g, str) for g in all_cytokines)

        # Should contain genes from all categories
        expected_count = sum(len(genes) for genes in CYTOKINE_TARGETS.values())
        # Note: get_all_cytokines removes duplicates, so may be slightly less
        assert len(all_cytokines) <= expected_count
        assert len(all_cytokines) > 0

    def test_get_category_genes(self):
        """Test get_category_genes returns correct genes."""
        from splisosm.gene_lists import get_category_genes, CYTOKINE_TARGETS

        for category in CYTOKINE_TARGETS:
            genes = get_category_genes(category)
            assert genes == CYTOKINE_TARGETS[category]

    def test_get_category_genes_invalid(self):
        """Test get_category_genes raises on invalid category."""
        from splisosm.gene_lists import get_category_genes

        with pytest.raises(ValueError):
            get_category_genes('invalid_category')

    def test_filter_genes_by_category(self):
        """Test filter_genes_by_category."""
        from splisosm.gene_lists import filter_genes_by_category

        genes = ['IL33', 'VEGFA', 'NOTACYTOKINE', 'CXCL12']

        # Filter to interleukins
        filtered = filter_genes_by_category(genes, category='interleukins')
        assert filtered == ['IL33']

        # Filter to growth factors
        filtered = filter_genes_by_category(genes, category='growth_factors')
        assert filtered == ['VEGFA']

    def test_filter_genes_all_cytokines(self):
        """Test filter_genes_by_category with all cytokines."""
        from splisosm.gene_lists import filter_genes_by_category

        genes = ['IL33', 'VEGFA', 'NOTACYTOKINE', 'CXCL12']

        filtered = filter_genes_by_category(genes, include_all_cytokines=True)
        assert 'IL33' in filtered
        assert 'VEGFA' in filtered
        assert 'CXCL12' in filtered
        assert 'NOTACYTOKINE' not in filtered

    def test_get_rbp_regulators(self):
        """Test get_rbp_regulators returns copy."""
        from splisosm.gene_lists import get_rbp_regulators, RBP_REGULATORS

        rbps = get_rbp_regulators()
        assert rbps == RBP_REGULATORS

        # Should be a copy, not the same object
        rbps.append('TEST')
        assert 'TEST' not in RBP_REGULATORS

    def test_get_hypoxia_markers(self):
        """Test get_hypoxia_markers."""
        from splisosm.gene_lists import get_hypoxia_markers

        markers = get_hypoxia_markers()
        assert 'HIF1A' in markers
        assert 'LDHA' in markers
        assert 'CA9' in markers

    def test_get_stromal_markers(self):
        """Test get_stromal_markers."""
        from splisosm.gene_lists import get_stromal_markers

        markers = get_stromal_markers()
        assert 'ACTA2' in markers
        assert 'FAP' in markers

    def test_get_immune_markers_all(self):
        """Test get_immune_markers returns all markers."""
        from splisosm.gene_lists import get_immune_markers

        markers = get_immune_markers()
        assert 'CD3D' in markers  # T cells
        assert 'CD68' in markers  # Macrophages

    def test_get_immune_markers_specific(self):
        """Test get_immune_markers for specific cell type."""
        from splisosm.gene_lists import get_immune_markers

        t_cell_markers = get_immune_markers('t_cells')
        assert 'CD3D' in t_cell_markers
        assert 'CD68' not in t_cell_markers

    def test_get_immune_markers_invalid(self):
        """Test get_immune_markers raises on invalid cell type."""
        from splisosm.gene_lists import get_immune_markers

        with pytest.raises(ValueError):
            get_immune_markers('invalid_cell_type')


class TestGenesInData:
    """Tests for genes_in_data function."""

    def test_genes_in_data(self):
        """Test filtering genes to those in data."""
        from splisosm.gene_lists import genes_in_data

        desired = ['IL33', 'VEGFA', 'NOTINDATA']
        available = ['IL33', 'VEGFA', 'OTHER']

        present = genes_in_data(desired, available)
        assert present == ['IL33', 'VEGFA']

    def test_genes_in_data_empty(self):
        """Test with no overlap."""
        from splisosm.gene_lists import genes_in_data

        present = genes_in_data(['A', 'B'], ['C', 'D'])
        assert present == []


class TestSummarize:
    """Tests for summarize_gene_lists."""

    def test_summarize_gene_lists(self):
        """Test summarize_gene_lists returns counts."""
        from splisosm.gene_lists import summarize_gene_lists

        summary = summarize_gene_lists()

        assert isinstance(summary, dict)
        assert 'rbp_regulators' in summary
        assert 'hypoxia_markers' in summary
        assert 'cytokine_interleukins' in summary

        # All counts should be positive
        for key, count in summary.items():
            assert count > 0, f"Empty count for {key}"

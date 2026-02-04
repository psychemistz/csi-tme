"""SCALPEL validation for SPLISOSM spatial isoform patterns.

SCALPEL (Single Cell ALternative PolyAdenylation ELucidation) enables
isoform-level quantification from standard 3' scRNA-seq data.

This module provides utilities for validating SPLISOSM spatial findings
with SCALPEL single-cell data.

References:
    Ake et al. (2025) SCALPEL: Single cell alternative polyadenylation
    elucidation. Nat Commun 16:6402.
    GitHub: https://github.com/plasslab/SCALPEL
"""

from dataclasses import dataclass
from typing import Optional, List, Dict, Tuple, Any

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

try:
    from scipy import stats as scipy_stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def _check_dependencies():
    """Check required dependencies."""
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas required. Install with: pip install pandas")
    if not ANNDATA_AVAILABLE:
        raise ImportError("anndata required. Install with: pip install anndata")
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy required. Install with: pip install scipy")


@dataclass
class ValidationResult:
    """Result of SCALPEL validation for a gene.

    Attributes
    ----------
    gene : str
        Gene name
    splisosm_pvalue : float
        SPLISOSM spatial variability p-value
    concordance_statistic : float
        Mann-Whitney U statistic for cell type differences
    concordance_pvalue : float
        P-value for concordance test
    cell_type_1 : str
        First cell type compared
    cell_type_2 : str
        Second cell type compared
    mean_ratio_diff : float
        Difference in mean isoform ratios between cell types
    """
    gene: str
    splisosm_pvalue: float
    concordance_statistic: float
    concordance_pvalue: float
    cell_type_1: str
    cell_type_2: str
    mean_ratio_diff: float


def load_scalpel_results(
    apadge_path: str,
    gene_column: str = 'gene',
    isoform_column: str = 'isoform'
) -> "ad.AnnData":
    """Load SCALPEL APADGE (Alternative PolyAdenylation DifferentGene Expression) matrix.

    SCALPEL outputs isoform-level expression quantified from 3' scRNA-seq.
    The APADGE matrix contains per-cell isoform usage estimates.

    Parameters
    ----------
    apadge_path : str
        Path to SCALPEL output file (.h5ad, .csv, or .tsv)
    gene_column : str
        Column name for gene identifiers
    isoform_column : str
        Column name for isoform identifiers

    Returns
    -------
    adata : AnnData
        AnnData object with:
        - X: isoform expression matrix (cells x isoforms)
        - var: isoform metadata (gene, isoform name)
        - obs: cell metadata
    """
    _check_dependencies()

    if apadge_path.endswith('.h5ad'):
        adata = ad.read_h5ad(apadge_path)
    elif apadge_path.endswith('.csv'):
        df = pd.read_csv(apadge_path, index_col=0)
        adata = ad.AnnData(X=df.values, obs=pd.DataFrame(index=df.index))
        adata.var_names = df.columns.tolist()
    elif apadge_path.endswith('.tsv'):
        df = pd.read_csv(apadge_path, sep='\t', index_col=0)
        adata = ad.AnnData(X=df.values, obs=pd.DataFrame(index=df.index))
        adata.var_names = df.columns.tolist()
    else:
        raise ValueError(f"Unsupported file format: {apadge_path}")

    return adata


def extract_gene_isoforms(
    scalpel_adata: "ad.AnnData",
    gene: str,
    isoform_separator: str = '-'
) -> Tuple[NDArray, List[str]]:
    """Extract isoform expression for a specific gene.

    Parameters
    ----------
    scalpel_adata : AnnData
        SCALPEL output
    gene : str
        Gene name
    isoform_separator : str
        Separator between gene and isoform ID

    Returns
    -------
    isoform_expr : ndarray
        Isoform expression matrix (cells x isoforms)
    isoform_names : list
        Isoform names
    """
    # Find isoforms for this gene
    var_names = scalpel_adata.var_names.tolist()
    isoform_mask = [v.startswith(f"{gene}{isoform_separator}") for v in var_names]
    isoform_indices = [i for i, m in enumerate(isoform_mask) if m]

    if len(isoform_indices) == 0:
        # Try without separator (gene might be the full name)
        isoform_mask = [v == gene for v in var_names]
        isoform_indices = [i for i, m in enumerate(isoform_mask) if m]

    if len(isoform_indices) == 0:
        raise ValueError(f"Gene '{gene}' not found in SCALPEL data")

    isoform_names = [var_names[i] for i in isoform_indices]

    if hasattr(scalpel_adata.X, 'toarray'):
        isoform_expr = scalpel_adata.X[:, isoform_indices].toarray()
    else:
        isoform_expr = np.array(scalpel_adata.X[:, isoform_indices])

    return isoform_expr, isoform_names


def compute_isoform_ratios(
    isoform_expr: NDArray,
    min_total: float = 1.0
) -> Tuple[NDArray, NDArray]:
    """Compute isoform ratios from expression.

    Parameters
    ----------
    isoform_expr : ndarray
        Isoform expression (cells x isoforms)
    min_total : float
        Minimum total expression to compute ratios

    Returns
    -------
    ratios : ndarray
        Isoform ratios (cells x isoforms)
    valid_mask : ndarray
        Boolean mask for cells with valid ratios
    """
    totals = isoform_expr.sum(axis=1)
    valid_mask = totals >= min_total

    ratios = np.zeros_like(isoform_expr, dtype=float)
    ratios[valid_mask] = isoform_expr[valid_mask] / totals[valid_mask, np.newaxis]

    return ratios, valid_mask


def validate_svp_genes(
    splisosm_results: "pd.DataFrame",
    scalpel_adata: "ad.AnnData",
    cell_type_column: str = 'cell_type',
    spatial_region_mapping: Optional[Dict[str, str]] = None,
    min_cells_per_group: int = 10,
    isoform_separator: str = '-'
) -> "pd.DataFrame":
    """Validate SPLISOSM findings with SCALPEL single-cell data.

    For genes with significant spatial isoform variation in SPLISOSM,
    tests whether cell types enriched in different spatial regions
    show differential isoform usage in single-cell data.

    Parameters
    ----------
    splisosm_results : DataFrame
        Results from SplisosmNP.test_spatial_variability()
        Must have columns: gene, pvalue_ir (or padj_ir)
    scalpel_adata : AnnData
        SCALPEL output with cell type annotations
    cell_type_column : str
        Column in scalpel_adata.obs containing cell type labels
    spatial_region_mapping : dict, optional
        Maps spatial regions to cell types expected there.
        E.g., {'tumor_core': ['Malignant', 'TAM'], 'stroma': ['CAF', 'T_cell']}
        If None, tests all pairwise cell type comparisons.
    min_cells_per_group : int
        Minimum cells per group for valid comparison
    isoform_separator : str
        Separator between gene and isoform ID

    Returns
    -------
    validation : DataFrame
        Validation results with columns:
        - gene: Gene name
        - splisosm_pvalue: Original SPLISOSM p-value
        - cell_type_1, cell_type_2: Cell types compared
        - concordance_statistic: Mann-Whitney U statistic
        - concordance_pvalue: P-value for differential usage
        - mean_ratio_diff: Mean isoform ratio difference
        - validated: True if concordance_pvalue < 0.05
    """
    _check_dependencies()

    if cell_type_column not in scalpel_adata.obs.columns:
        raise ValueError(f"Cell type column '{cell_type_column}' not found")

    cell_types = scalpel_adata.obs[cell_type_column].unique().tolist()

    results = []

    # Get significant genes from SPLISOSM
    pval_col = 'padj_ir' if 'padj_ir' in splisosm_results.columns else 'pvalue_ir'
    sig_genes = splisosm_results[splisosm_results[pval_col] < 0.05]['gene'].tolist()

    for gene in sig_genes:
        splisosm_pval = splisosm_results[splisosm_results['gene'] == gene][pval_col].values[0]

        try:
            isoform_expr, isoform_names = extract_gene_isoforms(
                scalpel_adata, gene, isoform_separator
            )
        except ValueError:
            continue

        if len(isoform_names) < 2:
            continue

        ratios, valid_mask = compute_isoform_ratios(isoform_expr)

        # Get cell type comparisons
        if spatial_region_mapping is not None:
            # Use provided mapping
            comparisons = []
            regions = list(spatial_region_mapping.keys())
            for i in range(len(regions)):
                for j in range(i + 1, len(regions)):
                    for ct1 in spatial_region_mapping[regions[i]]:
                        for ct2 in spatial_region_mapping[regions[j]]:
                            if ct1 in cell_types and ct2 in cell_types:
                                comparisons.append((ct1, ct2))
        else:
            # All pairwise comparisons
            comparisons = []
            for i in range(len(cell_types)):
                for j in range(i + 1, len(cell_types)):
                    comparisons.append((cell_types[i], cell_types[j]))

        for ct1, ct2 in comparisons:
            # Get cells for each type
            mask1 = (scalpel_adata.obs[cell_type_column] == ct1).values & valid_mask
            mask2 = (scalpel_adata.obs[cell_type_column] == ct2).values & valid_mask

            n1 = mask1.sum()
            n2 = mask2.sum()

            if n1 < min_cells_per_group or n2 < min_cells_per_group:
                continue

            # Compare primary isoform ratio (first isoform)
            ratios_1 = ratios[mask1, 0]
            ratios_2 = ratios[mask2, 0]

            # Mann-Whitney U test
            stat, pval = scipy_stats.mannwhitneyu(
                ratios_1, ratios_2, alternative='two-sided'
            )

            mean_diff = np.mean(ratios_1) - np.mean(ratios_2)

            results.append({
                'gene': gene,
                'splisosm_pvalue': splisosm_pval,
                'cell_type_1': ct1,
                'cell_type_2': ct2,
                'concordance_statistic': stat,
                'concordance_pvalue': pval,
                'mean_ratio_diff': mean_diff,
                'n_cells_1': n1,
                'n_cells_2': n2,
                'validated': pval < 0.05
            })

    if len(results) == 0:
        return pd.DataFrame(columns=[
            'gene', 'splisosm_pvalue', 'cell_type_1', 'cell_type_2',
            'concordance_statistic', 'concordance_pvalue', 'mean_ratio_diff',
            'n_cells_1', 'n_cells_2', 'validated'
        ])

    return pd.DataFrame(results).sort_values('concordance_pvalue')


def compute_isoform_concordance(
    splisosm_results: "pd.DataFrame",
    scalpel_results: "pd.DataFrame",
    pvalue_threshold: float = 0.05
) -> Dict[str, Any]:
    """Compute overall concordance between SPLISOSM and SCALPEL.

    Parameters
    ----------
    splisosm_results : DataFrame
        SPLISOSM results
    scalpel_results : DataFrame
        Validation results from validate_svp_genes()
    pvalue_threshold : float
        P-value threshold for significance

    Returns
    -------
    concordance : dict
        Concordance metrics:
        - n_tested: Number of genes tested in both
        - n_validated: Number of genes validated
        - validation_rate: Fraction validated
        - fisher_pvalue: Fisher's exact test p-value
    """
    _check_dependencies()

    if len(scalpel_results) == 0:
        return {
            'n_tested': 0,
            'n_validated': 0,
            'validation_rate': 0.0,
            'fisher_pvalue': 1.0
        }

    # Get unique genes
    tested_genes = scalpel_results['gene'].unique()
    validated_genes = scalpel_results[scalpel_results['validated']]['gene'].unique()

    n_tested = len(tested_genes)
    n_validated = len(validated_genes)
    validation_rate = n_validated / n_tested if n_tested > 0 else 0.0

    # Fisher's exact test against random expectation (5%)
    # Contingency table: [[validated, not_validated], [expected_validated, expected_not]]
    expected_validated = int(n_tested * 0.05)
    expected_not = n_tested - expected_validated
    actual_not = n_tested - n_validated

    contingency = [
        [n_validated, actual_not],
        [expected_validated, expected_not]
    ]

    try:
        _, fisher_pvalue = scipy_stats.fisher_exact(contingency, alternative='greater')
    except ValueError:
        fisher_pvalue = 1.0

    return {
        'n_tested': n_tested,
        'n_validated': n_validated,
        'validation_rate': validation_rate,
        'fisher_pvalue': fisher_pvalue,
        'validated_genes': list(validated_genes)
    }

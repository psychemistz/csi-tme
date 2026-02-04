"""Multi-sample meta-analysis for combining SPLISOSM results.

Implements Fisher's method and Stouffer's Z-score for combining
p-values across 1000+ samples.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy import stats


@dataclass
class MetaAnalysisResult:
    """Result of meta-analysis for a single gene.

    Attributes
    ----------
    gene : str
        Gene name
    combined_pvalue : float
        Combined p-value across samples
    n_samples : int
        Number of samples with valid results
    method : str
        Method used ('fisher' or 'stouffer')
    statistic : float
        Test statistic (chi-square for Fisher, Z for Stouffer)
    individual_pvalues : ndarray
        P-values from individual samples
    """
    gene: str
    combined_pvalue: float
    n_samples: int
    method: str
    statistic: float
    individual_pvalues: Optional[NDArray] = None


def fisher_combine(pvalues: NDArray, min_pvalue: float = 1e-300) -> tuple:
    """Combine p-values using Fisher's method.

    χ² = -2 Σ log(pᵢ) ~ χ²_{2k}

    Parameters
    ----------
    pvalues : ndarray
        Array of p-values
    min_pvalue : float
        Floor for p-values (for numerical stability)

    Returns
    -------
    combined_pvalue : float
    statistic : float
        Chi-square statistic
    """
    # Clamp p-values
    pvalues = np.clip(pvalues, min_pvalue, 1.0)

    # Chi-square statistic
    chi2 = -2 * np.sum(np.log(pvalues))

    # Degrees of freedom
    df = 2 * len(pvalues)

    # Combined p-value
    combined_pvalue = stats.chi2.sf(chi2, df)

    return combined_pvalue, chi2


def stouffer_combine(
    pvalues: NDArray,
    weights: Optional[NDArray] = None,
    min_pvalue: float = 1e-300
) -> tuple:
    """Combine p-values using Stouffer's Z-score method.

    Z = Σᵢ wᵢ Φ⁻¹(1-pᵢ) / √(Σᵢ wᵢ²)

    More robust to heterogeneity than Fisher's method.

    Parameters
    ----------
    pvalues : ndarray
        Array of p-values
    weights : ndarray, optional
        Weights for each p-value (e.g., sample size)
    min_pvalue : float
        Floor for p-values

    Returns
    -------
    combined_pvalue : float
    statistic : float
        Z-score
    """
    # Clamp p-values
    pvalues = np.clip(pvalues, min_pvalue, 1 - min_pvalue)

    # Convert to Z-scores
    z_scores = stats.norm.isf(pvalues)

    if weights is None:
        weights = np.ones(len(pvalues))

    weights = np.array(weights)

    # Weighted sum
    z_combined = np.sum(weights * z_scores) / np.sqrt(np.sum(weights ** 2))

    # Combined p-value
    combined_pvalue = stats.norm.sf(z_combined)

    return combined_pvalue, z_combined


def hierarchical_fdr(
    pvalues_matrix: NDArray,
    within_sample_fdr: float = 0.1,
    meta_fdr: float = 0.05
) -> NDArray:
    """Apply hierarchical FDR correction.

    1. Within-sample BH correction
    2. Meta-analysis on filtered p-values

    Parameters
    ----------
    pvalues_matrix : ndarray
        Matrix of p-values (n_genes, n_samples)
    within_sample_fdr : float
        FDR threshold for within-sample correction
    meta_fdr : float
        FDR threshold for meta-analysis

    Returns
    -------
    significant_mask : ndarray
        Boolean mask of significant genes
    """
    n_genes, n_samples = pvalues_matrix.shape

    # Within-sample BH correction
    corrected_pvalues = np.zeros_like(pvalues_matrix)

    for j in range(n_samples):
        pvals = pvalues_matrix[:, j]
        valid_mask = ~np.isnan(pvals)

        if valid_mask.sum() == 0:
            continue

        # BH correction
        sorted_idx = np.argsort(pvals[valid_mask])
        sorted_pvals = pvals[valid_mask][sorted_idx]
        n_valid = len(sorted_pvals)

        thresholds = within_sample_fdr * np.arange(1, n_valid + 1) / n_valid
        adjusted = np.minimum.accumulate(sorted_pvals[::-1] * n_valid /
                                         np.arange(n_valid, 0, -1))[::-1]
        adjusted = np.minimum(adjusted, 1.0)

        # Map back to original indices
        corrected = np.ones(n_genes)
        corrected[valid_mask] = np.ones(n_valid)
        valid_indices = np.where(valid_mask)[0]
        for i, idx in enumerate(sorted_idx):
            corrected[valid_indices[idx]] = adjusted[i]

        corrected_pvalues[:, j] = corrected

    # Meta-analysis on genes passing within-sample filter
    meta_pvalues = np.ones(n_genes)

    for i in range(n_genes):
        gene_pvals = corrected_pvalues[i, :]
        valid_pvals = gene_pvals[~np.isnan(gene_pvals) & (gene_pvals < 1)]

        if len(valid_pvals) >= 2:
            meta_pvalues[i], _ = stouffer_combine(valid_pvals)

    # Final BH correction on meta p-values
    sorted_idx = np.argsort(meta_pvalues)
    sorted_meta = meta_pvalues[sorted_idx]
    thresholds = meta_fdr * np.arange(1, n_genes + 1) / n_genes

    significant_mask = np.zeros(n_genes, dtype=bool)
    if np.any(sorted_meta <= thresholds):
        k = np.max(np.where(sorted_meta <= thresholds)[0]) + 1
        significant_mask[sorted_idx[:k]] = True

    return significant_mask


def meta_analyze_results(
    sample_results: list,
    test: str = 'ir',
    method: str = 'stouffer'
) -> dict:
    """Run meta-analysis across multiple sample results.

    Parameters
    ----------
    sample_results : list
        List of SampleResult objects
    test : str
        Which test to use: 'ir', 'gc', or 'ic'
    method : str
        Combination method: 'fisher' or 'stouffer'

    Returns
    -------
    results : dict
        Maps gene -> MetaAnalysisResult
    """
    pvalue_key = f'pvalue_{test}'

    # Collect p-values for each gene across samples
    gene_pvalues = {}

    for sample_result in sample_results:
        for gene_result in sample_result.gene_results:
            gene = gene_result.gene
            pval = getattr(gene_result, pvalue_key)

            if gene not in gene_pvalues:
                gene_pvalues[gene] = []
            gene_pvalues[gene].append(pval)

    # Combine p-values
    combine_func = fisher_combine if method == 'fisher' else stouffer_combine
    results = {}

    for gene, pvalues in gene_pvalues.items():
        pvalues = np.array(pvalues)
        valid_pvalues = pvalues[~np.isnan(pvalues)]

        if len(valid_pvalues) < 2:
            continue

        combined_pvalue, statistic = combine_func(valid_pvalues)

        results[gene] = MetaAnalysisResult(
            gene=gene,
            combined_pvalue=combined_pvalue,
            n_samples=len(valid_pvalues),
            method=method,
            statistic=statistic,
            individual_pvalues=valid_pvalues
        )

    return results


def results_to_dataframe(meta_results: dict):
    """Convert meta-analysis results to DataFrame.

    Parameters
    ----------
    meta_results : dict
        Output of meta_analyze_results

    Returns
    -------
    df : DataFrame
    """
    import pandas as pd

    records = []
    for gene, result in meta_results.items():
        records.append({
            'gene': result.gene,
            'combined_pvalue': result.combined_pvalue,
            'n_samples': result.n_samples,
            'statistic': result.statistic,
            'method': result.method,
        })

    df = pd.DataFrame(records)
    df = df.sort_values('combined_pvalue')

    return df

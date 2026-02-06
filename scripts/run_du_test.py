#!/usr/bin/env python3
"""Run conditional HSIC (DU test) to identify RBP regulators of spatial isoform patterns.

Tests whether RBP expression correlates with isoform usage after removing
spatial confounding via residualization. This is the SPLISOSM "DU test"
(Differential Usage conditional on spatial location).

Target genes: TIER1 validated + top TIER2 spatial-only genes
RBP covariates: 19 RBPs from splisosm/gene_lists.py
Output: reports/du_test_rbp_results.csv, reports/du_test_rbp_summary.csv
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

from splisosm.api import SplisosmNP, benjamini_hochberg
from splisosm.gene_lists import RBP_REGULATORS
from splisosm.statistics.conditional import conditional_hsic
from splisosm.preprocessing.imputation import impute_with_counts

# Data directories
SR_DATA_DIR = project_root / "data/zenodo_16905935/human_glioma_sr"
ONT_DATA_DIR = project_root / "data/zenodo_16905935/human_glioma_ont"
REPORTS_DIR = project_root / "reports"

# TIER1 (validated) + top TIER2 (strongest spatial-only patterns)
TARGET_GENES = [
    # TIER1
    'SPP1', 'CLU', 'HMGB1', 'CD74', 'APP', 'SPARCL1', 'VCAN',
    # Top TIER2 (by best IR p-value)
    'SPARC', 'B2M', 'PTGDS', 'SERPINE2', 'NCAM1', 'C3',
]


def build_isoform_mapping_from_var(adata, gene_col=None):
    """Build gene -> isoform index mapping from adata.var."""
    if gene_col is None:
        for col in ['gene_symbol', 'gene', 'gene_name']:
            if col in adata.var.columns:
                gene_col = col
                break

    if gene_col is None:
        raise ValueError("No gene column found in adata.var")

    mapping = {}
    for i, (iso_name, row) in enumerate(adata.var.iterrows()):
        gene = row[gene_col]
        if gene not in mapping:
            mapping[gene] = []
        mapping[gene].append(i)

    # Keep only genes with 2+ isoforms
    return {g: idx for g, idx in mapping.items() if len(idx) >= 2}


def get_expression_vector(adata, gene_name):
    """Extract expression vector for a single gene (summing isoforms if needed)."""
    # Check var index directly
    if gene_name in adata.var_names:
        idx = list(adata.var_names).index(gene_name)
        x = adata.X[:, idx]
        if hasattr(x, 'toarray'):
            return x.toarray().flatten()
        return np.asarray(x).flatten()

    # Check gene_symbol column for summing isoforms
    for col in ['gene_symbol', 'gene', 'gene_name']:
        if col in adata.var.columns:
            mask = adata.var[col] == gene_name
            if mask.any():
                x = adata.X[:, mask.values]
                if hasattr(x, 'toarray'):
                    x = x.toarray()
                return np.asarray(x).sum(axis=1).flatten()

    return None


def run_du_test_sample(
    adata,
    sample_id,
    target_genes,
    rbp_list,
    min_spots=30,
    rbf_smoothing=1.0,
):
    """Run DU test for all target gene × RBP pairs in one sample.

    Returns DataFrame with columns:
        sample, gene, rbp, statistic, pvalue, n_spots, n_isoforms, rbp_expressed_spots
    """
    # Build isoform mapping
    isoform_mapping = build_isoform_mapping_from_var(adata)

    # Get spatial coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    elif 'X_spatial' in adata.obsm:
        coords = adata.obsm['X_spatial']
    elif 'array_row' in adata.obs and 'array_col' in adata.obs:
        coords = adata.obs[['array_row', 'array_col']].values.astype(float)
    else:
        print(f"  No spatial coordinates in {sample_id}")
        return pd.DataFrame()

    # Get count matrix as dense array
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    # Pre-extract RBP expression vectors
    rbp_vectors = {}
    for rbp in rbp_list:
        vec = get_expression_vector(adata, rbp)
        if vec is not None and (vec > 0).sum() >= min_spots:
            rbp_vectors[rbp] = vec

    if not rbp_vectors:
        print(f"  {sample_id}: No RBPs with sufficient expression")
        return pd.DataFrame()

    results = []

    for gene in target_genes:
        if gene not in isoform_mapping:
            continue

        isoform_indices = isoform_mapping[gene]
        isoform_counts = X[:, isoform_indices]
        gene_totals = isoform_counts.sum(axis=1)

        n_expressed = (gene_totals > 0).sum()
        if n_expressed < min_spots:
            continue

        # Compute imputed isoform proportions
        proportions = impute_with_counts(isoform_counts, gene_totals)

        for rbp, rbp_expr in rbp_vectors.items():
            try:
                result = conditional_hsic(
                    proportions, rbp_expr, coords,
                    method='rbf', rbf_smoothing=rbf_smoothing
                )
                results.append({
                    'sample': sample_id,
                    'gene': gene,
                    'rbp': rbp,
                    'statistic': result.statistic,
                    'pvalue': result.pvalue,
                    'n_spots': adata.n_obs,
                    'n_isoforms': len(isoform_indices),
                    'gene_expressed_spots': int(n_expressed),
                    'rbp_expressed_spots': int((rbp_expr > 0).sum()),
                })
            except Exception as e:
                print(f"    ERROR {gene}×{rbp} in {sample_id}: {e}")
                continue

    return pd.DataFrame(results)


def main():
    print("=" * 70)
    print("Conditional HSIC (DU Test): RBP Regulator Identification")
    print("=" * 70)

    print(f"\nTarget genes ({len(TARGET_GENES)}): {', '.join(TARGET_GENES)}")
    print(f"RBP candidates ({len(RBP_REGULATORS)}): {', '.join(RBP_REGULATORS)}")
    print(f"Max tests per sample: {len(TARGET_GENES) * len(RBP_REGULATORS)}")

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    # Discover samples from both platforms
    samples = []
    if SR_DATA_DIR.exists():
        sr_samples = sorted([
            d for d in SR_DATA_DIR.iterdir()
            if d.is_dir() and not d.name.startswith('.')
        ])
        samples.extend([(s, 'SR') for s in sr_samples])
        print(f"\nSR samples: {len(sr_samples)}")

    if ONT_DATA_DIR.exists():
        ont_samples = sorted([
            d for d in ONT_DATA_DIR.iterdir()
            if d.is_dir() and not d.name.startswith('.')
        ])
        samples.extend([(s, 'ONT') for s in ont_samples])
        print(f"ONT samples: {len(ont_samples)}")

    if not samples:
        print("ERROR: No sample directories found")
        return

    print(f"Total samples: {len(samples)}")

    # Process each sample
    all_results = []

    print("\n" + "-" * 70)
    print("Running DU tests...")
    print("-" * 70)

    for sample_dir, platform in tqdm(samples, desc="Samples"):
        iso_file = sample_dir / "iso.quant.h5ad"
        if not iso_file.exists():
            continue

        sample_id = f"{platform}_{sample_dir.name}"

        try:
            adata = sc.read_h5ad(iso_file)
            df = run_du_test_sample(
                adata, sample_id, TARGET_GENES, RBP_REGULATORS,
                min_spots=30, rbf_smoothing=1.0,
            )
            if len(df) > 0:
                df['platform'] = platform
                all_results.append(df)
                n_sig = (df['pvalue'] < 0.05).sum()
                print(f"  {sample_id}: {len(df)} tests, {n_sig} nominal sig (p<0.05)")
        except Exception as e:
            print(f"  ERROR in {sample_id}: {e}")
            continue

    if not all_results:
        print("\nNo results generated!")
        return

    # Combine and correct
    combined = pd.concat(all_results, ignore_index=True)
    combined['padj'] = benjamini_hochberg(combined['pvalue'].values)

    # Save full results
    output_file = REPORTS_DIR / "du_test_rbp_results.csv"
    combined.to_csv(output_file, index=False)
    print(f"\nSaved {len(combined)} results to: {output_file}")

    # Generate summary: aggregate across samples per gene-RBP pair
    summary = combined.groupby(['gene', 'rbp']).agg(
        n_samples=('sample', 'count'),
        n_sig_nominal=('pvalue', lambda x: (x < 0.05).sum()),
        min_pvalue=('pvalue', 'min'),
        median_pvalue=('pvalue', 'median'),
        mean_statistic=('statistic', 'mean'),
        median_gene_spots=('gene_expressed_spots', 'median'),
        median_rbp_spots=('rbp_expressed_spots', 'median'),
    ).reset_index()

    # Fisher's method for combining p-values across samples
    from scipy.stats import combine_pvalues
    fisher_results = []
    for (gene, rbp), group in combined.groupby(['gene', 'rbp']):
        pvals = group['pvalue'].values
        # Clamp p-values away from 0 to avoid -inf in log
        pvals = np.clip(pvals, 1e-300, 1.0)
        if len(pvals) >= 2:
            fisher_stat, fisher_p = combine_pvalues(pvals, method='fisher')
        else:
            fisher_stat, fisher_p = np.nan, pvals[0]
        fisher_results.append({
            'gene': gene, 'rbp': rbp,
            'fisher_statistic': fisher_stat, 'fisher_pvalue': fisher_p
        })

    fisher_df = pd.DataFrame(fisher_results)
    fisher_df['fisher_padj'] = benjamini_hochberg(fisher_df['fisher_pvalue'].values)

    summary = summary.merge(fisher_df, on=['gene', 'rbp'])
    summary = summary.sort_values('fisher_pvalue')

    summary_file = REPORTS_DIR / "du_test_rbp_summary.csv"
    summary.to_csv(summary_file, index=False)
    print(f"Saved summary to: {summary_file}")

    # Print results
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"\nTotal tests: {len(combined)}")
    print(f"Gene-RBP pairs tested: {len(summary)}")
    print(f"Nominal significant (p<0.05): {(combined['pvalue'] < 0.05).sum()}")
    print(f"FDR significant (padj<0.05): {(combined['padj'] < 0.05).sum()}")

    # Top hits by Fisher's combined p-value
    print("\nTop 20 gene-RBP associations (Fisher's combined p-value):")
    print("-" * 80)
    top = summary.head(20)
    for _, row in top.iterrows():
        sig_marker = "***" if row['fisher_padj'] < 0.01 else "**" if row['fisher_padj'] < 0.05 else "*" if row['fisher_padj'] < 0.1 else ""
        print(f"  {row['gene']:10s} × {row['rbp']:8s}  "
              f"Fisher p={row['fisher_pvalue']:.2e}  "
              f"padj={row['fisher_padj']:.2e}  "
              f"sig in {row['n_sig_nominal']}/{row['n_samples']} samples  "
              f"{sig_marker}")

    # Per-gene summary
    print("\nSignificant RBPs per gene (Fisher padj < 0.05):")
    print("-" * 60)
    sig_summary = summary[summary['fisher_padj'] < 0.05]
    if len(sig_summary) > 0:
        for gene in TARGET_GENES:
            gene_hits = sig_summary[sig_summary['gene'] == gene]
            if len(gene_hits) > 0:
                rbps = ', '.join(gene_hits['rbp'].tolist())
                print(f"  {gene}: {rbps}")
    else:
        print("  No significant associations at FDR < 0.05")
        # Show best hits per gene at relaxed threshold
        print("\nBest RBP per gene (uncorrected):")
        for gene in TARGET_GENES:
            gene_rows = summary[summary['gene'] == gene]
            if len(gene_rows) > 0:
                best = gene_rows.iloc[0]
                print(f"  {gene}: {best['rbp']} (Fisher p={best['fisher_pvalue']:.2e})")

    print("\n" + "=" * 70)
    print("DU test complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()

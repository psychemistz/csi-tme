#!/usr/bin/env python3
"""Run full SPLISOSM analysis (HSIC-GC, HSIC-IR, HSIC-IC) on Visium-SR glioma samples.

This script analyzes isoform data from Visium short-read sequencing.
Note: SR 3' capture detects APA (Alternative PolyAdenylation) patterns,
not internal splice isoforms.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

from splisosm.api import SplisosmNP, benjamini_hochberg
from splisosm.analysis.single_sample import SplisomAnalyzer, build_isoform_mapping

# SR data directory
SR_DATA_DIR = project_root / "data/zenodo_16905935/human_glioma_sr"
REPORTS_DIR = project_root / "reports"

# SecAct high-confidence targets
SECACT_TARGETS = [
    # Pilot validated
    'HMGB1', 'SPP1', 'VCAN', 'C3', 'IGFBP5', 'APP',
    # High expression targets
    'CLU', 'SPARCL1', 'ITM2B', 'AGT', 'EDIL3', 'DKK3',
    'NCAM1', 'CTSB', 'SPARC', 'CCL2',
    # Additional cytokines/secreted proteins
    'FN1', 'COL1A1', 'COL1A2', 'COL3A1', 'TIMP1', 'TIMP2',
    'SERPINE1', 'SERPINE2', 'LGALS1', 'LGALS3', 'MIF',
    'CXCL12', 'CD74', 'B2M', 'PTGDS', 'APOE', 'CST3',
]


def build_isoform_mapping_from_var(adata, gene_col=None):
    """Build isoform mapping from gene symbol column in var."""
    if gene_col is None:
        for col in ['gene_symbol', 'gene', 'gene_name']:
            if col in adata.var.columns:
                gene_col = col
                break

    if gene_col is None:
        raise ValueError("No gene column found in adata.var")

    # Build mapping: gene -> list of isoform indices
    mapping = {}
    for i, (iso_name, row) in enumerate(adata.var.iterrows()):
        gene = row[gene_col]
        if gene not in mapping:
            mapping[gene] = []
        mapping[gene].append(i)

    # Keep only genes with 2+ isoforms
    mapping = {g: indices for g, indices in mapping.items() if len(indices) >= 2}

    return mapping


def analyze_sr_sample(
    sample_dir: Path,
    target_genes: list = None,
    k_neighbors: int = 15,
    rho: float = 0.99,
    min_spots: int = 30
) -> pd.DataFrame:
    """Run full SPLISOSM analysis on a single SR sample.

    Returns results for all three tests: HSIC-GC, HSIC-IR, HSIC-IC.
    """
    sample_id = sample_dir.name
    iso_file = sample_dir / "iso.quant.h5ad"

    if not iso_file.exists():
        print(f"  No iso.quant.h5ad in {sample_id}")
        return pd.DataFrame()

    # Load isoform data
    adata = sc.read_h5ad(iso_file)

    # Ensure spatial coords exist
    if 'spatial' not in adata.obsm:
        if 'X_spatial' in adata.obsm:
            adata.obsm['spatial'] = adata.obsm['X_spatial']
        elif 'array_row' in adata.obs and 'array_col' in adata.obs:
            coords = adata.obs[['array_row', 'array_col']].values.astype(float)
            adata.obsm['spatial'] = coords
        else:
            print(f"  No spatial coordinates in {sample_id}")
            return pd.DataFrame()

    # Build isoform mapping from gene_symbol column
    isoform_mapping = build_isoform_mapping_from_var(adata)

    # Filter to target genes if specified
    if target_genes:
        isoform_mapping = {g: idx for g, idx in isoform_mapping.items() if g in target_genes}

    if not isoform_mapping:
        print(f"  No testable genes in {sample_id}")
        return pd.DataFrame()

    print(f"  {sample_id}: {adata.n_obs} spots, {len(isoform_mapping)} genes with isoforms")

    # Create SplisosmNP model
    model = SplisosmNP(
        adata,
        spatial_key='spatial',
        k_neighbors=k_neighbors,
        rho=rho,
        min_spots_expressed=min_spots,
        use_gpu=False
    )

    # Override the isoform mapping (since var names might not use standard separator)
    model._isoform_mapping = isoform_mapping

    # Run analysis - this returns all three HSIC statistics
    results = model.test_spatial_variability(
        test_type='all',
        genes=list(isoform_mapping.keys()),
        show_progress=False
    )

    if len(results) > 0:
        results['sample'] = sample_id
        results['n_spots'] = adata.n_obs

    return results


def main():
    print("=" * 70)
    print("Full SPLISOSM Analysis (HSIC-GC, HSIC-IR, HSIC-IC) on Visium-SR")
    print("=" * 70)

    # Discover samples
    samples = sorted([d for d in SR_DATA_DIR.iterdir() if d.is_dir() and not d.name.startswith('.')])
    print(f"\nFound {len(samples)} SR samples:")
    for s in samples:
        print(f"  - {s.name}")

    # Process each sample
    all_results = []

    print("\n" + "-" * 70)
    print("Running SPLISOSM analysis...")
    print("-" * 70)

    for sample_dir in tqdm(samples, desc="Samples"):
        try:
            results = analyze_sr_sample(
                sample_dir,
                target_genes=SECACT_TARGETS,
                k_neighbors=15,
                rho=0.99,
                min_spots=30
            )

            if len(results) > 0:
                all_results.append(results)
                n_sig_ir = (results['pvalue_ir'] < 0.05).sum()
                n_sig_gc = (results['pvalue_gc'] < 0.05).sum()
                n_sig_ic = (results['pvalue_ic'] < 0.05).sum()
                print(f"    {sample_dir.name}: Tested {len(results)} genes | "
                      f"Sig: IR={n_sig_ir}, GC={n_sig_gc}, IC={n_sig_ic}")

        except Exception as e:
            print(f"    ERROR in {sample_dir.name}: {e}")
            continue

    if not all_results:
        print("\nNo results generated!")
        return pd.DataFrame()

    # Combine results
    combined = pd.concat(all_results, ignore_index=True)

    # Apply global FDR correction
    combined['padj_ir_global'] = benjamini_hochberg(combined['pvalue_ir'].values)
    combined['padj_gc_global'] = benjamini_hochberg(combined['pvalue_gc'].values)
    combined['padj_ic_global'] = benjamini_hochberg(combined['pvalue_ic'].values)

    # Save results
    output_file = REPORTS_DIR / "sr_splisosm_full_results.csv"
    combined.to_csv(output_file, index=False)
    print(f"\nSaved full results to: {output_file}")

    # Generate summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"Total tests: {len(combined)}")
    print(f"Samples: {combined['sample'].nunique()}")
    print(f"Genes tested: {combined['gene'].nunique()}")

    # Significance summary per test type
    print("\nSignificance by test type (nominal p < 0.05):")
    print("-" * 50)
    for test, col in [('HSIC-IR (isoform ratio)', 'pvalue_ir'),
                      ('HSIC-GC (gene expression)', 'pvalue_gc'),
                      ('HSIC-IC (combined)', 'pvalue_ic')]:
        n_sig = (combined[col] < 0.05).sum()
        n_genes = combined[combined[col] < 0.05]['gene'].nunique()
        print(f"  {test}: {n_sig} tests ({n_genes} genes)")

    # Top genes by each test
    print("\nTop 10 genes by HSIC-GC (spatial gene expression variation):")
    print("-" * 60)
    best_gc = combined.loc[combined.groupby('gene')['pvalue_gc'].idxmin()]
    best_gc = best_gc.sort_values('pvalue_gc').head(10)
    print(best_gc[['gene', 'sample', 'pvalue_gc', 'padj_gc']].to_string(index=False))

    print("\nTop 10 genes by HSIC-IR (spatial isoform ratio variation):")
    print("-" * 60)
    best_ir = combined.loc[combined.groupby('gene')['pvalue_ir'].idxmin()]
    best_ir = best_ir.sort_values('pvalue_ir').head(10)
    print(best_ir[['gene', 'sample', 'pvalue_ir', 'padj_ir']].to_string(index=False))

    print("\nTop 10 genes by HSIC-IC (combined variation):")
    print("-" * 60)
    best_ic = combined.loc[combined.groupby('gene')['pvalue_ic'].idxmin()]
    best_ic = best_ic.sort_values('pvalue_ic').head(10)
    print(best_ic[['gene', 'sample', 'pvalue_ic', 'padj_ic']].to_string(index=False))

    # Create summary per gene
    summary = combined.groupby('gene').agg({
        'pvalue_ir': ['min', 'mean'],
        'pvalue_gc': ['min', 'mean'],
        'pvalue_ic': ['min', 'mean'],
        'sample': 'count'
    }).reset_index()
    summary.columns = ['gene', 'min_pval_ir', 'mean_pval_ir',
                      'min_pval_gc', 'mean_pval_gc',
                      'min_pval_ic', 'mean_pval_ic', 'n_samples']

    # Count significant samples per gene
    for test, col in [('ir', 'pvalue_ir'), ('gc', 'pvalue_gc'), ('ic', 'pvalue_ic')]:
        sig_counts = combined[combined[col] < 0.05].groupby('gene').size()
        summary[f'n_sig_{test}'] = summary['gene'].map(sig_counts).fillna(0).astype(int)

    summary_file = REPORTS_DIR / "sr_splisosm_gene_summary.csv"
    summary.to_csv(summary_file, index=False)
    print(f"\nSaved gene summary to: {summary_file}")

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)

    return combined


if __name__ == "__main__":
    results = main()

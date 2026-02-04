#!/usr/bin/env python3
"""Run SPLISOSM HSIC-IR on ONT glioma samples for SecAct target genes.

This script analyzes full isoform data from Visium-ONT long-read sequencing,
which can detect true splice isoforms (unlike 3' SR which detects APA only).
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
from splisosm.analysis.single_sample import build_isoform_mapping


# ONT data directory
ONT_DATA_DIR = project_root / "data/zenodo_16905935/human_glioma_ont"

# SecAct high-confidence targets (validated in pilot)
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


def load_ont_sample(sample_dir: Path) -> sc.AnnData:
    """Load ONT isoform data for a sample."""
    iso_file = sample_dir / "iso.quant.h5ad"
    if not iso_file.exists():
        raise FileNotFoundError(f"No iso.quant.h5ad in {sample_dir}")

    adata = sc.read_h5ad(iso_file)
    return adata


def build_isoform_mapping_from_var(adata: sc.AnnData) -> dict:
    """Build isoform mapping using gene_symbol column in var.

    This handles both known isoforms (GENE-201) and novel isoforms (GENE_Iso_1).
    """
    gene_col = None
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


def analyze_ont_sample(
    adata: sc.AnnData,
    sample_id: str,
    target_genes: list = None,
    k_neighbors: int = 15,
    rho: float = 0.99,
    min_spots: int = 30
) -> pd.DataFrame:
    """Run SPLISOSM on a single ONT sample.

    Parameters
    ----------
    adata : AnnData
        Isoform-level counts with spatial coordinates
    sample_id : str
        Sample identifier
    target_genes : list, optional
        Genes to analyze (None = all genes with isoforms)
    k_neighbors : int
        ICAR kernel neighbors
    rho : float
        ICAR autocorrelation
    min_spots : int
        Minimum spots for testing

    Returns
    -------
    results : DataFrame
        HSIC-IR results for tested genes
    """
    from splisosm.analysis.single_sample import SplisomAnalyzer
    from splisosm.io.anndata import to_splisosm_format

    # Ensure spatial coords exist
    if 'spatial' not in adata.obsm:
        if 'X_spatial' in adata.obsm:
            adata.obsm['spatial'] = adata.obsm['X_spatial']
        else:
            if 'array_row' in adata.obs and 'array_col' in adata.obs:
                coords = adata.obs[['array_row', 'array_col']].values.astype(float)
                adata.obsm['spatial'] = coords
            else:
                raise ValueError(f"No spatial coordinates found for {sample_id}")

    # Build isoform mapping from gene_symbol column
    isoform_mapping = build_isoform_mapping_from_var(adata)

    # Filter to target genes if specified
    if target_genes:
        isoform_mapping = {g: idx for g, idx in isoform_mapping.items() if g in target_genes}

    if not isoform_mapping:
        return pd.DataFrame()

    # Extract counts and coordinates
    counts, coords, genes = to_splisosm_format(adata, spatial_key='spatial')

    # Create analyzer
    analyzer = SplisomAnalyzer(
        k_neighbors=k_neighbors,
        rho=rho,
        min_spots_expressed=min_spots,
        use_gpu=False
    )

    # Run analysis
    result = analyzer.analyze_sample(
        counts=counts,
        coords=coords,
        genes=genes,
        isoform_mapping=isoform_mapping,
        sample_id=sample_id,
        show_progress=False
    )

    # Convert to DataFrame
    df = result.to_dataframe()

    if len(df) > 0:
        df['sample_id'] = sample_id
        df['sample_type'] = 'GBM' if 'GBM' in sample_id else 'DMG'
        # Add BH-corrected p-values
        df['padj_ir'] = benjamini_hochberg(df['pvalue_ir'].values)
        df['padj_gc'] = benjamini_hochberg(df['pvalue_gc'].values)
        df['padj_ic'] = benjamini_hochberg(df['pvalue_ic'].values)

    return df


def aggregate_existing_results(ont_dir: Path) -> pd.DataFrame:
    """Aggregate pre-computed sv_results.csv files."""
    results = []

    for sample_dir in sorted(ont_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        sv_file = sample_dir / "sv_results.csv"
        if sv_file.exists():
            df = pd.read_csv(sv_file)
            df['sample_id'] = sample_dir.name
            df['sample_type'] = 'GBM' if 'GBM' in sample_dir.name else 'DMG'
            results.append(df)

    if results:
        return pd.concat(results, ignore_index=True)
    return pd.DataFrame()


def main():
    print("=" * 60)
    print("SPLISOSM Analysis of ONT Glioma Samples")
    print("=" * 60)

    # Discover samples
    samples = sorted([d for d in ONT_DATA_DIR.iterdir() if d.is_dir()])
    print(f"\nFound {len(samples)} ONT samples:")
    for s in samples:
        print(f"  - {s.name}")

    # Part 1: Aggregate existing full-genome results
    print("\n" + "-" * 60)
    print("Part 1: Aggregating existing pre-computed results")
    print("-" * 60)

    existing_results = aggregate_existing_results(ONT_DATA_DIR)
    if len(existing_results) > 0:
        print(f"Loaded {len(existing_results)} gene results from {existing_results['sample_id'].nunique()} samples")

        # Filter to SecAct targets
        secact_existing = existing_results[existing_results['gene'].isin(SECACT_TARGETS)]
        print(f"SecAct targets in pre-computed results: {secact_existing['gene'].nunique()} genes")

        # Show top hits
        if len(secact_existing) > 0:
            secact_existing = secact_existing.sort_values('pvalue_hsic-ir')
            print("\nTop 10 SecAct genes by HSIC-IR p-value:")
            top_cols = ['gene', 'sample_id', 'n_iso', 'pvalue_hsic-ir', 'padj_hsic-ir']
            print(secact_existing[top_cols].head(10).to_string(index=False))

    # Part 2: Run fresh analysis on SecAct targets
    print("\n" + "-" * 60)
    print("Part 2: Running SPLISOSM on SecAct target genes")
    print("-" * 60)

    all_results = []

    for sample_dir in tqdm(samples, desc="Analyzing samples"):
        sample_id = sample_dir.name

        try:
            # Load data
            adata = load_ont_sample(sample_dir)
            n_spots = adata.n_obs
            n_isoforms = adata.n_vars

            # Extract gene names from isoform IDs
            if 'gene_symbol' in adata.var.columns:
                gene_col = 'gene_symbol'
            elif 'gene' in adata.var.columns:
                gene_col = 'gene'
            elif 'gene_name' in adata.var.columns:
                gene_col = 'gene_name'
            else:
                gene_col = None

            n_genes = adata.var[gene_col].nunique() if gene_col else 'unknown'
            print(f"\n{sample_id}: {n_spots} spots, {n_isoforms} isoforms, {n_genes} genes")

            # Run analysis
            results = analyze_ont_sample(
                adata,
                sample_id,
                target_genes=SECACT_TARGETS,
                k_neighbors=15,
                rho=0.99,
                min_spots=30
            )

            if len(results) > 0:
                all_results.append(results)
                n_sig = (results['pvalue_ir'] < 0.05).sum()
                print(f"  Tested {len(results)} genes, {n_sig} significant at p<0.05")
            else:
                print(f"  No testable genes found")

        except Exception as e:
            print(f"  ERROR: {e}")
            continue

    # Combine results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)

        # Apply FDR correction across all tests
        combined['padj_global'] = benjamini_hochberg(combined['pvalue_ir'].values)

        # Sort by p-value
        combined = combined.sort_values('pvalue_ir')

        # Save results
        output_file = ONT_DATA_DIR / "ont_secact_splisosm_results.csv"
        combined.to_csv(output_file, index=False)
        print(f"\nSaved results to: {output_file}")

        # Summary statistics
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"Total tests: {len(combined)}")
        print(f"Samples: {combined['sample_id'].nunique()}")
        print(f"Genes tested: {combined['gene'].nunique()}")

        # Significant at different thresholds
        for thresh in [0.05, 0.01, 0.001, 1e-5]:
            n_sig = (combined['pvalue_ir'] < thresh).sum()
            n_genes = combined[combined['pvalue_ir'] < thresh]['gene'].nunique()
            print(f"Significant at p<{thresh}: {n_sig} tests, {n_genes} unique genes")

        # Top genes by frequency of significance
        sig_genes = combined[combined['pvalue_ir'] < 0.05].groupby('gene').size()
        sig_genes = sig_genes.sort_values(ascending=False)

        print("\nGenes significant in most samples (p<0.05):")
        print(sig_genes.head(15).to_string())

        # Show best results per gene
        print("\nBest p-value per gene:")
        best_per_gene = combined.loc[combined.groupby('gene')['pvalue_ir'].idxmin()]
        best_per_gene = best_per_gene.sort_values('pvalue_ir')
        display_cols = ['gene', 'sample_id', 'n_isoforms', 'n_spots_expressed',
                       'pvalue_ir', 'padj_ir', 'pvalue_gc']
        # Filter to columns that exist
        display_cols = [c for c in display_cols if c in best_per_gene.columns]
        print(best_per_gene[display_cols].head(20).to_string(index=False))

    else:
        print("\nNo results generated!")

    return combined if all_results else pd.DataFrame()


if __name__ == "__main__":
    results = main()

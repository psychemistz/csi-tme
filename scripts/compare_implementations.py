#!/usr/bin/env python3
"""Compare local SPLISOSM implementation with official precomputed results.

This script identifies exact discrepancies between implementations.
"""

import sys
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path
from scipy import stats
from sklearn.neighbors import kneighbors_graph
from scipy.sparse.linalg import splu
from scipy import sparse

# Add local splisosm to path
sys.path.insert(0, '/vf/users/parks34/projects/3csitme')

from splisosm.kernels.icar import ICARKernel
from splisosm.statistics.hsic import center_kernel, compute_hsic_statistic, get_kernel_eigenvalues
from splisosm.statistics.pvalue import liu_sf, hsic_pvalue, compute_composite_eigenvalues


def load_sample_data(sample_dir: Path):
    """Load sample data and reference results."""
    iso_path = sample_dir / 'iso.quant.h5ad'
    ref_path = sample_dir / 'sv_results.csv'

    if not iso_path.exists() or not ref_path.exists():
        return None, None

    adata = ad.read_h5ad(iso_path)
    ref_df = pd.read_csv(ref_path)

    return adata, ref_df


def extract_coords(adata):
    """Extract spatial coordinates."""
    if 'spatial' in adata.obsm:
        return np.array(adata.obsm['spatial'])
    elif 'array_row' in adata.obs and 'array_col' in adata.obs:
        return np.column_stack([adata.obs['array_row'], adata.obs['array_col']])
    else:
        raise ValueError("Cannot find spatial coordinates")


def build_icar_kernel(coords, k_neighbors=6, rho=0.99, regularization=1e-6):
    """Build ICAR spatial kernel matching official implementation."""
    n = coords.shape[0]
    k = min(k_neighbors, n - 1)

    # Build k-NN adjacency
    W = kneighbors_graph(coords, n_neighbors=k, mode='connectivity')
    W = W + W.T  # Symmetrize
    W.data = np.ones_like(W.data)  # Binary

    # Degree matrix and precision
    degrees = np.array(W.sum(axis=1)).flatten()
    D = sparse.diags(degrees, format='csr')
    Q = D - rho * W + sparse.eye(n, format='csr') * regularization

    # Invert via LU factorization
    lu = splu(Q.tocsc())
    L = np.zeros((n, n))
    I = np.eye(n)
    for i in range(n):
        L[:, i] = lu.solve(I[:, i])

    L = (L + L.T) / 2  # Ensure symmetry
    return L


def compute_gene_hsic_ir(ratios, L_sp, n):
    """Compute HSIC-IR exactly as official implementation.

    Official SPLISOSM:
    1. y = ratios (centered per column already in ratios_obs)
    2. K_y = y @ y.T (feature kernel)
    3. hsic_scaled = tr(K_sp @ K_y) where K_sp is centered spatial kernel
    4. hsic_statistic = hsic_scaled / (n-1)^2
    5. p-value: liu_sf(hsic_scaled * n, eigenvalues)
    """
    # Center ratios (mean=0 per column)
    y = ratios - ratios.mean(axis=0, keepdims=True)

    # Feature kernel K_y = y @ y.T
    K_y = y @ y.T

    # Center spatial kernel
    L_sp_centered = center_kernel(L_sp)

    # Center feature kernel (official centers both)
    K_y_centered = center_kernel(K_y)

    # HSIC trace (on centered kernels)
    hsic_trace = np.sum(L_sp_centered * K_y_centered)

    # Normalized statistic
    hsic_statistic = hsic_trace / ((n - 1) ** 2)

    # Eigenvalues for p-value
    eig_sp = np.linalg.eigvalsh(L_sp_centered)
    eig_y = np.linalg.eigvalsh(K_y_centered)

    # Filter small eigenvalues
    eig_sp = eig_sp[np.abs(eig_sp) > 1e-10]
    eig_y = eig_y[np.abs(eig_y) > 1e-10]

    # Composite eigenvalues (outer product)
    lambda_xy = np.outer(eig_sp, eig_y).flatten()

    # P-value via Liu's method
    # Official: liu_sf(hsic_scaled * n, eigenvalues)
    # hsic_scaled = hsic_trace (before dividing by (n-1)^2)
    scaled_stat_for_pvalue = hsic_trace * n / ((n-1)**2)  # = hsic_statistic * n * (n-1)^2 / (n-1)^2

    pvalue, df, delta = liu_sf(hsic_trace * n, lambda_xy)

    return hsic_statistic, pvalue, hsic_trace


def compute_gene_hsic_ir_local(ratios, L_sp, n):
    """Compute HSIC-IR using local implementation."""
    # Center ratios
    y = ratios - ratios.mean(axis=0, keepdims=True)

    # Feature kernel
    K_y = y @ y.T

    # Use local implementation
    hsic_stat = compute_hsic_statistic(K_y, L_sp, centered=False, unbiased=True)

    # Get eigenvalues
    L_sp_centered = center_kernel(L_sp)
    K_y_centered = center_kernel(K_y)
    eig_sp = get_kernel_eigenvalues(L_sp_centered)
    eig_y = get_kernel_eigenvalues(K_y_centered)

    # P-value
    pval_result = hsic_pvalue(hsic_stat, eig_y, eig_sp, n, method='liu')

    return hsic_stat, pval_result.pvalue


def analyze_gene_by_gene(adata, ref_df, k_neighbors=6, rho=0.99, max_genes=50):
    """Analyze gene by gene to find discrepancies."""
    coords = extract_coords(adata)
    n_spots = len(coords)

    print(f"Sample has {n_spots} spots")
    print(f"Building ICAR kernel with k={k_neighbors}, rho={rho}")

    # Build spatial kernel
    L_sp = build_icar_kernel(coords, k_neighbors=k_neighbors, rho=rho)

    # Get ratios from layer
    if 'ratios_obs' in adata.layers:
        ratios_all = np.array(adata.layers['ratios_obs'])
    else:
        raise ValueError("No ratios_obs layer found")

    # Get gene mapping (var names contain gene info)
    var_df = adata.var.copy()
    if 'Gene' in var_df.columns:
        gene_col = 'Gene'
    elif 'gene_symbol' in var_df.columns:
        gene_col = 'gene_symbol'
    else:
        raise ValueError("Cannot find gene column")

    # Build gene -> isoform indices mapping
    gene_to_isoforms = {}
    for idx, gene in enumerate(var_df[gene_col]):
        if gene not in gene_to_isoforms:
            gene_to_isoforms[gene] = []
        gene_to_isoforms[gene].append(idx)

    # Filter to genes with multiple isoforms
    gene_to_isoforms = {g: v for g, v in gene_to_isoforms.items() if len(v) >= 2}

    # Match with reference
    ref_genes = set(ref_df['gene'].values)
    common_genes = [g for g in gene_to_isoforms.keys() if g in ref_genes][:max_genes]

    print(f"Testing {len(common_genes)} genes")

    results = []
    for gene in common_genes:
        iso_indices = gene_to_isoforms[gene]
        ratios = ratios_all[:, iso_indices]

        # Skip if all NaN or constant
        if np.isnan(ratios).all() or ratios.std() < 1e-10:
            continue

        # Replace NaN with column means
        col_means = np.nanmean(ratios, axis=0)
        for j in range(ratios.shape[1]):
            ratios[np.isnan(ratios[:, j]), j] = col_means[j]

        # Compute our way
        stat_local, pval_local = compute_gene_hsic_ir_local(ratios, L_sp, n_spots)

        # Compute matching official (as close as we can)
        stat_official, pval_official, hsic_trace = compute_gene_hsic_ir(ratios, L_sp, n_spots)

        # Get reference
        ref_row = ref_df[ref_df['gene'] == gene]
        if len(ref_row) == 0:
            continue
        ref_pval = ref_row['pvalue_hsic-ir'].values[0]

        results.append({
            'gene': gene,
            'n_isoforms': len(iso_indices),
            'pval_ref': ref_pval,
            'pval_local': pval_local,
            'pval_official_match': pval_official,
            'stat_local': stat_local,
            'stat_official': stat_official,
            'hsic_trace': hsic_trace,
            'log10_diff_local': np.log10(pval_local + 1e-300) - np.log10(ref_pval + 1e-300),
            'log10_diff_official': np.log10(pval_official + 1e-300) - np.log10(ref_pval + 1e-300),
        })

    return pd.DataFrame(results)


def main():
    # Test on MGH258 sample
    sample_dir = Path('/data/parks34/projects/3csitme/data/zenodo_16905935/human_glioma_sr/MGH258')

    print("Loading data...")
    adata, ref_df = load_sample_data(sample_dir)
    if adata is None:
        print("Failed to load data")
        return

    print(f"Loaded {adata.shape[0]} spots, {adata.shape[1]} isoforms")
    print(f"Reference has {len(ref_df)} genes")

    # Analyze
    results_df = analyze_gene_by_gene(adata, ref_df, k_neighbors=6, rho=0.99, max_genes=100)

    print("\n=== Comparison Results ===")
    print(f"Genes tested: {len(results_df)}")

    # Correlation analysis
    log_ref = -np.log10(results_df['pval_ref'].values + 1e-300)
    log_local = -np.log10(results_df['pval_local'].values + 1e-300)
    log_official = -np.log10(results_df['pval_official_match'].values + 1e-300)

    r_local, _ = stats.spearmanr(log_ref, log_local)
    r_official, _ = stats.spearmanr(log_ref, log_official)

    print(f"\nLocal implementation Spearman r: {r_local:.4f}")
    print(f"Official-match Spearman r: {r_official:.4f}")

    print(f"\nMean log10 diff (local): {results_df['log10_diff_local'].mean():.4f}")
    print(f"Mean log10 diff (official): {results_df['log10_diff_official'].mean():.4f}")

    print("\n=== Top 10 genes by discrepancy ===")
    results_df['abs_diff'] = np.abs(results_df['log10_diff_local'])
    top_diff = results_df.nlargest(10, 'abs_diff')
    print(top_diff[['gene', 'n_isoforms', 'pval_ref', 'pval_local', 'log10_diff_local']].to_string())

    # Save results
    output_path = Path('/vf/users/parks34/projects/3csitme/scripts/comparison_results.csv')
    results_df.to_csv(output_path, index=False)
    print(f"\nResults saved to {output_path}")

    return results_df


if __name__ == '__main__':
    results = main()

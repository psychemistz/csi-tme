#!/usr/bin/env python3
"""Test different ratio layers to match official SPLISOSM."""

import numpy as np
import anndata as ad
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.neighbors import kneighbors_graph
from scipy import sparse
from scipy.sparse.linalg import splu
import sys

sys.path.insert(0, '/vf/users/parks34/projects/3csitme')
from splisosm.statistics.hsic import center_kernel


def liu_sf(t, eigenvalues):
    """Liu's approximation."""
    eigenvalues = np.asarray(eigenvalues).flatten()
    eigenvalues = eigenvalues[np.abs(eigenvalues) > 1e-10]

    if len(eigenvalues) == 0:
        return 1.0

    n_terms = len(eigenvalues)
    dofs = np.ones(n_terms)
    deltas = np.zeros(n_terms)

    c1 = np.sum(eigenvalues * (dofs + deltas))
    c2 = np.sum(eigenvalues**2 * (dofs + 2 * deltas))
    c3 = np.sum(eigenvalues**3 * (dofs + 3 * deltas))
    c4 = np.sum(eigenvalues**4 * (dofs + 4 * deltas))

    if c2 <= 0:
        return 1.0

    s1 = c3 / (c2 ** 1.5)
    s2 = c4 / (c2 ** 2)

    if s1 ** 2 > s2:
        a = 1.0 / (s1 - np.sqrt(s1**2 - s2))
        delta_approx = s1 * a**3 - a**2
        df_approx = a**2 - 2 * delta_approx
    else:
        a = 1.0 / s1 if abs(s1) > 1e-10 else 1.0
        delta_approx = 0.0
        df_approx = 1.0 / (s1**2) if abs(s1) > 1e-10 else 1.0

    df_approx = max(df_approx, 0.5)
    delta_approx = max(delta_approx, 0.0)

    mu_x = c1
    sigma_x = np.sqrt(2 * c2)
    mu_chi = df_approx + delta_approx
    sigma_chi = np.sqrt(2 * (df_approx + 2 * delta_approx))

    if sigma_x > 0:
        t_star = (t - mu_x) / sigma_x
        t_final = t_star * sigma_chi + mu_chi
    else:
        t_final = t

    if delta_approx < 1e-10:
        pvalue = stats.chi2.sf(t_final, df_approx)
    else:
        pvalue = stats.ncx2.sf(t_final, df_approx, delta_approx)

    return np.clip(pvalue, 1e-300, 1.0)


def build_spatial_kernel(coords, k=6, rho=0.99, reg=1e-6):
    """Build ICAR spatial kernel."""
    n = len(coords)
    k = min(k, n - 1)

    W = kneighbors_graph(coords, n_neighbors=k, mode='connectivity')
    W = W + W.T
    W.data = np.ones_like(W.data)
    degrees = np.array(W.sum(axis=1)).flatten()
    D = sparse.diags(degrees, format='csr')
    Q = D - rho * W + sparse.eye(n, format='csr') * reg
    lu = splu(Q.tocsc())
    L = np.zeros((n, n))
    I = np.eye(n)
    for i in range(n):
        L[:, i] = lu.solve(I[:, i])
    L = (L + L.T) / 2
    return L


def compute_hsic_pvalue(ratios, K_sp, n):
    """Compute HSIC p-value."""
    # Center ratios
    y = ratios - ratios.mean(axis=0, keepdims=True)

    # HSIC trace
    hsic_trace = np.trace(y.T @ K_sp @ y)

    # Eigenvalues
    eig_sp = np.linalg.eigvalsh(K_sp)
    eig_sp = eig_sp[np.abs(eig_sp) > 1e-5]

    gram_y = y.T @ y
    eig_y = np.linalg.eigvalsh(gram_y)
    eig_y = eig_y[np.abs(eig_y) > 1e-10]

    # Composite eigenvalues
    lambda_xy = np.outer(eig_sp, eig_y).flatten()

    # P-value
    pvalue = liu_sf(hsic_trace * n, lambda_xy)

    return pvalue


def main():
    # Load data
    sample_dir = Path('/data/parks34/projects/3csitme/data/zenodo_16905935/human_glioma_sr/MGH258')
    adata = ad.read_h5ad(sample_dir / 'iso.quant.h5ad')
    ref_df = pd.read_csv(sample_dir / 'sv_results.csv')

    coords_all = np.array(adata.obsm['spatial'])
    n = len(coords_all)

    # Load different layers
    def get_layer(name):
        layer = adata.layers[name]
        if hasattr(layer, 'toarray'):
            return layer.toarray()
        return np.array(layer)

    ratios_obs = get_layer('ratios_obs')
    ratios_smoothed = get_layer('ratios_smoothed')

    # Build spatial kernel and center it
    L_sp = build_spatial_kernel(coords_all, k=6, rho=0.99)
    K_sp = center_kernel(L_sp)

    # Get gene info
    gene_col = 'Gene' if 'Gene' in adata.var.columns else 'gene_symbol'

    # Build gene mapping
    gene_to_isoforms = {}
    for idx, gene in enumerate(adata.var[gene_col]):
        if gene not in gene_to_isoforms:
            gene_to_isoforms[gene] = []
        gene_to_isoforms[gene].append(idx)
    gene_to_isoforms = {g: v for g, v in gene_to_isoforms.items() if len(v) >= 2}

    # Test genes
    test_genes = ['H3-3A', 'NCAM1', 'CALM1', 'SOX2-OT', 'PLP1', 'PTMA', 'GPM6B', 'RAC1', 'SOX2', 'MYL6']

    print(f"Total spots: {n}")
    print("\n" + "=" * 120)
    print("Testing different ratio layers: ratios_obs vs ratios_smoothed")
    print("=" * 120)
    header = f"{'gene':>10} {'p_ref':>15} {'p_obs':>15} {'p_smooth':>15} {'diff_obs':>10} {'diff_smooth':>10}"
    print(header)
    print("-" * 120)

    for gene in test_genes:
        if gene not in gene_to_isoforms:
            continue

        iso_indices = gene_to_isoforms[gene]
        n_iso = len(iso_indices)

        # Get ratios for this gene
        r_obs = ratios_obs[:, iso_indices].copy()
        r_smooth = ratios_smoothed[:, iso_indices].copy()

        # Handle NaN
        for ratios in [r_obs, r_smooth]:
            col_means = np.nanmean(ratios, axis=0)
            for j in range(n_iso):
                ratios[np.isnan(ratios[:, j]), j] = col_means[j]

        # Reference p-value
        ref_row = ref_df[ref_df['gene'] == gene]
        ref_pval = ref_row['pvalue_hsic-ir'].values[0] if len(ref_row) > 0 else np.nan

        # Compute p-values with different layers
        p_obs = compute_hsic_pvalue(r_obs, K_sp, n)
        p_smooth = compute_hsic_pvalue(r_smooth, K_sp, n)

        # Log differences
        diff_obs = np.log10(p_obs + 1e-300) - np.log10(ref_pval + 1e-300)
        diff_smooth = np.log10(p_smooth + 1e-300) - np.log10(ref_pval + 1e-300)

        print(f"{gene:>10} {ref_pval:>15.6e} {p_obs:>15.6e} {p_smooth:>15.6e} {diff_obs:>10.3f} {diff_smooth:>10.3f}")

    # Now test a much larger set of genes and compare correlations
    print("\n" + "=" * 120)
    print("Correlation analysis on all genes")
    print("=" * 120)

    results_obs = []
    results_smooth = []
    ref_pvals = []

    for gene in gene_to_isoforms:
        ref_row = ref_df[ref_df['gene'] == gene]
        if len(ref_row) == 0:
            continue

        iso_indices = gene_to_isoforms[gene]
        n_iso = len(iso_indices)

        r_obs = ratios_obs[:, iso_indices].copy()
        r_smooth = ratios_smoothed[:, iso_indices].copy()

        # Handle NaN
        for ratios in [r_obs, r_smooth]:
            col_means = np.nanmean(ratios, axis=0)
            for j in range(n_iso):
                ratios[np.isnan(ratios[:, j]), j] = col_means[j]

        ref_pval = ref_row['pvalue_hsic-ir'].values[0]
        p_obs = compute_hsic_pvalue(r_obs, K_sp, n)
        p_smooth = compute_hsic_pvalue(r_smooth, K_sp, n)

        ref_pvals.append(ref_pval)
        results_obs.append(p_obs)
        results_smooth.append(p_smooth)

    # Compute correlations
    log_ref = -np.log10(np.array(ref_pvals) + 1e-300)
    log_obs = -np.log10(np.array(results_obs) + 1e-300)
    log_smooth = -np.log10(np.array(results_smooth) + 1e-300)

    r_obs, _ = stats.spearmanr(log_ref, log_obs)
    r_smooth, _ = stats.spearmanr(log_ref, log_smooth)

    print(f"Number of genes tested: {len(ref_pvals)}")
    print(f"\nSpearman correlation (ratios_obs): {r_obs:.4f}")
    print(f"Spearman correlation (ratios_smoothed): {r_smooth:.4f}")

    # Mean absolute log difference
    mad_obs = np.mean(np.abs(log_obs - log_ref))
    mad_smooth = np.mean(np.abs(log_smooth - log_ref))
    print(f"\nMean abs log10 diff (ratios_obs): {mad_obs:.4f}")
    print(f"Mean abs log10 diff (ratios_smoothed): {mad_smooth:.4f}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Test if spot filtering (using only expressed spots) explains the discrepancy."""

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


def compute_hsic_pvalue(ratios, coords, k=6, rho=0.99):
    """Compute HSIC p-value for a gene."""
    n = len(coords)

    # Build spatial kernel and center it
    L_sp = build_spatial_kernel(coords, k=k, rho=rho)
    K_sp = center_kernel(L_sp)

    # Center ratios and build feature kernel
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

    return pvalue, hsic_trace


def main():
    # Load data
    sample_dir = Path('/data/parks34/projects/3csitme/data/zenodo_16905935/human_glioma_sr/MGH258')
    adata = ad.read_h5ad(sample_dir / 'iso.quant.h5ad')
    ref_df = pd.read_csv(sample_dir / 'sv_results.csv')

    coords_all = np.array(adata.obsm['spatial'])
    n_total = len(coords_all)

    # Get counts and ratios
    if 'counts' in adata.layers:
        counts_layer = adata.layers['counts']
        if hasattr(counts_layer, 'toarray'):
            counts_all = counts_layer.toarray()
        else:
            counts_all = np.array(counts_layer)
    else:
        if hasattr(adata.X, 'toarray'):
            counts_all = adata.X.toarray()
        else:
            counts_all = np.array(adata.X)

    ratios_layer = adata.layers['ratios_obs']
    if hasattr(ratios_layer, 'toarray'):
        ratios_all = ratios_layer.toarray()
    else:
        ratios_all = np.array(ratios_layer)

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
    test_genes = ['H3-3A', 'NCAM1', 'CALM1', 'SOX2-OT', 'PLP1', 'PTMA', 'GPM6B', 'RAC1']

    print(f"Total spots: {n_total}")
    print("\n" + "=" * 130)
    print("Comparing: all spots vs expressed spots only")
    print("=" * 130)
    header = f"{'gene':>10} {'n_iso':>6} {'pct_on':>8} {'n_on':>8} {'p_ref':>15} {'p_all':>15} {'p_on':>15} {'diff_all':>10} {'diff_on':>10}"
    print(header)
    print("-" * 130)

    for gene in test_genes:
        if gene not in gene_to_isoforms:
            continue

        iso_indices = gene_to_isoforms[gene]
        n_iso = len(iso_indices)

        # Get counts and ratios for this gene
        counts = counts_all[:, iso_indices]
        ratios = ratios_all[:, iso_indices].copy()

        # Gene totals
        gene_totals = counts.sum(axis=1)

        # Find spots with expression
        spots_on = gene_totals > 0
        n_on = spots_on.sum()
        pct_on = n_on / n_total

        # Handle NaN in ratios (for all spots)
        col_means = np.nanmean(ratios, axis=0)
        for j in range(n_iso):
            ratios[np.isnan(ratios[:, j]), j] = col_means[j]

        # Reference p-value and pct
        ref_row = ref_df[ref_df['gene'] == gene]
        ref_pval = ref_row['pvalue_hsic-ir'].values[0] if len(ref_row) > 0 else np.nan
        ref_pct = ref_row['pct_spot_on'].values[0] if len(ref_row) > 0 else np.nan

        # Method 1: Use all spots
        p_all, trace_all = compute_hsic_pvalue(ratios, coords_all, k=6, rho=0.99)

        # Method 2: Use only expressed spots
        if n_on > 50:  # Need enough spots
            coords_on = coords_all[spots_on]
            ratios_on = ratios[spots_on]
            p_on, trace_on = compute_hsic_pvalue(ratios_on, coords_on, k=6, rho=0.99)
        else:
            p_on = np.nan

        # Log differences
        diff_all = np.log10(p_all + 1e-300) - np.log10(ref_pval + 1e-300)
        diff_on = np.log10(p_on + 1e-300) - np.log10(ref_pval + 1e-300) if not np.isnan(p_on) else np.nan

        print(f"{gene:>10} {n_iso:>6} {pct_on:>8.4f} {n_on:>8} {ref_pval:>15.6e} {p_all:>15.6e} {p_on:>15.6e} {diff_all:>10.3f} {diff_on:>10.3f}")

        # Compare ref pct vs our pct
        if abs(pct_on - ref_pct) > 0.001:
            print(f"  WARNING: pct_spot_on mismatch: ours={pct_on:.6f}, ref={ref_pct:.6f}")

    print("\n" + "=" * 130)

    # Now let's check what reference parameters were used
    print("\nChecking reference parameters (from sv_results.csv for H3-3A):")
    row = ref_df[ref_df['gene'] == 'H3-3A']
    print(f"  n_iso: {row['n_iso'].values[0]}")
    print(f"  pct_spot_on: {row['pct_spot_on'].values[0]:.6f}")
    print(f"  count_avg: {row['count_avg'].values[0]:.4f}")
    print(f"  count_std: {row['count_std'].values[0]:.4f}")
    print(f"  perplexity: {row['perplexity'].values[0]:.4f}")
    print(f"  major_ratio_avg: {row['major_ratio_avg'].values[0]:.4f}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Test different eigenvalue computation methods to match official SPLISOSM."""

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


def main():
    # Load data
    sample_dir = Path('/data/parks34/projects/3csitme/data/zenodo_16905935/human_glioma_sr/MGH258')
    adata = ad.read_h5ad(sample_dir / 'iso.quant.h5ad')
    ref_df = pd.read_csv(sample_dir / 'sv_results.csv')

    # Extract coords
    coords = np.array(adata.obsm['spatial'])
    n = len(coords)
    k = 6
    rho = 0.99

    print(f"Sample size: n = {n}")

    # Build spatial kernel
    W = kneighbors_graph(coords, n_neighbors=k, mode='connectivity')
    W = W + W.T
    W.data = np.ones_like(W.data)
    degrees = np.array(W.sum(axis=1)).flatten()
    D = sparse.diags(degrees, format='csr')
    Q = D - rho * W + sparse.eye(n, format='csr') * 1e-6
    lu = splu(Q.tocsc())
    L_sp = np.zeros((n, n))
    I = np.eye(n)
    for i in range(n):
        L_sp[:, i] = lu.solve(I[:, i])
    L_sp = (L_sp + L_sp.T) / 2

    # Center spatial kernel (this is K_sp in official notation)
    K_sp = center_kernel(L_sp)
    eig_sp = np.linalg.eigvalsh(K_sp)
    eig_sp_filtered = eig_sp[np.abs(eig_sp) > 1e-5]

    print(f"Spatial kernel eigenvalues: {len(eig_sp_filtered)} non-trivial")

    # Get gene info
    gene_col = 'Gene' if 'Gene' in adata.var.columns else 'gene_symbol'
    ratios_all = np.array(adata.layers['ratios_obs'])

    # Build gene mapping
    gene_to_isoforms = {}
    for idx, gene in enumerate(adata.var[gene_col]):
        if gene not in gene_to_isoforms:
            gene_to_isoforms[gene] = []
        gene_to_isoforms[gene].append(idx)
    gene_to_isoforms = {g: v for g, v in gene_to_isoforms.items() if len(v) >= 2}

    # Test genes
    test_genes = ['H3-3A', 'NCAM1', 'CALM1', 'SOX2-OT', 'PLP1', 'PTMA', 'GPM6B', 'RAC1']

    print("\n" + "=" * 120)
    print("Testing different eigenvalue computation methods:")
    print("=" * 120)
    header = f"{'gene':>10} {'p_ref':>15} {'p_yTy':>15} {'p_Ky':>15} {'p_Ky_ctr':>15} {'diff_yTy':>10} {'diff_Ky':>10} {'diff_ctr':>10}"
    print(header)
    print("-" * 120)

    for gene in test_genes:
        if gene not in gene_to_isoforms:
            continue

        iso_indices = gene_to_isoforms[gene]
        ratios = ratios_all[:, iso_indices].copy()
        p = len(iso_indices)

        # Handle NaN
        col_means = np.nanmean(ratios, axis=0)
        for j in range(p):
            ratios[np.isnan(ratios[:, j]), j] = col_means[j]

        # Center ratios (column-wise, like official)
        y = ratios - ratios.mean(axis=0, keepdims=True)

        # Compute HSIC trace
        # tr(K_sp @ (y @ y.T)) = tr(y.T @ K_sp @ y)
        # This is scalar: sum of diagonal of (p x p) matrix
        hsic_trace = np.trace(y.T @ K_sp @ y)

        # Reference p-value
        ref_row = ref_df[ref_df['gene'] == gene]
        ref_pval = ref_row['pvalue_hsic-ir'].values[0] if len(ref_row) > 0 else np.nan

        # Method 1: Official style - eigenvalues of y.T @ y (p x p)
        # Note: y.T @ y is the sample covariance (uncentered gram matrix of columns)
        gram_y = y.T @ y  # (p, p)
        eig_y_yTy = np.linalg.eigvalsh(gram_y)
        eig_y_yTy_filtered = eig_y_yTy[np.abs(eig_y_yTy) > 1e-10]

        lambda_xy_yTy = np.outer(eig_sp_filtered, eig_y_yTy_filtered).flatten()
        p_yTy = liu_sf(hsic_trace * n, lambda_xy_yTy)

        # Method 2: Our style - eigenvalues of K_y = y @ y.T (n x n), NOT centered
        K_y = y @ y.T  # (n, n)
        eig_y_Ky = np.linalg.eigvalsh(K_y)
        eig_y_Ky_filtered = eig_y_Ky[np.abs(eig_y_Ky) > 1e-10]

        lambda_xy_Ky = np.outer(eig_sp_filtered, eig_y_Ky_filtered).flatten()
        p_Ky = liu_sf(hsic_trace * n, lambda_xy_Ky)

        # Method 3: Our current implementation - eigenvalues of centered K_y
        K_y_centered = center_kernel(K_y)
        eig_y_Ky_ctr = np.linalg.eigvalsh(K_y_centered)
        eig_y_Ky_ctr_filtered = eig_y_Ky_ctr[np.abs(eig_y_Ky_ctr) > 1e-10]

        lambda_xy_Ky_ctr = np.outer(eig_sp_filtered, eig_y_Ky_ctr_filtered).flatten()
        p_Ky_ctr = liu_sf(hsic_trace * n, lambda_xy_Ky_ctr)

        # Log differences
        diff_yTy = np.log10(p_yTy + 1e-300) - np.log10(ref_pval + 1e-300)
        diff_Ky = np.log10(p_Ky + 1e-300) - np.log10(ref_pval + 1e-300)
        diff_ctr = np.log10(p_Ky_ctr + 1e-300) - np.log10(ref_pval + 1e-300)

        print(f"{gene:>10} {ref_pval:>15.6e} {p_yTy:>15.6e} {p_Ky:>15.6e} {p_Ky_ctr:>15.6e} {diff_yTy:>10.3f} {diff_Ky:>10.3f} {diff_ctr:>10.3f}")

        # Debug info for the first gene
        if gene == 'H3-3A':
            print(f"\n  Debug for {gene}:")
            print(f"    n = {n}, p = {p}")
            print(f"    hsic_trace = {hsic_trace:.6e}")
            print(f"    hsic_trace * n = {hsic_trace * n:.6e}")
            print(f"    len(eig_sp_filtered) = {len(eig_sp_filtered)}")
            print(f"    len(eig_y_yTy_filtered) = {len(eig_y_yTy_filtered)}")
            print(f"    len(eig_y_Ky_filtered) = {len(eig_y_Ky_filtered)}")
            print(f"    len(eig_y_Ky_ctr_filtered) = {len(eig_y_Ky_ctr_filtered)}")
            print(f"    sum(eig_sp_filtered) = {np.sum(eig_sp_filtered):.6e}")
            print(f"    sum(eig_y_yTy_filtered) = {np.sum(eig_y_yTy_filtered):.6e}")
            print(f"    sum(eig_y_Ky_filtered) = {np.sum(eig_y_Ky_filtered):.6e}")
            print(f"    sum(eig_y_Ky_ctr_filtered) = {np.sum(eig_y_Ky_ctr_filtered):.6e}")
            print()

    print("\nKey observations:")
    print("- p_yTy: eigenvalues from y.T @ y (official approach)")
    print("- p_Ky: eigenvalues from K_y = y @ y.T (same non-zero values as y.T @ y)")
    print("- p_Ky_ctr: eigenvalues from centered K_y (our current approach)")


if __name__ == '__main__':
    main()

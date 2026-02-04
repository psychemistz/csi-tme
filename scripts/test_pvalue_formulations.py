#!/usr/bin/env python3
"""Test different p-value formulations to match official SPLISOSM exactly."""

import numpy as np
from scipy import stats
from scipy.stats import ncx2


def liu_sf_official(t, eigenvalues, dofs=None, deltas=None):
    """Liu's approximation exactly as implemented in official SPLISOSM.

    Based on likelihood.py in official repo.
    """
    eigenvalues = np.asarray(eigenvalues).flatten()
    eigenvalues = eigenvalues[np.abs(eigenvalues) > 1e-10]

    if len(eigenvalues) == 0:
        return 1.0

    n_terms = len(eigenvalues)

    if dofs is None:
        dofs = np.ones(n_terms)
    if deltas is None:
        deltas = np.zeros(n_terms)

    # Compute cumulants
    c1 = np.sum(eigenvalues * (dofs + deltas))
    c2 = np.sum(eigenvalues**2 * (dofs + 2 * deltas))
    c3 = np.sum(eigenvalues**3 * (dofs + 3 * deltas))
    c4 = np.sum(eigenvalues**4 * (dofs + 4 * deltas))

    if c2 <= 0:
        return 1.0

    s1 = c3 / (c2 ** 1.5)
    s2 = c4 / (c2 ** 2)

    # Liu's approximation
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
        pvalue = ncx2.sf(t_final, df_approx, delta_approx)

    return np.clip(pvalue, 1e-300, 1.0)


def compute_pvalue_official_style(hsic_trace, eig_sp, eig_y, n):
    """Compute p-value the way official SPLISOSM does.

    Official approach:
    1. No eigenvalue scaling
    2. Pass hsic_trace * n to liu_sf
    3. Composite eigenvalues = outer product of spatial and feature eigenvalues
    """
    lambda_xy = np.outer(eig_sp, eig_y).flatten()
    lambda_xy = lambda_xy[np.abs(lambda_xy) > 1e-10]
    return liu_sf_official(hsic_trace * n, lambda_xy)


def compute_pvalue_local_style(hsic_statistic, eig_sp, eig_y, n):
    """Compute p-value the way our local implementation does.

    Local approach:
    1. Scale eigenvalues by 1/n
    2. Pass n * hsic_statistic to liu_sf
    3. hsic_statistic = hsic_trace / (n-1)^2
    """
    eig_sp_scaled = eig_sp / n
    eig_y_scaled = eig_y / n
    lambda_xy = np.outer(eig_sp_scaled, eig_y_scaled).flatten()
    lambda_xy = lambda_xy[np.abs(lambda_xy) > 1e-10]
    return liu_sf_official(n * hsic_statistic, lambda_xy)


def test_equivalence():
    """Test if the two formulations are equivalent."""
    np.random.seed(42)

    # Generate test data
    n = 2212  # Like MGH258 sample
    p = 3  # Number of isoforms

    # Simulated spatial kernel eigenvalues (positive)
    eig_sp = np.abs(np.random.randn(n)) + 0.1
    eig_sp = eig_sp[eig_sp > 1e-5]

    # Simulated feature kernel eigenvalues
    eig_y = np.abs(np.random.randn(min(p, n))) + 0.1
    eig_y = eig_y[eig_y > 1e-5]

    # Test different hsic_trace values
    test_traces = [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]

    print("Testing p-value formulations:")
    print("=" * 80)
    print(f"{'hsic_trace':>15} {'hsic_stat':>15} {'p_official':>15} {'p_local':>15} {'log_diff':>15}")
    print("=" * 80)

    for hsic_trace in test_traces:
        hsic_statistic = hsic_trace / ((n - 1) ** 2)

        p_official = compute_pvalue_official_style(hsic_trace, eig_sp, eig_y, n)
        p_local = compute_pvalue_local_style(hsic_statistic, eig_sp, eig_y, n)

        log_diff = np.log10(p_local + 1e-300) - np.log10(p_official + 1e-300)

        print(f"{hsic_trace:>15.6f} {hsic_statistic:>15.2e} {p_official:>15.6e} {p_local:>15.6e} {log_diff:>15.4f}")

    print("\n" + "=" * 80)
    print("\nAnalysis:")

    # The difference comes from how eigenvalues are scaled
    # Official: lambda_xy (unscaled), statistic = trace * n
    # Local: lambda_xy / n^2, statistic = trace * n / (n-1)^2

    # For the null distribution:
    # E[T_official] = sum(lambda_xy) = sum_sp * sum_y
    # E[T_local] = sum(lambda_xy / n^2) = sum_sp * sum_y / n^2

    # Var[T_official] = 2 * sum(lambda_xy^2)
    # Var[T_local] = 2 * sum((lambda_xy / n^2)^2) = 2 * sum(lambda_xy^2) / n^4

    # For standardization:
    # official_z = (trace*n - E_off) / sqrt(Var_off)
    # local_z = (trace*n/(n-1)^2 - E_loc) / sqrt(Var_loc)
    #         = (trace*n/(n-1)^2 - E_off/n^2) / (sqrt(Var_off)/n^2)
    #         = n^2 * (trace*n/(n-1)^2 - E_off/n^2) / sqrt(Var_off)
    #         = (trace*n^3/(n-1)^2 - E_off) / sqrt(Var_off)

    # So local_z ≈ (trace*n - E_off) / sqrt(Var_off) only if n^3/(n-1)^2 ≈ n
    # n^3/(n-1)^2 = n^3/(n^2 - 2n + 1) ≈ n for large n
    # But for n=2212: n^3/(n-1)^2 ≈ 2212^3/2211^2 ≈ 2212.9 vs n=2212

    print(f"n = {n}")
    print(f"n^3/(n-1)^2 = {n**3 / (n-1)**2:.4f}")
    print(f"Ratio to n = {n**3 / ((n-1)**2 * n):.6f}")
    print("\nThe two formulations should give nearly identical results for large n.")


def analyze_real_discrepancy():
    """Analyze the real discrepancy on actual data."""
    import sys
    sys.path.insert(0, '/vf/users/parks34/projects/3csitme')

    import anndata as ad
    import pandas as pd
    from pathlib import Path
    from sklearn.neighbors import kneighbors_graph
    from scipy import sparse
    from scipy.sparse.linalg import splu

    from splisosm.statistics.hsic import center_kernel

    # Load data
    sample_dir = Path('/data/parks34/projects/3csitme/data/zenodo_16905935/human_glioma_sr/MGH258')
    adata = ad.read_h5ad(sample_dir / 'iso.quant.h5ad')
    ref_df = pd.read_csv(sample_dir / 'sv_results.csv')

    # Extract coords
    coords = np.array(adata.obsm['spatial'])
    n = len(coords)
    k = 6
    rho = 0.99

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

    # Center spatial kernel
    L_sp_centered = center_kernel(L_sp)
    eig_sp = np.linalg.eigvalsh(L_sp_centered)
    eig_sp = eig_sp[np.abs(eig_sp) > 1e-10]

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

    # Test on a few genes
    test_genes = ['H3-3A', 'NCAM1', 'CALM1', 'SOX2-OT', 'PLP1']

    print("\nComparing formulations on real data:")
    print("=" * 100)
    print(f"{'gene':>12} {'n_iso':>6} {'ref_pval':>15} {'official':>15} {'local':>15} {'no_scale':>15}")
    print("=" * 100)

    for gene in test_genes:
        if gene not in gene_to_isoforms:
            continue

        iso_indices = gene_to_isoforms[gene]
        ratios = ratios_all[:, iso_indices].copy()

        # Handle NaN
        col_means = np.nanmean(ratios, axis=0)
        for j in range(ratios.shape[1]):
            ratios[np.isnan(ratios[:, j]), j] = col_means[j]

        # Center ratios
        y = ratios - ratios.mean(axis=0, keepdims=True)

        # Feature kernel and its eigenvalues
        K_y = y @ y.T
        K_y_centered = center_kernel(K_y)
        eig_y = np.linalg.eigvalsh(K_y_centered)
        eig_y = eig_y[np.abs(eig_y) > 1e-10]

        # HSIC trace
        hsic_trace = np.sum(L_sp_centered * K_y_centered)
        hsic_statistic = hsic_trace / ((n - 1) ** 2)

        # P-values with different approaches
        p_official = compute_pvalue_official_style(hsic_trace, eig_sp, eig_y, n)
        p_local = compute_pvalue_local_style(hsic_statistic, eig_sp, eig_y, n)

        # No eigenvalue scaling, but use hsic_statistic * n
        lambda_xy_no_scale = np.outer(eig_sp, eig_y).flatten()
        p_no_scale = liu_sf_official(n * hsic_statistic, lambda_xy_no_scale)

        # Reference
        ref_row = ref_df[ref_df['gene'] == gene]
        ref_pval = ref_row['pvalue_hsic-ir'].values[0] if len(ref_row) > 0 else np.nan

        print(f"{gene:>12} {len(iso_indices):>6} {ref_pval:>15.6e} {p_official:>15.6e} {p_local:>15.6e} {p_no_scale:>15.6e}")

    print("\nConclusion:")
    print("- 'official' = hsic_trace * n, unscaled eigenvalues")
    print("- 'local' = n * hsic_statistic, eigenvalues / n^2")
    print("- 'no_scale' = n * hsic_statistic, unscaled eigenvalues")


if __name__ == '__main__':
    test_equivalence()
    print("\n" + "=" * 100 + "\n")
    analyze_real_discrepancy()

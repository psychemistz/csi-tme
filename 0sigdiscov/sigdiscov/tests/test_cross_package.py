#!/usr/bin/env python
"""
Cross-package equivalence test: verify compute_morans.py, sigdiscovpy, and
sigdiscov (R) produce algorithmically and numerically identical results
against precomputed reference data.

Tests three implementations:
  1. sigdiscovpy Python package (low-level API, replicating compute_morans.py algorithm)
  2. sigdiscov R package (via subprocess, using its own batch function)
  3. sigdiscov R package (via subprocess, using low-level functions matching Python algorithm)

Usage:
    python tests/test_cross_package.py [--gpu 0] [--n-factors 20] [--n-radii 3] [--ring-width 20]
"""

import sys
import os
import subprocess
import tempfile
import shutil
from pathlib import Path
import numpy as np
import h5py
import argparse
import time
import json

# Paths
SCRIPT_DIR = Path(__file__).parent.parent
DATA_PATH = SCRIPT_DIR / 'dataset' / 'coad_spatial.h5ad'
FACTOR_PATH = SCRIPT_DIR / 'dataset' / 'db' / 'AllSigFilteredBy_MoranI_TCGA_ICGC_0.25_ds3.tsv.gz'
REF_DIR = Path(__file__).parent / 'reference_data'

# Test parameters matching reference data
CELL_PAIRS = [
    ("tumor", "CAF"),
    ("tumor", "T CD8 memory"),
    ("CAF", "macrophage"),
    ("macrophage", "CAF"),
    ("Treg", "macrophage"),
]
ALL_RADII = [20, 100, 260, 500, 1000, 2000, 3000, 5000]


def load_reference(quantile_label: str):
    """Load reference H5 data for a given quantile label."""
    ref_path = REF_DIR / f'ref_{quantile_label}.h5'
    assert ref_path.exists(), f"Reference not found: {ref_path}"

    with h5py.File(ref_path, 'r') as f:
        data = {
            'matrix': f['morans_i_matrices'][:],
            'pairs': [x.decode() for x in f['pairs'][:]],
            'radii': f['radii'][:].tolist(),
            'factor_genes': [x.decode() for x in f['factor_genes'][:]],
            'target_genes': [x.decode() for x in f['target_genes'][:]],
        }
    return data


def compare_matrices(test_mat, ref_mat, label: str, verbose: bool = True):
    """Compare two matrices and return statistics."""
    # Handle NaN: treat NaN in test as 0 (matching compute_morans.py behavior)
    test_clean = np.where(np.isnan(test_mat), 0.0, test_mat)
    ref_clean = np.where(np.isnan(ref_mat), 0.0, ref_mat)

    abs_diff = np.abs(test_clean - ref_clean)
    total = abs_diff.size

    stats = {
        'label': label,
        'total_values': total,
        'max_abs_diff': float(abs_diff.max()),
        'mean_abs_diff': float(abs_diff.mean()),
        'rmse': float(np.sqrt(np.mean(abs_diff ** 2))),
        'exact_match_pct': 100 * float(np.sum(abs_diff == 0)) / total,
        'close_1e6_pct': 100 * float(np.sum(abs_diff < 1e-6)) / total,
        'close_1e4_pct': 100 * float(np.sum(abs_diff < 1e-4)) / total,
        'close_1e2_pct': 100 * float(np.sum(abs_diff < 1e-2)) / total,
    }

    # Max relative diff for non-trivial values
    nontrivial = np.abs(ref_clean) > 1e-4
    if nontrivial.any():
        rel_diff = abs_diff[nontrivial] / np.abs(ref_clean[nontrivial])
        stats['max_rel_diff'] = float(rel_diff.max())
        stats['mean_rel_diff'] = float(rel_diff.mean())
    else:
        stats['max_rel_diff'] = 0.0
        stats['mean_rel_diff'] = 0.0

    if verbose:
        print(f"\n  === {label} ===")
        print(f"  Values compared:   {stats['total_values']:,}")
        print(f"  Max abs diff:      {stats['max_abs_diff']:.2e}")
        print(f"  Mean abs diff:     {stats['mean_abs_diff']:.2e}")
        print(f"  RMSE:              {stats['rmse']:.2e}")
        print(f"  Max rel diff:      {stats['max_rel_diff']:.2e}")
        print(f"  Exact matches:     {stats['exact_match_pct']:.1f}%")
        print(f"  Within 1e-6:       {stats['close_1e6_pct']:.1f}%")
        print(f"  Within 1e-4:       {stats['close_1e4_pct']:.1f}%")

    return stats


# ============================================================================
# TEST 1: sigdiscovpy (Python package, low-level API)
# ============================================================================

def test_sigdiscovpy(ref_data, quantile_val, test_radii, ring_width, n_factors_limit=None):
    """Test sigdiscovpy by manually replicating compute_morans.py algorithm."""
    print("\n" + "=" * 70)
    print("TEST: sigdiscovpy (Python package)")
    print("=" * 70)

    try:
        import anndata as ad
        from scipy import sparse as sp_sparse
        from sigdiscovpy.core.weights import (
            _compute_distances_chunked, row_normalize_weights
        )
        from sigdiscovpy.core.normalization import standardize_matrix
        from sigdiscovpy.core.spatial_lag import compute_spatial_lag_batch
        from sigdiscovpy.core.metrics import compute_metric_batch
    except ImportError as e:
        print(f"  SKIP: sigdiscovpy not available: {e}")
        return None

    ref_pairs = ref_data['pairs']
    ref_radii = ref_data['radii']
    ref_factors = ref_data['factor_genes']
    ref_targets = ref_data['target_genes']

    # Limit factors for speed
    if n_factors_limit and n_factors_limit < len(ref_factors):
        factor_subset = ref_factors[:n_factors_limit]
        factor_slice = slice(0, n_factors_limit)
    else:
        factor_subset = ref_factors
        factor_slice = slice(None)
        n_factors_limit = len(ref_factors)

    # Select radii to test
    radii_indices = [ref_radii.index(r) for r in test_radii if r in ref_radii]
    test_radii_actual = [ref_radii[i] for i in radii_indices]

    print(f"  Pairs:   {len(ref_pairs)}")
    print(f"  Radii:   {test_radii_actual}")
    print(f"  Factors: {n_factors_limit} (of {len(ref_factors)})")
    print(f"  Targets: {len(ref_targets)}")
    print(f"  Quantile: {quantile_val}")
    print(f"  Ring width: {ring_width}")

    # Load data (same as compute_morans.py)
    print("  Loading data...")
    t0 = time.time()
    adata = ad.read_h5ad(DATA_PATH)
    print(f"  Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes ({time.time()-t0:.1f}s)")

    # Filter genes (matching compute_morans.py)
    n_cells = adata.shape[0]
    if sp_sparse.issparse(adata.X):
        X_csc = adata.X.tocsc()
        gene_counts = np.diff(X_csc.indptr)
    else:
        gene_counts = (adata.X > 0).sum(axis=0)
    min_cells = max(10, n_cells * 0.001)
    keep_genes = gene_counts >= min_cells
    adata_filtered = adata[:, keep_genes].copy()
    gene_names = adata_filtered.var_names.tolist()
    n_genes = len(gene_names)
    print(f"  Genes after filter: {n_genes}")

    # Verify gene lists match
    assert gene_names == ref_targets, \
        f"Gene list mismatch: {len(gene_names)} vs {len(ref_targets)}"

    # Find factor indices
    factor_indices = [gene_names.index(f) for f in factor_subset if f in gene_names]
    factor_names_found = [gene_names[i] for i in factor_indices]
    assert factor_names_found == list(factor_subset), "Factor gene mismatch"
    n_factors = len(factor_indices)

    # Dense + normalize (matching compute_morans.py: axis=0, population std)
    print("  Normalizing...")
    if sp_sparse.issparse(adata_filtered.X):
        X_dense = adata_filtered.X.toarray().astype(np.float32)
    else:
        X_dense = adata_filtered.X.astype(np.float32)

    # Compute mean/std per gene (axis=0 = across cells)
    X_mean = X_dense.mean(axis=0, keepdims=True)
    X_std = X_dense.std(axis=0, keepdims=True)
    X_normalized = (X_dense - X_mean) / (X_std + 1e-10)

    # Factor data
    X_factors_raw = X_dense[:, factor_indices]  # raw factor expression, ALL cells
    X_factors_norm = X_normalized[:, factor_indices]  # normalized factor expression, ALL cells

    # Quantile thresholds (on raw factor expression, ALL cells)
    quantile_thresholds = np.quantile(X_factors_raw, quantile_val, axis=0)

    # Cache cell data by type
    cell_coords_cache = {}
    cell_expr_cache = {}
    cell_factor_norm_cache = {}
    for pair_str in ref_pairs:
        sender, receiver = pair_str.split('->')
        for ct in [sender, receiver]:
            if ct not in cell_coords_cache:
                mask = (adata_filtered.obs['cell_type'] == ct).values
                cell_coords_cache[ct] = adata_filtered.obsm['spatial'][mask].astype(np.float64)
                cell_expr_cache[ct] = X_normalized[mask, :]
                cell_factor_norm_cache[ct] = X_factors_norm[mask, :]

    # Compute I_ND for each pair/radius
    all_results = {}
    total_t0 = time.time()

    for pair_str in ref_pairs:
        sender, receiver = pair_str.split('->')
        sender_coords = cell_coords_cache[sender]
        receiver_coords = cell_coords_cache[receiver]
        sender_factor_norm = cell_factor_norm_cache[sender]
        receiver_expr_norm = cell_expr_cache[receiver]
        n_senders = len(sender_coords)

        for ri, radius in zip(radii_indices, test_radii_actual):
            inner_radius = max(0, radius - ring_width)

            t1 = time.time()
            # Build weight matrix in FLOAT32 matching compute_morans.py GPU computation
            sigma = radius / 3.0
            gaussian_factor = np.float32(-1.0 / (2.0 * sigma * sigma))
            radius_sq = np.float32(radius * radius)
            inner_radius_sq = np.float32(inner_radius * inner_radius)

            n_receivers = len(receiver_coords)
            sender_f32 = sender_coords.astype(np.float32)
            receiver_f32 = receiver_coords.astype(np.float32)

            w_rows, w_cols, w_vals = [], [], []
            chunk_sz = 500
            for cs in range(0, n_senders, chunk_sz):
                ce = min(cs + chunk_sz, n_senders)
                dx = sender_f32[cs:ce, 0:1] - receiver_f32[:, 0].T
                dy = sender_f32[cs:ce, 1:2] - receiver_f32[:, 1].T
                dsq = dx * dx + dy * dy

                m = dsq <= radius_sq
                if inner_radius > 0:
                    m = m & (dsq > inner_radius_sq)

                wc = np.exp(dsq * gaussian_factor) * m
                wc[wc <= 1e-6] = 0

                # Inline row normalization (matching compute_morans.py and R C++)
                chunk_row_sums = wc.sum(axis=1, keepdims=True)
                chunk_row_sums = np.where(chunk_row_sums > 0, chunk_row_sums, np.float32(1.0))
                wc /= chunk_row_sums

                nz = wc > 0
                if nz.any():
                    ri, ci = np.where(nz)
                    w_rows.append(ri + cs)
                    w_cols.append(ci)
                    w_vals.append(wc[nz].astype(np.float32))

                del dx, dy, dsq, m, wc

            if not w_rows:
                all_results[(pair_str, radius)] = np.zeros((n_factors, n_genes), dtype=np.float32)
                print(f"    {pair_str:30s} R={radius:5.0f}  (no connections)")
                continue

            all_r = np.concatenate(w_rows)
            all_c = np.concatenate(w_cols)
            all_w = np.concatenate(w_vals)

            W = sp_sparse.csr_matrix(
                (all_w, (all_r, all_c)),
                shape=(n_senders, n_receivers), dtype=np.float32,
            )

            # Spatial lags: W @ receiver_expr_norm → (n_senders x n_genes)
            # Need receiver_expr in (n_receivers, n_genes) format
            spatial_lags = W @ receiver_expr_norm  # scipy sparse @ numpy → numpy

            # Compute I_ND for each factor
            result_matrix = np.zeros((n_factors, n_genes), dtype=np.float32)

            for f_idx in range(n_factors):
                # Quantile mask: compare NORMALIZED sender expression to RAW threshold
                # (This matches compute_morans.py's "bug")
                mask = sender_factor_norm[:, f_idx] > quantile_thresholds[f_idx]

                if mask.sum() < 10:
                    continue

                high_indices = np.where(mask)[0]
                factor_high = sender_factor_norm[high_indices, f_idx].astype(np.float32)
                lags_high = spatial_lags[high_indices, :].astype(np.float32)

                # I_ND = dot(factor_high/||factor_high||, lags_high) / ||lags_high_columns||
                # Use float32 matching compute_morans.py GPU computation
                factor_norm_val = np.sqrt(np.sum(factor_high ** 2))
                if factor_norm_val < 1e-10:
                    continue

                factor_normalized = factor_high / factor_norm_val
                correlations = factor_normalized @ lags_high  # (n_genes,)
                spatial_norms = np.sqrt(np.sum(lags_high ** 2, axis=0))  # (n_genes,)

                valid = spatial_norms > 0
                result_matrix[f_idx, valid] = correlations[valid] / spatial_norms[valid]

            all_results[(pair_str, radius)] = result_matrix
            elapsed = time.time() - t1
            print(f"    {pair_str:30s} R={radius:5.0f}  ({elapsed:.1f}s)")

    total_time = time.time() - total_t0
    print(f"  Total computation: {total_time:.1f}s")

    # Compare to reference
    print("\n  Comparing to reference...")
    all_stats = []
    for pair_str in ref_pairs:
        pi = ref_pairs.index(pair_str)
        for ri, radius in zip(radii_indices, test_radii_actual):
            test_mat = all_results[(pair_str, radius)]
            ref_mat = ref_data['matrix'][pi, ri, factor_slice, :]

            label = f"sigdiscovpy {pair_str} R={radius}"
            stats = compare_matrices(test_mat, ref_mat, label, verbose=False)
            all_stats.append(stats)

    # Summary
    max_abs_all = max(s['max_abs_diff'] for s in all_stats)
    max_rel_all = max(s['max_rel_diff'] for s in all_stats)
    mean_close = np.mean([s['close_1e6_pct'] for s in all_stats])

    print(f"\n  SIGDISCOVPY SUMMARY:")
    print(f"  Slices tested:     {len(all_stats)}")
    print(f"  Max abs diff:      {max_abs_all:.2e}")
    print(f"  Max rel diff:      {max_rel_all:.2e}")
    print(f"  Mean within 1e-6:  {mean_close:.1f}%")

    # Show worst slices
    worst = sorted(all_stats, key=lambda s: s['max_abs_diff'], reverse=True)[:5]
    print(f"\n  Worst 5 slices:")
    for s in worst:
        print(f"    {s['label']:50s}  max_abs={s['max_abs_diff']:.2e}  "
              f"within_1e-4={s['close_1e4_pct']:.1f}%")

    return {
        'package': 'sigdiscovpy',
        'stats': all_stats,
        'max_abs_diff': max_abs_all,
        'max_rel_diff': max_rel_all,
        'mean_close_1e6_pct': mean_close,
        'all_slices': all_results,
    }


# ============================================================================
# TEST 2: sigdiscov R package (low-level functions matching Python algorithm)
# ============================================================================

def _export_data_for_r(tmpdir, ref_data, factor_subset, n_factors_limit):
    """Export data from Python to HDF5 for R consumption.

    All numeric data is exported in float32 to match compute_morans.py GPU precision.
    Also exports pre-normalized expression for the R low-level test.
    """
    import anndata as ad
    from scipy import sparse as sp_sparse

    export_path = tmpdir / 'data_for_r.h5'
    if export_path.exists():
        return export_path

    print("  Exporting data for R...")
    t0 = time.time()
    adata = ad.read_h5ad(DATA_PATH)

    # Filter genes
    n_cells = adata.shape[0]
    if sp_sparse.issparse(adata.X):
        X_csc = adata.X.tocsc()
        gene_counts = np.diff(X_csc.indptr)
    else:
        gene_counts = (adata.X > 0).sum(axis=0)
    min_cells = max(10, n_cells * 0.001)
    keep_genes = gene_counts >= min_cells
    adata_f = adata[:, keep_genes].copy()
    gene_names = adata_f.var_names.tolist()

    # Float32 matching compute_morans.py GPU precision
    if sp_sparse.issparse(adata_f.X):
        X_dense = adata_f.X.toarray().astype(np.float32)
    else:
        X_dense = adata_f.X.astype(np.float32)

    coords = adata_f.obsm['spatial'].astype(np.float32)
    cell_types = adata_f.obs['cell_type'].values.astype(str)

    # Normalize in float32 (matching compute_morans.py)
    X_mean = X_dense.mean(axis=0)
    X_std = X_dense.std(axis=0)
    X_norm = (X_dense - X_mean) / (X_std + np.float32(1e-10))

    # rhdf5 reads in column-major (R) order, so transpose to get cells x genes in R
    with h5py.File(export_path, 'w') as f:
        f.create_dataset('expr', data=X_dense.T)  # raw float32
        f.create_dataset('expr_norm', data=X_norm.T)  # normalized float32
        f.create_dataset('coords', data=coords.T)  # float32
        f.create_dataset('cell_types', data=np.array(cell_types, dtype='S'))
        f.create_dataset('gene_names', data=np.array(gene_names, dtype='S'))

    print(f"  Exported: {X_dense.shape} ({time.time()-t0:.1f}s)")
    return export_path


def test_sigdiscov_r(ref_data, quantile_val, test_radii, ring_width,
                     n_factors_limit=None, tmpdir=None):
    """Test sigdiscov R package using low-level functions matching Python algorithm."""
    print("\n" + "=" * 70)
    print("TEST: sigdiscov R package (low-level, matching Python algorithm)")
    print("=" * 70)

    ref_pairs = ref_data['pairs']
    ref_radii = ref_data['radii']
    ref_factors = ref_data['factor_genes']
    ref_targets = ref_data['target_genes']

    if n_factors_limit and n_factors_limit < len(ref_factors):
        factor_subset = ref_factors[:n_factors_limit]
        factor_slice = slice(0, n_factors_limit)
        n_factors_test = n_factors_limit
    else:
        factor_subset = ref_factors
        factor_slice = slice(None)
        n_factors_test = len(ref_factors)

    radii_indices = [ref_radii.index(r) for r in test_radii if r in ref_radii]
    test_radii_actual = [ref_radii[i] for i in radii_indices]

    print(f"  Pairs:   {len(ref_pairs)}")
    print(f"  Radii:   {test_radii_actual}")
    print(f"  Factors: {n_factors_test}")
    print(f"  Quantile: {quantile_val}")
    print(f"  Ring width: {ring_width}")

    if tmpdir is None:
        tmpdir = Path(tempfile.mkdtemp(prefix='r_test_'))

    # Export data from Python for R
    data_path = _export_data_for_r(tmpdir, ref_data, factor_subset, n_factors_limit)

    r_output = tmpdir / 'r_results.h5'
    r_script = tmpdir / 'test_lowlevel.R'

    # Convert Python lists to R vector syntax
    def to_r_char_vec(lst):
        return 'c(' + ', '.join(f'"{x}"' for x in lst) + ')'
    def to_r_num_vec(lst):
        return 'c(' + ', '.join(str(x) for x in lst) + ')'

    r_code = f'''
library(sigdiscov)
library(rhdf5)

cat("Loading exported data (float32 precision)...\\n")
t0 <- proc.time()
expr_raw <- h5read("{data_path}", "expr")       # cells x genes (float32 precision)
coords <- h5read("{data_path}", "coords")        # cells x 2 (float32)
cell_types <- as.character(h5read("{data_path}", "cell_types"))
gene_names <- as.character(h5read("{data_path}", "gene_names"))
H5close()

n_cells <- nrow(expr_raw)
n_genes <- ncol(expr_raw)
cat(sprintf("Loaded: %d cells, %d genes (%.1fs)\\n",
    n_cells, n_genes, (proc.time() - t0)[3]))

# Normalize in R: population std (ddof=0), matching numpy's X.mean()/X.std()
cat("Normalizing (population std, matching numpy)...\\n")
t_norm <- proc.time()
expr_norm <- matrix(0, nrow = n_cells, ncol = n_genes)
for (g in seq_len(n_genes)) {{
    x <- expr_raw[, g]
    m <- mean(x)
    s <- sqrt(sum((x - m)^2) / n_cells)
    if (s > 1e-10) {{
        expr_norm[, g] <- (x - m) / s
    }}
}}
cat(sprintf("Normalization: %.1fs\\n", (proc.time() - t_norm)[3]))

# Factor setup
factor_genes <- {to_r_char_vec(factor_subset)}
factor_idx <- match(factor_genes, gene_names)
valid <- !is.na(factor_idx)
factor_idx <- factor_idx[valid]
factor_genes <- factor_genes[valid]
n_factors <- length(factor_idx)
cat(sprintf("Factors found: %d\\n", n_factors))

# Quantile thresholds (on RAW float32 factor expression, ALL cells)
quantile_prob <- {quantile_val}
quant_thresholds <- numeric(n_factors)
for (f in seq_len(n_factors)) {{
    quant_thresholds[f] <- quantile(expr_raw[, factor_idx[f]], probs = quantile_prob)
}}

# Test parameters
pairs_list <- {to_r_char_vec(ref_pairs)}
radii_list <- {to_r_num_vec(test_radii_actual)}
annular_width <- {ring_width}

# Output: pairs x radii x factors x genes
n_pairs <- length(pairs_list)
n_radii <- length(radii_list)
result_array <- array(0, dim = c(n_pairs, n_radii, n_factors, n_genes))

cat(sprintf("Computing I_ND for %d pairs x %d radii x %d factors x %d genes\\n",
    n_pairs, n_radii, n_factors, n_genes))

for (pi in seq_len(n_pairs)) {{
    pair_str <- pairs_list[pi]
    parts <- strsplit(pair_str, "->")[[1]]
    sender_type <- parts[1]
    receiver_type <- parts[2]

    sender_idx <- which(cell_types == sender_type)
    receiver_idx <- which(cell_types == receiver_type)

    sender_coords <- coords[sender_idx, , drop = FALSE]
    receiver_coords <- coords[receiver_idx, , drop = FALSE]

    # Float32-precision normalized expression
    sender_factor_norm <- expr_norm[sender_idx, factor_idx, drop = FALSE]
    receiver_expr_norm <- expr_norm[receiver_idx, , drop = FALSE]

    for (ri in seq_len(n_radii)) {{
        radius <- radii_list[ri]
        inner_radius <- max(0, radius - annular_width)

        t1 <- proc.time()

        # Build weight matrix in float32 (matching GPU precision)
        W <- create_gaussian_weights_cpp(
            sender_coords, receiver_coords,
            radius = radius,
            inner_radius = inner_radius,
            sigma = radius / 3.0,
            min_weight = 1e-6,
            use_float32 = TRUE
        )

        if (sum(W) == 0) {{
            cat(sprintf("  %s R=%d: no connections, skipping\\n", pair_str, radius))
            next
        }}

        # Compute I_ND using C++ float32 path (spatial lags + cosine sim)
        # Same function used by batch_compute_all_pairs_cpp
        ind_matrix <- compute_ind_matrix_cpp(
            sender_factor_norm,
            receiver_expr_norm,
            W,
            quant_thresholds,
            min_cells = 10L,
            use_float32 = TRUE
        )

        # Replace NA with 0
        ind_matrix[is.na(ind_matrix)] <- 0
        result_array[pi, ri, , ] <- ind_matrix

        elapsed <- (proc.time() - t1)[3]
        cat(sprintf("  %s R=%d (%.1fs)\\n", pair_str, as.integer(radius), elapsed))
    }}
}}

# Save to HDF5
cat(sprintf("Saving to %s\\n", "{r_output}"))
if (file.exists("{r_output}")) file.remove("{r_output}")
h5createFile("{r_output}")
h5write(result_array, "{r_output}", "morans_i_matrices")
h5write(pairs_list, "{r_output}", "pairs")
h5write(radii_list, "{r_output}", "radii")
h5write(factor_genes, "{r_output}", "factor_genes")
h5write(gene_names, "{r_output}", "target_genes")
H5close()
cat("Done.\\n")
'''

    r_script.write_text(r_code)

    # Run R script
    print("  Running R script...")
    t0 = time.time()
    env = os.environ.copy()
    result = subprocess.run(
        ['bash', '-c', f'module load R/4.3 2>/dev/null; Rscript {r_script}'],
        capture_output=True, text=True, timeout=3600, env=env
    )

    if result.returncode != 0:
        print(f"  R script failed (return code {result.returncode})")
        print(f"  STDOUT: {result.stdout[-2000:]}")
        print(f"  STDERR: {result.stderr[-2000:]}")
        return None

    print(f"  R script completed ({time.time()-t0:.1f}s)")
    print(f"  R output: {result.stdout[-500:]}")

    # Load R results and compare
    if not r_output.exists():
        print("  ERROR: R output file not produced")
        return None

    with h5py.File(r_output, 'r') as f:
        r_matrix_raw = f['morans_i_matrices'][:]
        r_pairs = [x.decode() if isinstance(x, bytes) else str(x) for x in f['pairs'][:]]
        r_radii = f['radii'][:].tolist()

    # R writes in column-major (Fortran) order; Python reads row-major
    # R array dim=(pairs, radii, factors, genes) → Python reads as (genes, factors, radii, pairs)
    # Transpose to get (pairs, radii, factors, genes)
    r_matrix = r_matrix_raw.transpose()
    print(f"  R output shape (after transpose): {r_matrix.shape}")

    # Compare to reference
    print("\n  Comparing R results to reference...")
    all_stats = []
    r_slices = {}
    for pair_str in ref_pairs:
        pi_ref = ref_pairs.index(pair_str)
        pi_r = r_pairs.index(pair_str)

        for ri_r, radius in enumerate(test_radii_actual):
            ri_ref = ref_radii.index(radius)

            test_mat = r_matrix[pi_r, ri_r, :, :]
            ref_mat = ref_data['matrix'][pi_ref, ri_ref, factor_slice, :]
            r_slices[(pair_str, radius)] = test_mat

            label = f"R-lowlevel {pair_str} R={radius}"
            stats = compare_matrices(test_mat, ref_mat, label, verbose=False)
            all_stats.append(stats)

    # Summary
    max_abs_all = max(s['max_abs_diff'] for s in all_stats)
    max_rel_all = max(s['max_rel_diff'] for s in all_stats)
    mean_close = np.mean([s['close_1e6_pct'] for s in all_stats])

    print(f"\n  SIGDISCOV R (LOW-LEVEL) SUMMARY:")
    print(f"  Slices tested:     {len(all_stats)}")
    print(f"  Max abs diff:      {max_abs_all:.2e}")
    print(f"  Max rel diff:      {max_rel_all:.2e}")
    print(f"  Mean within 1e-6:  {mean_close:.1f}%")

    worst = sorted(all_stats, key=lambda s: s['max_abs_diff'], reverse=True)[:5]
    print(f"\n  Worst 5 slices:")
    for s in worst:
        print(f"    {s['label']:50s}  max_abs={s['max_abs_diff']:.2e}  "
              f"within_1e-4={s['close_1e4_pct']:.1f}%")

    return {
        'package': 'sigdiscov_r_lowlevel',
        'stats': all_stats,
        'max_abs_diff': max_abs_all,
        'max_rel_diff': max_rel_all,
        'mean_close_1e6_pct': mean_close,
        'all_slices': r_slices,
    }


# ============================================================================
# TEST 3: sigdiscov R package (batch function, R's own algorithm)
# ============================================================================

def test_sigdiscov_r_batch(ref_data, quantile_val, test_radii, ring_width,
                           n_factors_limit=None, tmpdir=None):
    """Test sigdiscov R package using its batch function (may have algorithmic differences)."""
    print("\n" + "=" * 70)
    print("TEST: sigdiscov R package (batch function)")
    print("=" * 70)
    print("  Tests batch_compute_all_pairs_cpp which should match")
    print("  compute_morans.py algorithm (population std, fixed-width rings,")
    print("  all-cells quantile, normalized>raw comparison).")

    ref_pairs = ref_data['pairs']
    ref_radii = ref_data['radii']
    ref_factors = ref_data['factor_genes']
    ref_targets = ref_data['target_genes']

    if n_factors_limit and n_factors_limit < len(ref_factors):
        factor_subset = ref_factors[:n_factors_limit]
        factor_slice = slice(0, n_factors_limit)
    else:
        factor_subset = ref_factors
        factor_slice = slice(None)
        n_factors_limit = len(ref_factors)

    radii_indices = [ref_radii.index(r) for r in test_radii if r in ref_radii]
    test_radii_actual = [ref_radii[i] for i in radii_indices]

    if tmpdir is None:
        tmpdir = Path(tempfile.mkdtemp(prefix='r_batch_'))

    r_output = tmpdir / 'r_batch_results.h5'
    r_script = tmpdir / 'test_batch.R'

    pairs_json = json.dumps(ref_pairs)
    # Convert Python lists to R vector syntax
    def to_r_char_vec(lst):
        return 'c(' + ', '.join(f'"{x}"' for x in lst) + ')'
    def to_r_num_vec(lst):
        return 'c(' + ', '.join(str(x) for x in lst) + ')'

    # Export data from Python
    data_path = _export_data_for_r(tmpdir, ref_data, factor_subset, n_factors_limit)

    r_code = f'''
library(sigdiscov)
library(rhdf5)

cat("Loading exported data...\\n")
t0 <- proc.time()
expr_filtered <- h5read("{data_path}", "expr")  # cells x genes
coords <- h5read("{data_path}", "coords")  # cells x 2
cell_types <- as.character(h5read("{data_path}", "cell_types"))
gene_names <- as.character(h5read("{data_path}", "gene_names"))
H5close()

n_cells <- nrow(expr_filtered)
n_genes <- ncol(expr_filtered)
cat(sprintf("Loaded: %d cells, %d genes (%.1fs)\\n",
    n_cells, n_genes, (proc.time() - t0)[3]))

# Build SCData object (genes x cells format for batch_compute_all_pairs_cpp)
data_obj <- list(
    expr = t(expr_filtered),  # genes x cells
    coords = coords,
    cell_types = cell_types,
    gene_names = gene_names,
    n_genes = n_genes,
    n_cells = n_cells
)

# Factor genes
factor_genes <- {to_r_char_vec(factor_subset)}
factor_idx_r <- match(factor_genes, gene_names)
valid <- !is.na(factor_idx_r)
factor_idx_r <- factor_idx_r[valid]
factor_genes <- factor_genes[valid]
n_factors <- length(factor_genes)

cat(sprintf("Factors: %d\\n", n_factors))

# Build pairs dataframe
all_pairs <- {to_r_char_vec(ref_pairs)}
pairs_df <- data.frame(
    sender = sapply(strsplit(all_pairs, "->"), "[", 1),
    receiver = sapply(strsplit(all_pairs, "->"), "[", 2),
    stringsAsFactors = FALSE
)

radii <- {to_r_num_vec(test_radii_actual)}

cat(sprintf("Running batch computation: %d pairs x %d radii...\\n",
    nrow(pairs_df), length(radii)))
t0 <- proc.time()

# Call C++ batch function directly
result_list <- batch_compute_all_pairs_cpp(
    expr_matrix = data_obj$expr,
    coords = coords,
    cell_types = cell_types,
    factor_indices = factor_idx_r - 1L,  # 0-based
    radii = radii,
    pairs = pairs_df,
    quantile_prob = {quantile_val},
    min_cells = 10L,
    min_connections = 10L,
    ring_width = {ring_width},
    use_float32 = TRUE,
    verbose = TRUE
)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("Batch computation: %.1fs\\n", elapsed))

# Convert to 4D array
n_pairs <- nrow(pairs_df)
n_radii <- length(radii)
result_array <- array(NA_real_, dim = c(n_pairs, n_radii, n_factors, n_genes))
for (r in seq_len(n_radii)) {{
    result_array[, r, , ] <- result_list[[r]]
}}

# Replace NA with 0 for comparison
result_array[is.na(result_array)] <- 0

# Save
cat(sprintf("Saving to %s\\n", "{r_output}"))
if (file.exists("{r_output}")) file.remove("{r_output}")
h5createFile("{r_output}")
h5write(result_array, "{r_output}", "morans_i_matrices")
h5write(all_pairs, "{r_output}", "pairs")
h5write(radii, "{r_output}", "radii")
h5write(factor_genes, "{r_output}", "factor_genes")
h5write(gene_names, "{r_output}", "target_genes")
H5close()
cat("Done.\\n")
'''

    r_script.write_text(r_code)

    print("  Running R batch script...")
    t0 = time.time()
    result = subprocess.run(
        ['bash', '-c', f'module load R/4.3 2>/dev/null; Rscript {r_script}'],
        capture_output=True, text=True, timeout=3600,
    )

    if result.returncode != 0:
        print(f"  R script failed (return code {result.returncode})")
        print(f"  STDOUT: {result.stdout[-2000:]}")
        print(f"  STDERR: {result.stderr[-2000:]}")
        return None

    print(f"  R batch completed ({time.time()-t0:.1f}s)")
    print(f"  R output: {result.stdout[-500:]}")

    if not r_output.exists():
        print("  ERROR: R output file not produced")
        return None

    with h5py.File(r_output, 'r') as f:
        r_matrix_raw = f['morans_i_matrices'][:]
        r_pairs = [x.decode() if isinstance(x, bytes) else str(x) for x in f['pairs'][:]]

    r_matrix = r_matrix_raw.transpose()
    print(f"  R output shape (after transpose): {r_matrix.shape}")

    # Compare
    print("\n  Comparing R batch results to reference...")
    all_stats = []
    r_slices = {}
    for pair_str in ref_pairs:
        pi_ref = ref_pairs.index(pair_str)
        pi_r = r_pairs.index(pair_str)

        for ri_r, radius in enumerate(test_radii_actual):
            ri_ref = ref_radii.index(radius)

            test_mat = r_matrix[pi_r, ri_r, :, :]
            ref_mat = ref_data['matrix'][pi_ref, ri_ref, factor_slice, :]
            r_slices[(pair_str, radius)] = test_mat

            label = f"R-batch {pair_str} R={radius}"
            stats = compare_matrices(test_mat, ref_mat, label, verbose=False)
            all_stats.append(stats)

    max_abs_all = max(s['max_abs_diff'] for s in all_stats)
    max_rel_all = max(s['max_rel_diff'] for s in all_stats)
    mean_close = np.mean([s['close_1e6_pct'] for s in all_stats])

    print(f"\n  SIGDISCOV R (BATCH) SUMMARY:")
    print(f"  Slices tested:     {len(all_stats)}")
    print(f"  Max abs diff:      {max_abs_all:.2e}")
    print(f"  Max rel diff:      {max_rel_all:.2e}")
    print(f"  Mean within 1e-6:  {mean_close:.1f}%")

    worst = sorted(all_stats, key=lambda s: s['max_abs_diff'], reverse=True)[:5]
    print(f"\n  Worst 5 slices:")
    for s in worst:
        print(f"    {s['label']:50s}  max_abs={s['max_abs_diff']:.2e}  "
              f"within_1e-4={s['close_1e4_pct']:.1f}%")

    return {
        'package': 'sigdiscov_r_batch',
        'stats': all_stats,
        'max_abs_diff': max_abs_all,
        'max_rel_diff': max_rel_all,
        'mean_close_1e6_pct': mean_close,
        'all_slices': r_slices,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='Cross-package equivalence test')
    parser.add_argument('--gpu', type=str, default='0', help='GPU device ID')
    parser.add_argument('--n-factors', type=int, default=20,
                        help='Number of factors to test (default: 20, use 0 for all)')
    parser.add_argument('--n-radii', type=int, default=3,
                        help='Number of radii to test (default: 3)')
    parser.add_argument('--quantile', type=str, default='q0',
                        choices=['q0', 'q25', 'q50'],
                        help='Quantile level to test (default: q0)')
    parser.add_argument('--ring-width', type=float, default=20.0,
                        help='Annular ring width in um (default: 20.0)')
    parser.add_argument('--keep-output', action='store_true')
    parser.add_argument('--skip-python', action='store_true',
                        help='Skip sigdiscovpy test')
    parser.add_argument('--skip-r-lowlevel', action='store_true',
                        help='Skip R low-level test')
    parser.add_argument('--skip-r-batch', action='store_true',
                        help='Skip R batch test')
    args = parser.parse_args()

    # Validate prerequisites
    assert DATA_PATH.exists(), f"Data file not found: {DATA_PATH}"
    assert FACTOR_PATH.exists(), f"Factor file not found: {FACTOR_PATH}"

    quantile_map = {'q0': 0.0, 'q25': 0.25, 'q50': 0.50}
    quantile_val = quantile_map[args.quantile]

    # Load reference
    print(f"Loading reference data ({args.quantile})...")
    ref_data = load_reference(args.quantile)
    print(f"  Shape: {ref_data['matrix'].shape}")

    # Select radii subset
    if args.n_radii >= len(ALL_RADII):
        test_radii = ALL_RADII
    else:
        # Pick evenly spaced radii including first and last
        indices = np.linspace(0, len(ALL_RADII) - 1, args.n_radii, dtype=int)
        test_radii = [ALL_RADII[i] for i in indices]

    n_factors = args.n_factors if args.n_factors > 0 else None
    ring_width = args.ring_width

    print("=" * 70)
    print("CROSS-PACKAGE EQUIVALENCE TEST")
    print("=" * 70)
    print(f"Reference:  ref_{args.quantile}.h5")
    print(f"Quantile:   {quantile_val}")
    print(f"Radii:      {test_radii}")
    print(f"Ring width: {ring_width}")
    print(f"Factors:    {n_factors or 'all'}")
    print(f"Pairs:      {len(CELL_PAIRS)}")

    tmpdir = Path(tempfile.mkdtemp(prefix='xpkg_test_'))
    results = {}

    try:
        # Test 1: sigdiscovpy
        if not args.skip_python:
            r = test_sigdiscovpy(ref_data, quantile_val, test_radii, ring_width, n_factors)
            if r:
                results['sigdiscovpy'] = r

        # Test 2: R low-level (matching Python algorithm)
        if not args.skip_r_lowlevel:
            r = test_sigdiscov_r(ref_data, quantile_val, test_radii, ring_width,
                                n_factors, tmpdir)
            if r:
                results['sigdiscov_r_lowlevel'] = r

        # Test 3: R batch (R's own algorithm)
        if not args.skip_r_batch:
            r = test_sigdiscov_r_batch(ref_data, quantile_val, test_radii, ring_width,
                                      n_factors, tmpdir)
            if r:
                results['sigdiscov_r_batch'] = r

    finally:
        if args.keep_output:
            print(f"\nOutput kept at: {tmpdir}")
        else:
            shutil.rmtree(tmpdir, ignore_errors=True)

    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"{'Package':<30s} {'Max abs diff':>12s} {'Max rel diff':>12s} "
          f"{'Within 1e-6':>12s} {'Verdict':>10s}")
    print("-" * 80)

    for name, r in results.items():
        # IDENTICAL: ~1e-5 (float32 GPU operation ordering noise)
        # CLOSE: ~1e-3 (float64 vs float32 precision, same algorithm)
        # DIFFERENT: >1e-2 (algorithmic differences)
        verdict = "IDENTICAL" if r['max_abs_diff'] < 1e-4 else \
                  "CLOSE" if r['max_abs_diff'] < 1e-2 else "DIFFERENT"
        print(f"  {name:<28s} {r['max_abs_diff']:>12.2e} {r['max_rel_diff']:>12.2e} "
              f"{r['mean_close_1e6_pct']:>11.1f}% {verdict:>10s}")

    # Cross-comparison: R vs sigdiscovpy (bypasses old reference)
    if 'sigdiscovpy' in results and len(results) > 1:
        py_slices = results['sigdiscovpy']['all_slices']
        print("\n" + "=" * 70)
        print("CROSS-COMPARISON (R vs sigdiscovpy, same-run)")
        print("=" * 70)
        for name, r in results.items():
            if name == 'sigdiscovpy':
                continue
            max_diff = 0
            n_compared = 0
            for key, r_slice in r['all_slices'].items():
                if key in py_slices:
                    d = np.abs(r_slice - py_slices[key])
                    max_diff = max(max_diff, np.nanmax(d))
                    n_compared += 1
            verdict = "IDENTICAL" if max_diff < 1e-4 else \
                      "CLOSE" if max_diff < 1e-2 else "DIFFERENT"
            print(f"  {name:<28s} vs sigdiscovpy: max_abs={max_diff:.2e}  "
                  f"({n_compared} slices)  {verdict}")

    return 0 if all(r['max_abs_diff'] < 1e-2 for r in results.values()) else 1


if __name__ == '__main__':
    sys.exit(main())

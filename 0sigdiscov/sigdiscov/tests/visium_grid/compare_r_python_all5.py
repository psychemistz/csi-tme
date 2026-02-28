#!/usr/bin/env python
"""
Compare R sigdiscov and Python sigdiscovpy grid-weight pairwise Moran's I
for all 5 Visium datasets against C/MKL reference.
"""

import sys
import os
import subprocess
import time
import numpy as np
import pandas as pd

sys.path.insert(0, '/vf/users/parks34/projects/0sigdiscov/pkg_dev/sigdiscovpy')

from sigdiscovpy import create_grid_weights, pairwise_moran
from sigdiscovpy.core.weights import get_weight_sum
from sigdiscovpy.io.loaders import parse_spot_names

DATA_DIR = '/data/Jiang_Lab/Data/Beibei/paired_wise_moranI'
N_GENES = 100  # subset for speed

# R environment setup
R_ENV = os.environ.copy()
R_MODULES = "module load gcc/11.3.0 R/4.3 2>&1 | tail -1"
RSCRIPT = "/usr/local/apps/R/4.3/4.3.0/bin/Rscript"
GCC_LIB = subprocess.run(
    "module load gcc/11.3.0 2>&1 && gcc -print-file-name=libstdc++.so",
    shell=True, capture_output=True, text=True
).stdout.strip()
GCC_LIB_DIR = os.path.dirname(GCC_LIB) if GCC_LIB else ""


def run_r(vst_file, n_genes, output_file):
    """Run R script to compute pairwise Moran's I with grid weights."""
    cmd = (
        f"source /etc/profile.d/modules.sh 2>/dev/null; "
        f"module load gcc/11.3.0 R/4.3 2>/dev/null; "
        f"export LD_LIBRARY_PATH=$(dirname $(gcc -print-file-name=libstdc++.so)):$LD_LIBRARY_PATH; "
        f"{RSCRIPT} /tmp/run_r_grid.R {vst_file} {n_genes} {output_file}"
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        print(f"  R STDERR: {result.stderr[-500:]}")
        raise RuntimeError(f"R script failed: exit code {result.returncode}")
    print(result.stdout.strip())
    return output_file


def run_python(vst_file, n_genes):
    """Run Python to compute pairwise Moran's I with grid weights."""
    vst = pd.read_csv(vst_file, sep='\t', index_col=0)
    spot_names = list(vst.columns)
    spot_coords = parse_spot_names(spot_names)

    W = create_grid_weights(spot_coords, max_radius=3, same_spot=False, platform="visium")
    weight_sum = get_weight_sum(W)

    n = min(n_genes, vst.shape[0])
    expr = vst.values[:n, :].astype(np.float64)

    I_matrix = pairwise_moran(expr, coords=spot_coords, W=W, use_gpu=False)
    return I_matrix, weight_sum, W.nnz


def load_reference(ref_file, n_genes):
    """Load lower-triangular C/MKL reference and build symmetric matrix."""
    rows = []
    with open(ref_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= n_genes:
                break
            vals = [float(v) for v in line.strip().split('\t')]
            rows.append(vals)

    mat = np.zeros((n_genes, n_genes), dtype=np.float64)
    for i, row_vals in enumerate(rows):
        for j, v in enumerate(row_vals):
            mat[i, j] = v
            mat[j, i] = v
    return mat


def round_to_sigfigs(x, n=6):
    """Round array to n significant figures (matching C printf %g default)."""
    with np.errstate(divide='ignore', invalid='ignore'):
        magnitude = np.floor(np.log10(np.abs(x)))
        magnitude = np.where(np.isfinite(magnitude), magnitude, 0)
        factor = 10.0 ** (n - 1 - magnitude)
        return np.round(x * factor) / factor


def compare(name, a, b, sigfigs=None):
    """Compare two matrices and print stats."""
    if sigfigs:
        a = round_to_sigfigs(a, sigfigs)
        b = round_to_sigfigs(b, sigfigs)

    diff = np.abs(a - b)
    max_d = diff.max()
    mean_d = diff.mean()

    nontrivial = np.abs(b) > 1e-4
    max_rel = (diff[nontrivial] / np.abs(b[nontrivial])).max() if nontrivial.any() else 0.0

    pct_exact = 100 * (diff == 0).sum() / diff.size

    sfx = f" ({sigfigs} sig figs)" if sigfigs else ""
    print(f"  {name:30s}  max_abs={max_d:.2e}  mean_abs={mean_d:.2e}  "
          f"max_rel={max_rel:.2e}  exact={pct_exact:.1f}%{sfx}")
    return max_d


print("=" * 90)
print(f"Cross-package comparison: Python vs R vs C/MKL reference (first {N_GENES} genes)")
print("=" * 90)

all_results = []

for ds in range(1, 6):
    vst_file = f"{DATA_DIR}/{ds}_vst.tsv"
    ref_file = f"{DATA_DIR}/{ds}_spatial_sig_vst_inhouse_s0_r3.tsv"
    r_out = f"/tmp/r_grid_ds{ds}.bin"

    print(f"\n{'─' * 90}")
    print(f"Dataset {ds}: {vst_file}")
    print(f"{'─' * 90}")

    # Python
    print("\n[Python]")
    t0 = time.time()
    py_mat, py_ws, py_nnz = run_python(vst_file, N_GENES)
    py_time = time.time() - t0
    print(f"  Weight sum: {py_ws:.6f}, NNZ: {py_nnz}, Time: {py_time:.1f}s")

    # R
    print("\n[R]")
    t0 = time.time()
    run_r(vst_file, N_GENES, r_out)
    r_time = time.time() - t0

    n = min(N_GENES, py_mat.shape[0])
    r_mat = np.fromfile(r_out, dtype=np.float64).reshape(n, n)
    r_ws = np.fromfile(r_out + ".ws", dtype=np.float64)[0]
    print(f"  R weight sum: {r_ws:.6f}, Time: {r_time:.1f}s")

    # Reference
    ref_mat = load_reference(ref_file, n)

    # Compare
    print(f"\n  Weight sum comparison: Python={py_ws:.6f}  R={r_ws:.6f}  diff={abs(py_ws - r_ws):.2e}")

    print("\n  Full precision:")
    d_py_ref = compare("Python vs C/MKL ref", py_mat, ref_mat)
    d_r_ref = compare("R vs C/MKL ref", r_mat, ref_mat)
    d_py_r = compare("Python vs R", py_mat, r_mat)

    print("\n  6 significant figures:")
    d_py_ref_6 = compare("Python vs C/MKL ref", py_mat, ref_mat, sigfigs=6)
    d_r_ref_6 = compare("R vs C/MKL ref", r_mat, ref_mat, sigfigs=6)
    d_py_r_6 = compare("Python vs R", py_mat, r_mat, sigfigs=6)

    all_results.append({
        'dataset': ds,
        'py_vs_ref': d_py_ref,
        'r_vs_ref': d_r_ref,
        'py_vs_r': d_py_r,
        'py_vs_ref_6': d_py_ref_6,
        'r_vs_ref_6': d_r_ref_6,
        'py_vs_r_6': d_py_r_6,
        'ws_diff': abs(py_ws - r_ws),
    })

# Summary
print(f"\n{'=' * 90}")
print("SUMMARY")
print(f"{'=' * 90}")
print(f"\n  Full precision:")
print(f"  {'DS':>3s}  {'Py vs Ref':>12s}  {'R vs Ref':>12s}  {'Py vs R':>12s}")
for r in all_results:
    print(f"  {r['dataset']:3d}  {r['py_vs_ref']:12.2e}  {r['r_vs_ref']:12.2e}  {r['py_vs_r']:12.2e}")

print(f"\n  6 significant figures:")
print(f"  {'DS':>3s}  {'Py vs Ref':>12s}  {'R vs Ref':>12s}  {'Py vs R':>12s}")
for r in all_results:
    print(f"  {r['dataset']:3d}  {r['py_vs_ref_6']:12.2e}  {r['r_vs_ref_6']:12.2e}  {r['py_vs_r_6']:12.2e}")

max_py_r = max(r['py_vs_r'] for r in all_results)
max_ref_6 = max(r['py_vs_ref_6'] for r in all_results)
print(f"\nPython vs R (full precision): max diff = {max_py_r:.2e}")
print(f"Python vs Ref (6 sig figs):  max diff = {max_ref_6:.2e}")
if max_ref_6 == 0:
    print("PASS: All three implementations identical at 6 significant figures")

#!/usr/bin/env python
"""Check if the 1e-7 difference is from TSV decimal truncation."""
import numpy as np

REF = '/data/Jiang_Lab/Data/Beibei/paired_wise_moranI/1_spatial_sig_vst_inhouse_s0_r3.tsv'

# Count significant digits in reference values
max_digits = 0
digit_counts = []
with open(REF, 'r') as f:
    for i, line in enumerate(f):
        if i >= 10:
            break
        for val_str in line.strip().split('\t'):
            # Strip sign and leading zeros
            s = val_str.lstrip('-')
            # Remove leading "0."
            if '.' in s:
                integer_part, decimal_part = s.split('.')
                if integer_part == '0':
                    # Count significant digits in decimal part (strip leading zeros)
                    sig = decimal_part.lstrip('0')
                    n_digits = len(sig)
                else:
                    n_digits = len(integer_part) + len(decimal_part)
            else:
                n_digits = len(s.lstrip('0'))
            digit_counts.append(n_digits)

print(f"Significant digits in reference values (first 10 rows):")
print(f"  Min: {min(digit_counts)}")
print(f"  Max: {max(digit_counts)}")
print(f"  Mean: {np.mean(digit_counts):.1f}")
print(f"  Median: {np.median(digit_counts):.1f}")

# Distribution
from collections import Counter
c = Counter(digit_counts)
for k in sorted(c):
    print(f"  {k} digits: {c[k]} values ({100*c[k]/len(digit_counts):.1f}%)")

# Quantify truncation error
# If values have 6 significant digits, relative precision is ~1e-6
# For a value like 0.00258679 (6 sig digits), precision is ~1e-8
# But the max_abs_diff we see is 4.58e-7
print(f"\nExpected precision for 6 significant digits:")
print(f"  Value ~0.01: precision ~0.01 * 1e-6 = 1e-8")
print(f"  Value ~0.05: precision ~0.05 * 1e-6 = 5e-8")
print(f"  Value ~0.5:  precision ~0.5 * 1e-6  = 5e-7")
print(f"  Value ~1.0:  precision ~1.0 * 1e-6  = 1e-6")
print(f"\nObserved max_abs_diff: 4.58e-07")

# Check the largest values in first 50 rows to see what magnitude causes the max diff
import sys
sys.path.insert(0, '/vf/users/parks34/projects/0sigdiscov/pkg_dev/sigdiscovpy')
from sigdiscovpy import create_grid_weights
from sigdiscovpy.core.weights import get_weight_sum
from sigdiscovpy.core.normalization import standardize_matrix
from sigdiscovpy.io.loaders import parse_spot_names
import pandas as pd

vst = pd.read_csv('/data/Jiang_Lab/Data/Beibei/paired_wise_moranI/1_vst.tsv', sep='\t', index_col=0)
spot_coords = parse_spot_names(list(vst.columns))
W = create_grid_weights(spot_coords, max_radius=3, same_spot=False, platform="visium")
expr = vst.values[:50, :].astype(np.float64)
Z = standardize_matrix(expr, axis=1, use_gpu=False)
ZW = Z @ W
I_py = (ZW @ Z.T) / get_weight_sum(W)

# Load ref
ref_rows = []
with open(REF, 'r') as f:
    for i, line in enumerate(f):
        if i >= 50:
            break
        ref_rows.append([float(v) for v in line.strip().split('\t')])
ref_mat = np.zeros((50, 50), dtype=np.float64)
for i, rv in enumerate(ref_rows):
    for j, v in enumerate(rv):
        ref_mat[i, j] = v
        ref_mat[j, i] = v

diff = np.abs(I_py - ref_mat)
# Find where max diff occurs
i_max, j_max = np.unravel_index(diff.argmax(), diff.shape)
print(f"\nMax diff location: ({i_max}, {j_max})")
print(f"  Python value:    {I_py[i_max, j_max]:.15f}")
print(f"  Reference value: {ref_mat[i_max, j_max]:.15f}")
print(f"  Abs diff:        {diff[i_max, j_max]:.2e}")
print(f"  Ref sig digits:  ", end="")

# Get the raw string from the ref file
with open(REF, 'r') as f:
    for i, line in enumerate(f):
        if i == i_max:
            vals = line.strip().split('\t')
            if j_max < len(vals):
                print(f"'{vals[j_max]}' ({len(vals[j_max].replace('-','').replace('.','').lstrip('0'))} sig digits)")
            elif i_max < len(vals):
                # Lower triangle: use (j_max, i_max)
                pass
            break

# Show relative error vs magnitude
print(f"\nRelative error vs value magnitude:")
nontrivial = np.abs(ref_mat) > 1e-4
magnitudes = np.abs(ref_mat[nontrivial])
rel_errors = diff[nontrivial] / magnitudes
for pct in [50, 90, 95, 99, 100]:
    p = np.percentile(rel_errors, pct) if pct < 100 else rel_errors.max()
    print(f"  {pct}th percentile relative error: {p:.2e}")

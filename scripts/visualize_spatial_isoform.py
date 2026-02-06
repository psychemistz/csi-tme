#!/usr/bin/env python3
"""
Generic spatial isoform visualization for any gene/sample.
5-panel figure:
  A) Gene expressing spots
  B) Tumor (Malignant) proportion from SpaCET
  C) Macrophage proportion from SpaCET
  D) Isoform 1 fraction
  E) Cell-type vs isoform 1 correlation (bar chart, major lineage types)
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import anndata as ad
import scipy.sparse as sp
from scipy.stats import spearmanr
from pathlib import Path

PROJECT_ROOT = Path("/vf/users/parks34/projects/3csitme")
SR_DIR = PROJECT_ROOT / "data/zenodo_16905935/human_glioma_sr"
REPORTS_DIR = PROJECT_ROOT / "reports"
PROP_DIR = REPORTS_DIR / "spacet_proportions"

# Major lineage cell types from SpaCET (Level 1)
MAJOR_LINEAGE = [
    "Malignant", "CAF", "Endothelial", "Plasma", "B cell",
    "T CD4", "T CD8", "NK", "cDC", "pDC",
    "Macrophage", "Mast", "Neutrophil"
]


def load_spatial_and_isoforms(sample_id, gene):
    """Load spatial coords and isoform data from h5ad."""
    iso_h5ad = SR_DIR / sample_id / "iso.quant.h5ad"
    adata = ad.read_h5ad(iso_h5ad)

    # Spatial coordinates
    if 'array_row' in adata.obs.columns and 'array_col' in adata.obs.columns:
        coords = np.column_stack([
            adata.obs['array_col'].values.astype(float),
            adata.obs['array_row'].values.astype(float)
        ])
    elif 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    else:
        raise ValueError("No spatial coordinates found")

    obs_names = adata.obs_names.tolist()

    # Find gene isoforms
    gene_col = None
    for col in ['gene_symbol', 'gene', 'gene_name']:
        if col in adata.var.columns:
            gene_col = col
            break

    if gene_col:
        mask = adata.var[gene_col] == gene
        iso_names = adata.var_names[mask].tolist()
    else:
        iso_names = [v for v in adata.var_names if v.startswith(gene)]

    if len(iso_names) < 2:
        raise ValueError(f"Only {len(iso_names)} isoforms found for {gene}")

    # Get isoform counts
    iso_idx = [adata.var_names.tolist().index(n) for n in iso_names]
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    iso_counts = X[:, iso_idx]

    # Isoform 1 fraction
    total = iso_counts.sum(axis=1)
    iso1_frac = np.where(total > 0, iso_counts[:, 0] / total, np.nan)

    return coords, obs_names, iso1_frac, total, iso_names


def load_spacet_proportions(sample_id):
    """Load SpaCET per-spot cell type proportions."""
    prop_file = PROP_DIR / f"{sample_id}_proportions.csv"
    if not prop_file.exists():
        raise FileNotFoundError(f"SpaCET proportions not found: {prop_file}")
    df = pd.read_csv(prop_file)
    df = df.set_index('spot')
    return df


def create_figure(gene, sample_id, output_path=None):
    """Create 5-panel spatial isoform visualization."""

    print(f"Loading {gene} data for {sample_id}...")
    coords, obs_names, iso1_frac, gene_total, iso_names = \
        load_spatial_and_isoforms(sample_id, gene)
    print(f"  {len(obs_names)} spots, {len(iso_names)} {gene} isoforms")
    for i, name in enumerate(iso_names):
        print(f"    iso{i+1}: {name}")

    props = load_spacet_proportions(sample_id)
    print(f"  SpaCET: {props.shape[0]} spots x {props.shape[1]} cell types")

    # Align spots
    common = [s for s in obs_names if s in props.index]
    spot_idx = [obs_names.index(s) for s in common]
    props = props.loc[common]

    coords_c = coords[spot_idx]
    iso1_c = iso1_frac[spot_idx]
    total_c = gene_total[spot_idx]

    # Expressed mask
    expressed = ~np.isnan(iso1_c) & (total_c > 0)
    n_expr = expressed.sum()
    print(f"  {gene}-expressing spots: {n_expr}/{len(common)}")

    # Use only major lineage cell types present in data
    available_major = [ct for ct in MAJOR_LINEAGE if ct in props.columns]

    # Compute Spearman correlations for expressed spots
    cor_data = []
    for ct in available_major:
        ct_vals = props[ct].values[expressed]
        iso_vals = iso1_c[expressed]
        valid = ~np.isnan(iso_vals)
        if valid.sum() > 10:
            r, p = spearmanr(ct_vals[valid], iso_vals[valid])
            cor_data.append({'cell_type': ct, 'r': r, 'p': p})
    cor_df = pd.DataFrame(cor_data).sort_values('r', ascending=True)

    # ---- Figure: 2 rows ----
    fig = plt.figure(figsize=(20, 10))
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1, 1],
                  hspace=0.3, wspace=0.3)

    pt = 10
    alpha_bg = 0.3
    alpha_fg = 0.85

    # Row 1: Spatial maps (A-D)

    # A) Gene expressing spots
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.scatter(coords_c[:, 0], coords_c[:, 1], c='#e0e0e0', s=pt,
                 alpha=alpha_bg, linewidths=0)
    ax_a.scatter(coords_c[expressed, 0], coords_c[expressed, 1],
                 c='#2171b5', s=pt, alpha=alpha_fg, linewidths=0)
    ax_a.set_title(f'A) {gene} Expression\n({n_expr}/{len(common)} spots)',
                   fontsize=12, fontweight='bold')
    ax_a.invert_yaxis()
    ax_a.set_aspect('equal')
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    for spine in ax_a.spines.values():
        spine.set_visible(True)
        spine.set_color('black')

    # B) Tumor (Malignant) proportion
    ax_b = fig.add_subplot(gs[0, 1])
    tumor_vals = props['Malignant'].values if 'Malignant' in props.columns \
        else np.zeros(len(common))
    sc_b = ax_b.scatter(coords_c[:, 0], coords_c[:, 1],
                        c=tumor_vals, cmap='Reds', s=pt, alpha=alpha_fg,
                        vmin=0, vmax=max(tumor_vals.max(), 0.01), linewidths=0)
    ax_b.set_title('B) Tumor Proportion\n(SpaCET)', fontsize=12, fontweight='bold')
    ax_b.invert_yaxis()
    ax_b.set_aspect('equal')
    ax_b.set_xticks([])
    ax_b.set_yticks([])
    for spine in ax_b.spines.values():
        spine.set_visible(True)
    plt.colorbar(sc_b, ax=ax_b, fraction=0.046, pad=0.04, shrink=0.8)

    # C) Macrophage proportion
    ax_c = fig.add_subplot(gs[0, 2])
    macro_vals = props['Macrophage'].values if 'Macrophage' in props.columns \
        else np.zeros(len(common))
    sc_c = ax_c.scatter(coords_c[:, 0], coords_c[:, 1],
                        c=macro_vals, cmap='Purples', s=pt, alpha=alpha_fg,
                        vmin=0, vmax=max(macro_vals.max(), 0.01), linewidths=0)
    ax_c.set_title('C) Macrophage Proportion\n(SpaCET)', fontsize=12, fontweight='bold')
    ax_c.invert_yaxis()
    ax_c.set_aspect('equal')
    ax_c.set_xticks([])
    ax_c.set_yticks([])
    for spine in ax_c.spines.values():
        spine.set_visible(True)
    plt.colorbar(sc_c, ax=ax_c, fraction=0.046, pad=0.04, shrink=0.8)

    # D) Isoform 1 fraction
    ax_d = fig.add_subplot(gs[0, 3])
    ax_d.scatter(coords_c[~expressed, 0], coords_c[~expressed, 1],
                 c='#e0e0e0', s=pt, alpha=alpha_bg, linewidths=0)
    sc_d = ax_d.scatter(coords_c[expressed, 0], coords_c[expressed, 1],
                        c=iso1_c[expressed], cmap='RdYlBu_r', s=pt,
                        alpha=alpha_fg, vmin=0, vmax=1, linewidths=0)
    ax_d.set_title(f'D) {gene} Isoform 1 Fraction', fontsize=12, fontweight='bold')
    ax_d.invert_yaxis()
    ax_d.set_aspect('equal')
    ax_d.set_xticks([])
    ax_d.set_yticks([])
    for spine in ax_d.spines.values():
        spine.set_visible(True)
    plt.colorbar(sc_d, ax=ax_d, fraction=0.046, pad=0.04, shrink=0.8,
                 label='Isoform 1 Fraction')

    # Row 2: E) Correlation bar chart spanning full width
    ax_e = fig.add_subplot(gs[1, :])

    n_ct = len(cor_df)
    colors = ['#e74c3c' if r > 0 else '#3498db' for r in cor_df['r']]
    x_pos = np.arange(n_ct)

    bars = ax_e.bar(x_pos, cor_df['r'].values, color=colors,
                    edgecolor='black', linewidth=0.8, width=0.7)
    ax_e.set_xticks(x_pos)
    ax_e.set_xticklabels(cor_df['cell_type'].values, rotation=45, ha='right',
                         fontsize=11)
    ax_e.axhline(0, color='black', linewidth=0.5)
    ax_e.set_ylabel(f'Spearman r (vs {gene} Isoform 1 Fraction)', fontsize=12)
    ax_e.set_title(f'E) Cell Type Proportion vs {gene} Isoform 1 Correlation',
                   fontsize=12, fontweight='bold')
    ax_e.set_xlim(-0.8, n_ct - 0.2)

    # Auto-scale y axis
    r_min = cor_df['r'].min()
    r_max = cor_df['r'].max()
    y_pad = max(abs(r_min), abs(r_max)) * 0.3
    ax_e.set_ylim(r_min - y_pad, r_max + y_pad)

    # Add significance markers
    for i, (_, row) in enumerate(cor_df.iterrows()):
        sig = '***' if row['p'] < 0.001 else (
            '**' if row['p'] < 0.01 else ('*' if row['p'] < 0.05 else ''))
        if sig:
            y_pos = row['r'] + (0.015 if row['r'] >= 0 else -0.015)
            ax_e.text(i, y_pos, sig, va='bottom' if row['r'] >= 0 else 'top',
                      ha='center', fontsize=11, fontweight='bold')

    # Add r values inside bars
    for i, (_, row) in enumerate(cor_df.iterrows()):
        if abs(row['r']) > 0.04:
            ax_e.text(i, row['r'] / 2, f"{row['r']:.2f}",
                      ha='center', va='center', fontsize=8,
                      fontweight='bold', color='white')

    ax_e.spines['top'].set_visible(False)
    ax_e.spines['right'].set_visible(False)

    fig.suptitle(f'{sample_id} — {gene} Spatial Isoform Switching in Glioma TME',
                 fontsize=15, fontweight='bold', y=1.01)

    if output_path is None:
        output_path = REPORTS_DIR / f"{gene.lower()}_validation_{sample_id}.png"

    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", default="VCAN")
    parser.add_argument("--sample", default="ZH1019T1")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    create_figure(args.gene, args.sample, args.output)

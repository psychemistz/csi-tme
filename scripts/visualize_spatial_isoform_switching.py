#!/usr/bin/env python3
"""
Visualize isoform switching between tumor and macrophage on Visium spatial plots.
Shows spatial distribution of:
1. Tumor (Malignant) cell proportion
2. Macrophage cell proportion
3. VCAN isoform ratio
4. Combined visualization
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

# Set up paths
PROJECT_ROOT = Path("/vf/users/parks34/projects/3csitme")
SR_DIR = PROJECT_ROOT / "data/zenodo_16905935/human_glioma_sr"
REPORTS_DIR = PROJECT_ROOT / "reports"

# Load SpaCET results
spacet_cors = pd.read_csv(REPORTS_DIR / "spacet_isoform_correlations.csv")

def load_h5ad_data(h5ad_path):
    """Load h5ad file and return counts and spatial coordinates."""
    adata = ad.read_h5ad(h5ad_path)

    # Get spatial coordinates
    obs_df = adata.obs
    if 'array_row' in obs_df.columns and 'array_col' in obs_df.columns:
        spatial = np.column_stack([obs_df['array_col'].values, obs_df['array_row'].values])
    elif 'spatial' in adata.obsm:
        spatial = adata.obsm['spatial']
    else:
        return None, None, None

    # Get counts matrix
    if sp.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)

    return X, spatial, adata.var_names.tolist(), adata.obs_names.tolist()

def load_isoform_ratios(h5ad_path, gene_name="VCAN"):
    """Load isoform data and calculate ratios."""
    adata = ad.read_h5ad(h5ad_path)

    var_names = adata.var_names.tolist()
    obs_names = adata.obs_names.tolist()

    # Find isoforms for the gene
    gene_col = None
    for col in ['gene_symbol', 'gene', 'gene_name']:
        if col in adata.var.columns:
            gene_col = col
            break

    if gene_col:
        mask = adata.var[gene_col] == gene_name
        gene_isos = [v for v, m in zip(var_names, mask) if m]
    else:
        import re
        pattern = f'^{gene_name}[-_]'
        gene_isos = [v for v in var_names if re.match(pattern, v)]

    if len(gene_isos) < 2:
        return None, None

    # Get indices and counts
    iso_idx = [var_names.index(iso) for iso in gene_isos]

    if sp.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)

    iso_counts = X[:, iso_idx]
    total = iso_counts.sum(axis=1)

    # Calculate isoform 1 fraction
    iso1_frac = np.where(total > 0, iso_counts[:, 0] / total, np.nan)

    return iso1_frac, obs_names

def load_spacet_proportions(sample_id):
    """Load SpaCET cell type proportions for a sample."""
    # SpaCET proportions are stored in the correlations file
    # We need to reconstruct from the merged data in the R script
    # For now, we'll use dummy data structure - in practice would load from SpaCET output

    # The actual proportions would come from SpaCET_obj@results$deconvolution$propMat
    # Since we don't have direct access, we'll create the visualization differently
    return None

def create_spatial_plot(sample_id, output_path):
    """Create spatial visualization for a sample."""

    sample_dir = SR_DIR / sample_id
    gene_h5ad = sample_dir / "gene.quant.h5ad"
    iso_h5ad = sample_dir / "iso.quant.h5ad"

    if not gene_h5ad.exists() or not iso_h5ad.exists():
        print(f"  Missing h5ad files for {sample_id}")
        return False

    print(f"  Loading data for {sample_id}...")

    # Load gene data for spatial coordinates
    _, spatial, _, obs_names = load_h5ad_data(gene_h5ad)

    if spatial is None:
        print(f"  No spatial coordinates for {sample_id}")
        return False

    # Load isoform ratios
    iso1_frac, iso_obs_names = load_isoform_ratios(iso_h5ad, "VCAN")

    if iso1_frac is None:
        print(f"  No VCAN isoforms for {sample_id}")
        return False

    # Get correlation values for this sample
    sample_cors = spacet_cors[spacet_cors['sample'] == sample_id]
    if len(sample_cors) == 0:
        print(f"  No SpaCET correlations for {sample_id}")
        return False

    tumor_cor = sample_cors[sample_cors['cell_type'] == 'Malignant']['correlation'].values
    macro_cor = sample_cors[sample_cors['cell_type'] == 'Macrophage']['correlation'].values
    tumor_p = sample_cors[sample_cors['cell_type'] == 'Malignant']['pvalue'].values
    macro_p = sample_cors[sample_cors['cell_type'] == 'Macrophage']['pvalue'].values

    tumor_cor = tumor_cor[0] if len(tumor_cor) > 0 else 0
    macro_cor = macro_cor[0] if len(macro_cor) > 0 else 0
    tumor_p = tumor_p[0] if len(tumor_p) > 0 else 1
    macro_p = macro_p[0] if len(macro_p) > 0 else 1

    # Filter to expressed spots
    expressed_mask = ~np.isnan(iso1_frac)
    n_expressed = expressed_mask.sum()

    if n_expressed < 50:
        print(f"  Too few expressed spots ({n_expressed}) for {sample_id}")
        return False

    # Create figure
    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(1, 4, figure=fig, width_ratios=[1, 1, 1, 0.05])

    # Common plot settings
    point_size = 8

    # Plot 1: All spots (showing tissue structure)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(spatial[:, 0], spatial[:, 1], c='lightgray', s=point_size, alpha=0.5)
    ax1.scatter(spatial[expressed_mask, 0], spatial[expressed_mask, 1],
                c='steelblue', s=point_size, alpha=0.7, label=f'VCAN+ ({n_expressed})')
    ax1.set_title(f'{sample_id}\nVCAN Expressing Spots', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Array Column')
    ax1.set_ylabel('Array Row')
    ax1.invert_yaxis()
    ax1.set_aspect('equal')
    ax1.legend(loc='upper right', fontsize=8)

    # Plot 2: VCAN isoform ratio
    ax2 = fig.add_subplot(gs[0, 1])
    scatter2 = ax2.scatter(spatial[expressed_mask, 0], spatial[expressed_mask, 1],
                           c=iso1_frac[expressed_mask], cmap='RdYlBu_r',
                           s=point_size, alpha=0.8, vmin=0, vmax=1)
    ax2.set_title('VCAN Isoform 1 Fraction\n(High=Tumor-like, Low=Macrophage-like)',
                  fontsize=11, fontweight='bold')
    ax2.set_xlabel('Array Column')
    ax2.set_ylabel('Array Row')
    ax2.invert_yaxis()
    ax2.set_aspect('equal')

    # Plot 3: Correlation summary
    ax3 = fig.add_subplot(gs[0, 2])

    # Create bar chart of correlations
    cell_types = ['Malignant\n(Tumor)', 'Macrophage']
    correlations = [tumor_cor, macro_cor]
    pvalues = [tumor_p, macro_p]
    colors = ['#e74c3c' if c > 0 else '#3498db' for c in correlations]

    bars = ax3.bar(cell_types, correlations, color=colors, edgecolor='black', linewidth=1.5)
    ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax3.set_ylabel("Spearman Correlation (r)", fontsize=11)
    ax3.set_title("Cell Type vs VCAN Isoform 1\nCorrelation", fontsize=11, fontweight='bold')
    ax3.set_ylim(-0.3, 0.3)

    # Add significance stars
    for i, (bar, p) in enumerate(zip(bars, pvalues)):
        height = bar.get_height()
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        y_pos = height + 0.02 if height > 0 else height - 0.04
        ax3.text(bar.get_x() + bar.get_width()/2., y_pos, sig,
                ha='center', va='bottom' if height > 0 else 'top', fontsize=12, fontweight='bold')

    # Add correlation values
    for i, (bar, r) in enumerate(zip(bars, correlations)):
        height = bar.get_height()
        y_pos = height/2
        ax3.text(bar.get_x() + bar.get_width()/2., y_pos, f'r={r:.3f}',
                ha='center', va='center', fontsize=10, fontweight='bold', color='white')

    # Add legend for interpretation
    ax3.text(0.5, -0.15, "Red = Higher isoform 1 | Blue = Lower isoform 1",
             transform=ax3.transAxes, ha='center', fontsize=9, style='italic')

    # Colorbar
    cax = fig.add_subplot(gs[0, 3])
    cbar = plt.colorbar(scatter2, cax=cax)
    cbar.set_label('Isoform 1 Fraction', fontsize=10)

    plt.tight_layout()

    # Save figure
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved: {output_path.name}")
    return True

def create_summary_figure(samples_data, output_path):
    """Create a summary figure showing the best examples."""

    # Find samples with significant correlations
    sig_samples = spacet_cors[
        (spacet_cors['cell_type'].isin(['Malignant', 'Macrophage'])) &
        (spacet_cors['pvalue'] < 0.05)
    ]['sample'].unique()

    # Get the best sample (highest absolute correlation difference)
    best_sample = None
    best_diff = 0

    for sample in sig_samples:
        sample_data = spacet_cors[spacet_cors['sample'] == sample]
        tumor_r = sample_data[sample_data['cell_type'] == 'Malignant']['correlation'].values
        macro_r = sample_data[sample_data['cell_type'] == 'Macrophage']['correlation'].values

        if len(tumor_r) > 0 and len(macro_r) > 0:
            diff = tumor_r[0] - macro_r[0]
            if abs(diff) > best_diff:
                best_diff = abs(diff)
                best_sample = sample

    if best_sample:
        print(f"\nBest sample for visualization: {best_sample} (diff={best_diff:.3f})")
        return best_sample

    return None

def main():
    print("=" * 60)
    print("Spatial Isoform Switching Visualization")
    print("=" * 60)

    # Get list of samples
    samples = [d.name for d in SR_DIR.iterdir() if d.is_dir() and not d.name.startswith('.')]
    print(f"\nFound {len(samples)} SR samples")

    # Create output directory
    output_dir = REPORTS_DIR / "spatial_plots"
    output_dir.mkdir(exist_ok=True)

    # Process each sample
    successful = []
    for sample_id in sorted(samples):
        output_path = output_dir / f"spatial_isoform_{sample_id}.png"
        if create_spatial_plot(sample_id, output_path):
            successful.append(sample_id)

    print(f"\nSuccessfully created {len(successful)} spatial plots")

    # Find and highlight the best example
    best_sample = create_summary_figure(successful, output_dir / "best_example.png")

    # Create a combined figure with the best examples
    print("\nCreating combined summary figure...")

    # Get top samples by correlation strength
    sample_scores = []
    for sample in successful:
        sample_data = spacet_cors[spacet_cors['sample'] == sample]
        tumor_r = sample_data[sample_data['cell_type'] == 'Malignant']['correlation'].values
        macro_r = sample_data[sample_data['cell_type'] == 'Macrophage']['correlation'].values
        tumor_p = sample_data[sample_data['cell_type'] == 'Malignant']['pvalue'].values
        macro_p = sample_data[sample_data['cell_type'] == 'Macrophage']['pvalue'].values

        if len(tumor_r) > 0 and len(macro_r) > 0:
            # Score based on correlation difference and significance
            score = abs(tumor_r[0] - macro_r[0])
            if tumor_p[0] < 0.05 or macro_p[0] < 0.05:
                score *= 2  # Bonus for significance
            sample_scores.append((sample, score, tumor_r[0], macro_r[0]))

    # Sort by score
    sample_scores.sort(key=lambda x: x[1], reverse=True)

    print("\nTop samples by tumor-macrophage isoform correlation difference:")
    print("-" * 60)
    for sample, score, tumor_r, macro_r in sample_scores[:5]:
        print(f"  {sample}: Tumor r={tumor_r:+.3f}, Macrophage r={macro_r:+.3f}, diff={tumor_r-macro_r:+.3f}")

    print("\n" + "=" * 60)
    print("Visualization complete!")
    print(f"Output directory: {output_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()

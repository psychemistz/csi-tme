#!/usr/bin/env python3
"""
Experiment v5 - CORRECTLY IMPLEMENTED Gaussian ANNULAR

Key correction: Gaussian is centered at SENDER (d=0), not at ring center.
- σ = outer_radius / 3 (matching v7.py)
- At small radius: steep Gaussian gradient → significant weight variation in ring
- At large radius: flat Gaussian tail → nearly uniform weights in ring
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import os
from core import calculate_lambda, solve_concentration_field_MM, compute_IND_ring, compute_IND_gaussian_annular

np.random.seed(42)

OUTPUT_FOLDER = './output_gaussian'
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# =============================================================================
# Plotting Functions
# =============================================================================

def plot_weight_gradient_analysis(bandwidth=100):
    """
    Show how weight gradient within annulus changes with radius.
    Gaussian centered at d=0 (sender), σ = outer_radius / 3
    """
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Test different outer radii
    outer_radii = [100, 200, 400, 800]
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(outer_radii)))
    
    # Panel 1: Full Gaussian curves with annular regions marked
    ax = axes[0, 0]
    distances = np.linspace(0, 1000, 1000)
    
    for outer_r, col in zip(outer_radii, colors):
        sigma = outer_r / 3.0
        inner_r = outer_r - bandwidth
        
        # Full Gaussian (centered at d=0)
        gaussian = np.exp(-distances**2 / (2 * sigma**2))
        
        ax.plot(distances, gaussian, '-', color=col, linewidth=2, 
                label=f'σ={sigma:.0f}μm (outer={outer_r}μm)', alpha=0.7)
        
        # Mark annular region
        ax.axvspan(inner_r, outer_r, alpha=0.15, color=col)
    
    ax.set_xlabel('Distance from Sender (μm)', fontsize=12)
    ax.set_ylabel('Gaussian Weight', fontsize=12)
    ax.set_title('Gaussian Centered at Sender (d=0)\nAnnular regions shaded', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1000)
    
    # Panel 2: Zoom into annular regions (normalized)
    ax = axes[0, 1]
    
    for outer_r, col in zip(outer_radii, colors):
        sigma = outer_r / 3.0
        inner_r = outer_r - bandwidth
        
        # Distances within annulus
        d_annulus = np.linspace(inner_r, outer_r, 100)
        
        # Gaussian weights in annulus
        w_annulus = np.exp(-d_annulus**2 / (2 * sigma**2))
        
        # Normalize to show relative variation
        w_normalized = w_annulus / w_annulus.max()
        
        # Plot as function of position within annulus (0 = inner, 1 = outer)
        x_normalized = (d_annulus - inner_r) / bandwidth
        
        ax.plot(x_normalized, w_normalized, '-', color=col, linewidth=2.5,
                label=f'outer={outer_r}μm')
    
    ax.set_xlabel('Position in Annulus (0=inner, 1=outer)', fontsize=12)
    ax.set_ylabel('Normalized Weight', fontsize=12)
    ax.set_title('Weight Gradient WITHIN Annulus\n(Small radius = steep, Large radius = flat)', 
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Panel 3: Weight ratio (inner/outer) vs radius
    ax = axes[1, 0]
    
    test_radii = np.arange(bandwidth + 10, 1000, 10)
    ratios = []
    inner_weights = []
    outer_weights = []
    
    for outer_r in test_radii:
        sigma = outer_r / 3.0
        inner_r = outer_r - bandwidth
        
        w_inner = np.exp(-inner_r**2 / (2 * sigma**2))
        w_outer = np.exp(-outer_r**2 / (2 * sigma**2))
        
        ratios.append(w_inner / w_outer)
        inner_weights.append(w_inner)
        outer_weights.append(w_outer)
    
    ax.plot(test_radii, ratios, 'b-', linewidth=2.5)
    ax.axhline(1.0, color='red', linestyle='--', linewidth=1.5, label='Uniform (ratio=1)')
    ax.set_xlabel('Outer Radius (μm)', fontsize=12)
    ax.set_ylabel('Weight Ratio (inner/outer)', fontsize=12)
    ax.set_title('Weight Gradient Strength vs Radius\n(Higher ratio = more gradient)', 
                fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.5, 6)
    
    # Panel 4: Summary table
    ax = axes[1, 1]
    ax.axis('off')
    
    # Create table data
    table_data = []
    for outer_r in [100, 200, 400, 800]:
        sigma = outer_r / 3.0
        inner_r = outer_r - bandwidth
        w_inner = np.exp(-inner_r**2 / (2 * sigma**2))
        w_outer = np.exp(-outer_r**2 / (2 * sigma**2))
        ratio = w_inner / w_outer
        variation = (w_inner - w_outer) / w_inner * 100
        
        table_data.append([
            f'{outer_r}',
            f'{inner_r}-{outer_r}',
            f'{sigma:.0f}',
            f'{w_inner:.4f}',
            f'{w_outer:.4f}',
            f'{ratio:.2f}x',
            f'{variation:.0f}%'
        ])
    
    table = ax.table(
        cellText=table_data,
        colLabels=['Outer (μm)', 'Annulus', 'σ (μm)', 'W(inner)', 'W(outer)', 'Ratio', 'Variation'],
        loc='center',
        cellLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)
    
    ax.set_title('Weight Statistics by Radius\n(bandwidth=100μm, σ=outer/3)', 
                fontsize=13, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/weight_gradient_analysis.png", dpi=300)
    plt.close()
    print(f"Saved: {OUTPUT_FOLDER}/weight_gradient_analysis.png")


def plot_results_comparison(results_gaussian, results_ring, fractions):
    """Plot both methods"""
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(fractions)))
    
    # Left: Gaussian ANNULAR
    ax = axes[0]
    for f, col in zip(fractions, colors):
        d = results_gaussian[f]
        radii = np.array([x['outer_radius'] for x in d['ind_curve']])
        vals = np.array([x['I_ND'] for x in d['ind_curve']])
        lambda_v = d['lambda']
        
        vals_smooth = gaussian_filter1d(vals, sigma=1.5)
        
        ax.plot(radii, vals_smooth, '-', color=col, linewidth=2.5, 
                label=f'{f*100:.0f}% (λ={lambda_v:.0f}μm)')
        ax.axvline(lambda_v, color=col, linestyle=':', linewidth=1, alpha=0.5)

    ax.axhline(0, color='black', linestyle='-', alpha=0.3)
    ax.set_xlabel("Outer Radius (μm)", fontsize=13)
    ax.set_ylabel("I_ND", fontsize=13)
    ax.set_title("Gaussian ANNULAR\n(Gaussian centered at sender, σ=radius/3)", 
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.2, 1.1)
    ax.set_xlim(0, 3000)
    
    # Right: Binary Ring
    ax = axes[1]
    for f, col in zip(fractions, colors):
        d = results_ring[f]
        radii = np.array([x['outer_radius'] for x in d['ind_curve']])
        vals = np.array([x['I_ND'] for x in d['ind_curve']])
        lambda_v = d['lambda']
        
        vals_smooth = gaussian_filter1d(vals, sigma=1.5)
        
        ax.plot(radii, vals_smooth, '-', color=col, linewidth=2.5, 
                label=f'{f*100:.0f}% (λ={lambda_v:.0f}μm)')
        ax.axvline(lambda_v, color=col, linestyle=':', linewidth=1, alpha=0.5)

    ax.axhline(0, color='black', linestyle='-', alpha=0.3)
    ax.set_xlabel("Outer Radius (μm)", fontsize=13)
    ax.set_ylabel("I_ND", fontsize=13)
    ax.set_title("Binary Ring (Uniform weights)\n(bw=100μm)", 
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.2, 1.1)
    ax.set_xlim(0, 3000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/ind_comparison.png", dpi=300)
    plt.close()
    print(f"Saved: {OUTPUT_FOLDER}/ind_comparison.png")


def plot_overlay(results_gaussian, results_ring, fractions):
    """Overlay to show difference"""
    
    plt.figure(figsize=(12, 8))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(fractions)))
    
    for f, col in zip(fractions, colors):
        # Gaussian ANNULAR
        dg = results_gaussian[f]
        radii_g = np.array([x['outer_radius'] for x in dg['ind_curve']])
        vals_g = np.array([x['I_ND'] for x in dg['ind_curve']])
        vals_g_smooth = gaussian_filter1d(vals_g, sigma=1.5)
        
        # Binary Ring
        dr = results_ring[f]
        radii_r = np.array([x['outer_radius'] for x in dr['ind_curve']])
        vals_r = np.array([x['I_ND'] for x in dr['ind_curve']])
        vals_r_smooth = gaussian_filter1d(vals_r, sigma=1.5)
        
        lambda_v = dg['lambda']
        
        plt.plot(radii_g, vals_g_smooth, '-', color=col, linewidth=2.5, 
                label=f'{f*100:.0f}% Gaussian')
        plt.plot(radii_r, vals_r_smooth, '--', color=col, linewidth=2.5, 
                alpha=0.7)
        plt.axvline(lambda_v, color=col, linestyle=':', linewidth=1, alpha=0.4)

    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.xlabel("Outer Radius (μm)", fontsize=14, fontweight='bold')
    plt.ylabel("I_ND", fontsize=14, fontweight='bold')
    plt.title("Gaussian ANNULAR (solid) vs Binary Ring (dashed)\n" + 
              "Gaussian centered at sender with σ=radius/3", 
              fontsize=16, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.2, 1.1)
    plt.xlim(0, 3000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/ind_overlay.png", dpi=300)
    plt.close()
    print(f"Saved: {OUTPUT_FOLDER}/ind_overlay.png")


def plot_difference(results_gaussian, results_ring, fractions):
    """Plot the difference between methods"""
    
    plt.figure(figsize=(12, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(fractions)))
    
    for f, col in zip(fractions, colors):
        dg = results_gaussian[f]
        dr = results_ring[f]
        
        radii = np.array([x['outer_radius'] for x in dg['ind_curve']])
        vals_g = np.array([x['I_ND'] for x in dg['ind_curve']])
        vals_r = np.array([x['I_ND'] for x in dr['ind_curve']])
        
        diff = vals_g - vals_r
        diff_smooth = gaussian_filter1d(diff, sigma=1.5)
        
        plt.plot(radii, diff_smooth, '-', color=col, linewidth=2, 
                label=f'{f*100:.0f}%')

    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.xlabel("Outer Radius (μm)", fontsize=13)
    plt.ylabel("I_ND(Gaussian) - I_ND(Binary)", fontsize=13)
    plt.title("Difference: Gaussian ANNULAR - Binary Ring\n" +
              "(Positive = Gaussian gives higher I_ND)", fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 3000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/ind_difference.png", dpi=300)
    plt.close()
    print(f"Saved: {OUTPUT_FOLDER}/ind_difference.png")


# =============================================================================
# Simulation
# =============================================================================

def run_experiment():
    print("="*70)
    print("EXPERIMENT: Correctly Implemented Gaussian ANNULAR")
    print("Gaussian centered at SENDER (d=0), σ = outer_radius / 3")
    print("="*70)
    
    # Domain
    center = np.array([0, 0])
    n_total = 100000
    max_radius = 5000
    domain_area = np.pi * max_radius**2
    
    # Bio-Parameters
    k_max = 10.0
    Kd = 30.0
    D_calculated = 100.0
    
    print(f"Diffusion Coefficient (D): {D_calculated:.2f} μm²/s")
    
    # Gene Params
    BASAL_FACTOR = 0.1
    HIGH_FACTOR = 10.0
    BASAL_RESPONSIVE = 0.1
    FOLD_CHANGE = 2.0
    BANDWIDTH = 100.0
    sigma_f = 0.1
    sigma_r = 0.1
    
    n_active = 20
    n_silent = 0
    n_senders = n_active + n_silent
    
    receiver_fractions = [0.02, 0.05, 0.10, 0.20, 0.30, 0.50]
    
    # Test outer radii
    test_radii = np.arange(BANDWIDTH + 10, 3000, 25)
    
    # Generate Positions
    angles = np.random.rand(n_total) * 2 * np.pi
    radii_pos = np.sqrt(np.random.rand(n_total)) * max_radius
    all_positions_orig = np.column_stack([
        center[0] + radii_pos * np.cos(angles),
        center[1] + radii_pos * np.sin(angles)
    ])
    
    results_gaussian = {}
    results_ring = {}
    
    # Plot weight gradient analysis
    plot_weight_gradient_analysis(bandwidth=BANDWIDTH)
    
    for frac in receiver_fractions:
        print(f"\nProcessing {frac*100:.0f}% Receivers...")

        all_positions = all_positions_orig.copy()

        # Assign Cell Types
        all_indices = np.arange(n_total)
        sender_indices = np.random.choice(all_indices, n_senders, replace=False)
        
        active_indices = sender_indices[:n_active]
        silent_indices = sender_indices[n_active:]
        
        # Fix senders at center
        all_positions[active_indices] = center
        all_positions[silent_indices] = center
        
        non_sender_indices = np.setdiff1d(all_indices, sender_indices)
        n_receivers = int(n_total * frac)
        receiver_indices = np.random.choice(non_sender_indices, n_receivers, replace=False)
        
        # Factor Expression
        factor_expr = BASAL_FACTOR * np.random.lognormal(0, sigma_f, n_total)
        factor_expr[active_indices] = HIGH_FACTOR * np.random.lognormal(0, sigma_f, n_active)
        factor_expr[silent_indices] = BASAL_FACTOR * np.random.lognormal(0, sigma_f, n_silent)  
        
        # Diffusion
        n_density = n_receivers / domain_area
        sender_pos_subset = all_positions[sender_indices]
        sender_expr_subset = factor_expr[sender_indices]
        
        concentrations, lambda_val = solve_concentration_field_MM(
            sender_pos_subset, sender_expr_subset, all_positions,
            n_density, D_calculated, k_max, Kd
        )
        
        # Responsive Expression
        responsive_expr = BASAL_RESPONSIVE * np.random.lognormal(0, sigma_r, n_total)
        activation = concentrations[receiver_indices] / (Kd + concentrations[receiver_indices])
        responsive_expr[receiver_indices] = BASAL_RESPONSIVE * (1 + FOLD_CHANGE * activation)
        
        # Compute I_ND with GAUSSIAN ANNULAR
        ind_curve_gaussian = []
        for r in test_radii:
            val, n_conn = compute_IND_gaussian_annular(
                sender_indices, receiver_indices,
                all_positions, factor_expr, responsive_expr,
                outer_radius=r, bandwidth=BANDWIDTH, sigma_fraction=3.0
            )
            ind_curve_gaussian.append({'outer_radius': r, 'I_ND': val, 'n_conn': n_conn})
            
        results_gaussian[frac] = {
            'lambda': lambda_val,
            'ind_curve': ind_curve_gaussian,
            'bandwidth': BANDWIDTH
        }
        
        # Compute I_ND with BINARY RING
        ind_curve_ring = []
        for r in test_radii:
            val, n_conn = compute_IND_ring(
                sender_indices, receiver_indices,
                all_positions, factor_expr, responsive_expr,
                distance=r - BANDWIDTH/2, bandwidth=BANDWIDTH
            )
            ind_curve_ring.append({'outer_radius': r, 'I_ND': val, 'n_conn': n_conn})
            
        results_ring[frac] = {
            'lambda': lambda_val,
            'ind_curve': ind_curve_ring,
            'bandwidth': BANDWIDTH
        }
        
        print(f"  λ = {lambda_val:.0f} μm")

    # Visualizations
    plot_results_comparison(results_gaussian, results_ring, receiver_fractions)
    plot_overlay(results_gaussian, results_ring, receiver_fractions)
    plot_difference(results_gaussian, results_ring, receiver_fractions)
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"{'Frac':<8} {'λ':<10} {'Gauss@200':<12} {'Ring@200':<12} {'Diff':<10}")
    print("-"*55)
    
    for frac in receiver_fractions:
        lambda_v = results_gaussian[frac]['lambda']
        
        # Find I_ND at radius=200
        g_curve = results_gaussian[frac]['ind_curve']
        r_curve = results_ring[frac]['ind_curve']
        
        g_radii = np.array([x['outer_radius'] for x in g_curve])
        g_vals = np.array([x['I_ND'] for x in g_curve])
        r_vals = np.array([x['I_ND'] for x in r_curve])
        
        idx_200 = np.argmin(np.abs(g_radii - 200))
        
        g_val = g_vals[idx_200]
        r_val = r_vals[idx_200]
        diff = g_val - r_val
        
        print(f"{frac*100:>5.0f}%  {lambda_v:>6.0f}μm   {g_val:>8.3f}     {r_val:>8.3f}     {diff:>+8.3f}")
    
    print("-"*55)
    print("\nKey insight:")
    print("- At SMALL radius: Gaussian weights inner edge MORE → higher I_ND")
    print("- At LARGE radius: Gaussian ≈ uniform → similar to Binary")
    print(f"\nResults saved to {OUTPUT_FOLDER}/")
    
    return results_gaussian, results_ring


if __name__ == "__main__":
    results_g, results_r = run_experiment()

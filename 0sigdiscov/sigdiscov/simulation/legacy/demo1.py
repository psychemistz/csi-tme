#!/usr/bin/env python3
"""
Experiment v5: Long-Range Niche (Lambda ~ 1000um)
Changes:
1. Calculated D dynamically to force Lambda=1000um at 10% receivers.
2. Fractions: 10, 20, 30, 40, 50%.
3. Output folder: ./output
4. Plotting range extended to 2500um to visualize the larger niche.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import os
from matplotlib.colors import LinearSegmentedColormap

from core import calculate_lambda, solve_concentration_field_MM, compute_IND_ring

# Custom colormaps
white_red = LinearSegmentedColormap.from_list("white_red", ["white", "red"])

np.random.seed(42)

OUTPUT_FOLDER = './output'
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# =============================================================================
# Plotting Functions
# =============================================================================

def plot_cell_distribution(all_positions, sender_indices, receiver_indices, center, frac_str):
    """Plot spatial distribution of cells"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot receivers
    other_indices = np.setdiff1d(np.arange(len(all_positions)), 
                                  np.concatenate([sender_indices, receiver_indices]))
    ax.scatter(all_positions[other_indices, 0], all_positions[other_indices, 1],
               s=1, c='lightgray', alpha=0.3, label='Other cells')
    ax.scatter(all_positions[receiver_indices, 0], all_positions[receiver_indices, 1],
               s=2, c='blue', alpha=0.5, label='Receivers')
    
    # Plot senders (at center)
    ax.scatter(all_positions[sender_indices, 0], all_positions[sender_indices, 1],
               s=100, c='red', marker='*', edgecolors='black', linewidths=0.5, 
               label='Senders', zorder=10)
    
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Cell Distribution ({frac_str} Receivers)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/cell_distribution_{frac_str}.png", dpi=300)
    plt.close()

def plot_concentration_gradient(all_positions, concentrations, center, lambda_val, frac_str):
    """Plot concentration as spatial gradient map"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Create concentration heatmap
    scatter = ax.scatter(all_positions[:, 0], all_positions[:, 1],
                        c=concentrations, s=1, cmap='hot', alpha=0.6,
                        norm=plt.matplotlib.colors.LogNorm(vmin=concentrations[concentrations > 0].min(), 
                                                           vmax=concentrations.max()))
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Concentration (a.u.)', fontsize=12)
    
    # Plot senders at center
    ax.scatter(center[0], center[1], s=200, c='cyan', marker='*', 
               edgecolors='black', linewidths=1, label='Senders', zorder=10)
    
    # Add lambda circle
    circle = plt.Circle(center, lambda_val, fill=False, edgecolor='cyan', 
                       linewidth=2, linestyle='--', label=f'λ = {lambda_val:.0f} μm')
    ax.add_patch(circle)
    
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Concentration Gradient ({frac_str} Receivers)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/concentration_gradient_{frac_str}.png", dpi=300)
    plt.close()

def plot_results_v5(results, fractions, center, max_radius):
    
    plt.figure(figsize=(10, 8))
    
    # Use a perceptually distinct colormap
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(fractions)))
    
    for f, col in zip(fractions, colors):
        d = results[f]
        dists = np.array([x['distance'] for x in d['ind_curve']])
        vals = np.array([x['I_ND'] for x in d['ind_curve']])
        lambda_v = d['lambda']
        
        # Smooth
        vals_smooth = gaussian_filter1d(vals, sigma=1.5)
        
        plt.plot(dists, vals_smooth, '-', color=col, linewidth=3, 
                 label=f'{f*100:.0f}% Receivers (λ={lambda_v:.0f}μm)')
        
        # Mark Lambda location on the curve
        plt.axvline(lambda_v, color=col, linestyle=':', linewidth=1.5, alpha=0.6)

    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.xlabel("Distance from Senders (μm)", fontsize=14, fontweight='bold')
    plt.ylabel("I_ND (Ring, bw=100μm)", fontsize=14, fontweight='bold')
    plt.title("Diffusion Coefficient (D = 100 μm²/s)", fontsize=16, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.ylim(-1.1, 1.1)
    plt.xlim(0, 5000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/long_range_niche.png", dpi=300)
    plt.close()

# =============================================================================
# Simulation
# =============================================================================

def run_experiment_v5():
    print("="*70)
    print("EXPERIMENT v5: Target Lambda ~ 1000um")
    print("="*70)
    
    # Domain
    domain_size = 5000
    center = np.array([0, 0])
    n_total = 100000
    max_radius = 5000
    domain_area = np.pi * max_radius**2
    
    # Fixed Bio-Parameters
    k_max = 10.0
    Kd = 30.0

    # --- CALCULATE REQUIRED D ---
    # Goal: Lambda = 1000 at 10% receivers
    target_lambda = 1000.0
    ref_fraction = 0.10
    ref_n = (n_total * ref_fraction) / domain_area
    
    # Formula: lambda^2 = (D * Kd) / (n * k_max)
    # D = (lambda^2 * n * k_max) / Kd
    # D_calculated = (target_lambda**2 * ref_n * k_max) / Kd
    D_calculated = 100.0
    
    print(f"Target Lambda: {target_lambda} um at {ref_fraction*100}% density")
    print(f"Calculated Diffusion Coefficient (D): {D_calculated:.2f} um^2/s")
    
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
    
    # Updated Fractions
    receiver_fractions = [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
    
    # Extended Distance Range (to see the drop off past 1000um)
    test_distances = np.arange(BANDWIDTH/2, 5000, 25)
    
    # Generate Positions once (re-used structure)
    angles = np.random.rand(n_total) * 2 * np.pi
    radii = np.sqrt(np.random.rand(n_total)) * max_radius
    all_positions_orig = np.column_stack([
        center[0] + radii * np.cos(angles),
        center[1] + radii * np.sin(angles)
    ])

    results = {}

    for frac in receiver_fractions:
        print(f"\nProcessing {frac*100:.0f}% Receivers...")

        # Fresh copy each iteration to avoid position corruption across iterations
        all_positions = all_positions_orig.copy()

        # 1. Assign Cell Types
        all_indices = np.arange(n_total)
        sender_indices = np.random.choice(all_indices, n_senders, replace=False)

        active_indices = sender_indices[:n_active]
        silent_indices = sender_indices[n_active:]

        # FIX POSITIONS: All senders to center
        all_positions[active_indices] = center
        all_positions[silent_indices] = center
        
        non_sender_indices = np.setdiff1d(all_indices, sender_indices)
        n_receivers = int(n_total * frac)
        receiver_indices = np.random.choice(non_sender_indices, n_receivers, replace=False)
        
        # 2. Gene 1: Factor Expression
        factor_expr = BASAL_FACTOR * np.random.lognormal(0, sigma_f, n_total)
        factor_expr[active_indices] = HIGH_FACTOR * np.random.lognormal(0, sigma_f, n_active)
        factor_expr[silent_indices] = BASAL_FACTOR * np.random.lognormal(0, sigma_f, n_silent)  
        
        # 3. Diffusion with HIGH D
        n_density = n_receivers / domain_area
        sender_pos_subset = all_positions[sender_indices]
        sender_expr_subset = factor_expr[sender_indices]
        
        concentrations, lambda_val = solve_concentration_field_MM(
            sender_pos_subset, sender_expr_subset, all_positions,
            n_density, D_calculated, k_max, Kd
        )
        
        # 4. Gene 2: Responsive Expression
        responsive_expr = BASAL_RESPONSIVE * np.random.lognormal(0, sigma_r, n_total)
        activation = concentrations[receiver_indices] / (Kd + concentrations[receiver_indices])
        responsive_expr[receiver_indices] = BASAL_RESPONSIVE * (1 + FOLD_CHANGE * activation)
        
        # 5. Compute I_ND
        ind_curve = []
        for d in test_distances:
            val, n_conn = compute_IND_ring(
                sender_indices, receiver_indices,
                all_positions, factor_expr, responsive_expr,
                d, bandwidth=BANDWIDTH
            )
            ind_curve.append({'distance': d, 'I_ND': val, 'n_conn': n_conn})
            
        results[frac] = {
            'lambda': lambda_val,
            'ind_curve': ind_curve,
            'bandwidth': BANDWIDTH
        }
        
        print(f"  Lambda: {lambda_val:.0f} um")
        
        # Generate plots for each fraction
        frac_str = f"{int(frac*100)}pct"
        plot_cell_distribution(all_positions, sender_indices, receiver_indices, center, frac_str)
        plot_concentration_gradient(all_positions, concentrations, center, lambda_val, frac_str)

    # Visualization
    plot_results_v5(results, receiver_fractions, center, max_radius)
    print(f"\nDone. Results in {OUTPUT_FOLDER}")

if __name__ == "__main__":
    run_experiment_v5()

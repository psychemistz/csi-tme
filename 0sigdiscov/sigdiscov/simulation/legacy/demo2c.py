#!/usr/bin/env python3
"""
Experiment v5d: Long-Range Niche with Multiple Active Sender Positions
Changes from v5c:
1. Active senders distributed randomly across N positions (configurable)
2. Positions can be RANDOM (anywhere in domain) or FIXED (W, E, N, S, C)
3. Creates multiple signal sources with varying magnitudes
4. Silent senders spread uniformly across domain (negative controls)
5. Tests if I_ND correctly detects signaling from multiple sources

Configuration:
- N_ACTIVE_POSITIONS: Number of distinct signal source locations
- USE_RANDOM_POSITIONS: True for random placement, False for fixed W/E/N/S/C
- MIN_SEPARATION: Minimum distance between random positions (μm)
- POSITION_OFFSET: Distance from center for fixed positions (μm)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import os
from core import (calculate_lambda, compute_IND_ring, solve_concentration_field_MM_multipos,
                  get_5_positions, get_random_positions, distribute_active_senders)
from matplotlib.colors import LinearSegmentedColormap

# Custom colormaps
white_red = LinearSegmentedColormap.from_list("white_red", ["white", "red"])

np.random.seed(42)

OUTPUT_FOLDER = './output'
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Fixed concentration normalization for comparable plots across sigma conditions
VMIN_CONC = 1e-4
VMAX_CONC = 1e4

# =============================================================================
# Plotting Functions
# =============================================================================

def plot_cell_distribution(all_positions, active_indices, silent_indices, receiver_indices, 
                           position_dict, active_assignments, frac_str):
    """Plot spatial distribution of cells with active positions marked"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Combine all sender indices
    sender_indices = np.concatenate([active_indices, silent_indices])
    
    # Plot other cells
    other_indices = np.setdiff1d(np.arange(len(all_positions)), 
                                  np.concatenate([sender_indices, receiver_indices]))
    ax.scatter(all_positions[other_indices, 0], all_positions[other_indices, 1],
               s=1, c='lightgray', alpha=0.3, label='Other cells')
    
    # Plot receivers
    ax.scatter(all_positions[receiver_indices, 0], all_positions[receiver_indices, 1],
               s=2, c='blue', alpha=0.5, label='Receivers')
    
    # Plot silent senders (distributed)
    ax.scatter(all_positions[silent_indices, 0], all_positions[silent_indices, 1],
               s=20, c='red', marker='D', alpha=0.6, 
               edgecolors='black', linewidths=0.3,
               label='Silent senders', zorder=9)
    
    # Dynamic colors and markers for positions
    n_pos = len(position_dict)
    color_cycle = plt.cm.tab10(np.linspace(0, 1, max(n_pos, 10)))
    marker_cycle = ['*', '^', 'v', '>', '<', 's', 'p', 'h', 'D', 'o']
    
    for i, (pos_name, coord) in enumerate(position_dict.items()):
        n_at_pos = active_assignments[pos_name]
        if n_at_pos > 0:
            # Size proportional to number of senders (signal magnitude)
            size = 100 + n_at_pos * 30
            ax.scatter(coord[0], coord[1], s=size, c=[color_cycle[i % len(color_cycle)]], 
                      marker=marker_cycle[i % len(marker_cycle)], edgecolors='black', linewidths=1,
                      label=f'{pos_name}: {n_at_pos} senders', zorder=10)
    
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Cell Distribution - {len(position_dict)} Active Positions ({frac_str} Receivers)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/cell_distribution_5pos_{frac_str}.png", dpi=300)
    plt.close()

def plot_concentration_gradient(all_positions, concentrations, active_indices, silent_indices, 
                                position_dict, active_assignments, lambda_val, frac_str, 
                                vmin_conc=1e-3, vmax_conc=1e3):
    """Plot concentration as spatial gradient map with source positions"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Create concentration heatmap
    scatter = ax.scatter(all_positions[:, 0], all_positions[:, 1],
                        c=concentrations, s=1, cmap='hot', alpha=0.6,
                        norm=plt.matplotlib.colors.LogNorm(vmin=vmin_conc, vmax=vmax_conc))                        
    
    # Plot silent senders at their distributed positions
    ax.scatter(all_positions[silent_indices, 0], all_positions[silent_indices, 1], 
               s=50, c='cyan', marker='D', alpha=0.8,
               edgecolors='black', linewidths=0.5, label='Silent senders', zorder=9)
    
    # Dynamic colors and markers for positions (bright colors for visibility on hot colormap)
    n_pos = len(position_dict)
    bright_colors = ['lime', 'yellow', 'cyan', 'magenta', 'white', 
                     'springgreen', 'gold', 'deepskyblue', 'hotpink', 'lavender']
    marker_cycle = ['*', '^', 'v', '>', '<', 's', 'p', 'h', 'D', 'o']
    
    for i, (pos_name, coord) in enumerate(position_dict.items()):
        n_at_pos = active_assignments[pos_name]
        if n_at_pos > 0:
            size = 150 + n_at_pos * 40
            ax.scatter(coord[0], coord[1], s=size, c=bright_colors[i % len(bright_colors)], 
                      marker=marker_cycle[i % len(marker_cycle)], edgecolors='black', linewidths=1.5,
                      label=f'{pos_name}: {n_at_pos} (mag={n_at_pos}x)', zorder=10)
            
            # Add lambda circle around each active position
            circle = plt.Circle(coord, lambda_val, fill=False, 
                               edgecolor=bright_colors[i % len(bright_colors)], 
                               linewidth=1.5, linestyle='--', alpha=0.6)
            ax.add_patch(circle)
    
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Concentration Gradient - {len(position_dict)} Sources (λ={lambda_val:.0f}μm, {frac_str} Receivers)', 
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/concentration_gradient_5pos_{frac_str}.png", dpi=300)
    plt.close()

# =============================================================================
# Simulation
# =============================================================================

def run_experiment_v5d():
    print("="*70)
    print("EXPERIMENT v5d: Multiple Active Sender Positions (W, E, N, S, C)")
    print("="*70)
    
    # Domain
    domain_size = 5000
    center = np.array([0, 0])
    n_total = 100000
    max_radius = 5000
    domain_area = np.pi * max_radius**2
    
    # Fixed Bio-Parameters
    Kd = 30.0
    k_max = 10.0
   
    # --- CALCULATE REQUIRED D ---
    target_lambda = 1000.0
    ref_fraction = 0.10
    ref_n = (n_total * ref_fraction) / domain_area
    D_calculated = 100.0
    
    print(f"Target Lambda: {target_lambda} um at {ref_fraction*100}% density")
    print(f"Calculated Diffusion Coefficient (D): {D_calculated:.2f} um^2/s")
    
    # Gene Params
    BASAL_FACTOR = 0.1
    HIGH_FACTOR = 1.5
    BASAL_RESPONSIVE = 0.1
    FOLD_CHANGE = 2.0
    BANDWIDTH = 100.0
    sigma_f = 0.1
    sigma_r = 0.1 
    
    n_active = 20
    n_silent = 40
    n_senders = n_active + n_silent
    
    # ========== POSITION CONFIGURATION ==========
    N_ACTIVE_POSITIONS = 20           # Number of distinct signal source locations
    USE_RANDOM_POSITIONS = True      # True: random positions, False: fixed W/E/N/S/C
    POSITION_OFFSET = 3000.0         # For fixed mode: distance from center (μm)
    MIN_SEPARATION = 800.0           # For random mode: min distance between positions (μm)

    # Secretion rate: scales expression -> concentration
    gamma_s = 100.0    # Calibrate so C_max ~ few * Kd

    # Get positions (random or fixed)
    if USE_RANDOM_POSITIONS:
        position_dict = get_random_positions(N_ACTIVE_POSITIONS, center, max_radius, MIN_SEPARATION)
        print(f"\n{N_ACTIVE_POSITIONS} RANDOM Active Sender Positions (min_sep={MIN_SEPARATION}μm):")
    else:
        position_dict = get_5_positions(center, POSITION_OFFSET)
        print(f"\n5 FIXED Active Sender Positions (offset={POSITION_OFFSET}μm):")
    
    for name, coord in position_dict.items():
        print(f"  {name}: ({coord[0]:.0f}, {coord[1]:.0f})")
    
    # Updated Fractions
    receiver_fractions = [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
    
    # Extended Distance Range
    test_distances = np.arange(BANDWIDTH/2, 5000, 25)
    
    # Generate Positions once (re-used structure)
    angles = np.random.rand(n_total) * 2 * np.pi
    radii = np.sqrt(np.random.rand(n_total)) * max_radius
    all_positions = np.column_stack([
        center[0] + radii * np.cos(angles),
        center[1] + radii * np.sin(angles)
    ])
    
    # =========================================================================
    # FIXED SENDER SETUP (constant across all receiver fractions)
    # =========================================================================
    all_indices = np.arange(n_total)
    sender_indices = np.random.choice(all_indices, n_senders, replace=False)
    
    active_indices = sender_indices[:n_active]
    silent_indices = sender_indices[n_active:]
    
    # Distribute active senders randomly across positions ONCE
    active_assignments, sender_pos_list = distribute_active_senders(n_active, position_dict)
    
    print(f"\nActive sender distribution (FIXED for all conditions):")
    for name, count in active_assignments.items():
        if count > 0:
            print(f"  {name}: {count} senders")
    
    # Set positions for active senders ONCE
    for i, (pos_name, coord) in enumerate(sender_pos_list):
        all_positions[active_indices[i]] = coord
    
    # Silent senders remain at their original random positions
    # =========================================================================
    
    results = {}
    
    for frac in receiver_fractions:
        print(f"\nProcessing {frac*100:.0f}% Receivers...")
        
        # Only receiver selection changes per fraction
        non_sender_indices = np.setdiff1d(all_indices, sender_indices)
        n_receivers = int(n_total * frac)
        receiver_indices = np.random.choice(non_sender_indices, n_receivers, replace=False)
        
        # Gene 1: Factor Expression (regenerate to add noise, but same structure)
        factor_expr = BASAL_FACTOR * np.random.lognormal(0, sigma_f, n_total)
        factor_expr[active_indices] = HIGH_FACTOR * np.random.lognormal(0, sigma_f, n_active)
        # Silent senders keep basal expression
        
        # Diffusion with multiple sources
        n_density = n_receivers / domain_area
        sender_pos_subset = all_positions[sender_indices]
        sender_expr_subset = factor_expr[sender_indices]
        
        concentrations, lambda_val = solve_concentration_field_MM_multipos(
            sender_pos_subset, sender_expr_subset, all_positions,
            n_density, D_calculated, k_max, Kd,
            position_dict=position_dict, active_assignments=active_assignments,
            secretion_rate=gamma_s
        )
        
        # Gene 2: Responsive Expression
        responsive_expr = BASAL_RESPONSIVE * np.random.lognormal(0, sigma_r, n_total)
        activation = concentrations[receiver_indices] / (Kd + concentrations[receiver_indices])
        responsive_expr[receiver_indices] = BASAL_RESPONSIVE * (1 + FOLD_CHANGE * activation)
        
        # Compute I_ND
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
            'bandwidth': BANDWIDTH,
            'active_assignments': active_assignments.copy()
        }
        
        print(f"  Lambda: {lambda_val:.0f} um")

        # Generate plots for each fraction
        frac_str = f"{int(frac*100)}pct"
        plot_cell_distribution(all_positions, active_indices, silent_indices, 
                               receiver_indices, position_dict, active_assignments, frac_str)
        plot_concentration_gradient(all_positions, concentrations, active_indices, silent_indices, 
                                    position_dict, active_assignments, lambda_val, frac_str, 
                                    vmin_conc=VMIN_CONC, vmax_conc=VMAX_CONC)        

    # Visualization
    plot_results_v5d(results, receiver_fractions, center, max_radius, position_dict, USE_RANDOM_POSITIONS)
    print(f"\nDone. Results in {OUTPUT_FOLDER}")

def plot_results_v5d(results, fractions, center, max_radius, position_dict, use_random):
    
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
    
    pos_type = "Random" if use_random else "Fixed (W,E,N,S,C)"
    plt.title(f"Niche Detection: {len(position_dict)} Active Positions ({pos_type})", fontsize=16, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-1.1, 1.1)
    plt.xlim(0, 5000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/long_range_niche_5positions.png", dpi=300)
    plt.close()
    
    # Additional plot: Show position distribution summary
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Position locations
    ax1 = axes[0]
    n_pos = len(position_dict)
    color_cycle = plt.cm.tab10(np.linspace(0, 1, max(n_pos, 10)))
    marker_cycle = ['*', '^', 'v', '>', '<', 's', 'p', 'h', 'D', 'o']
    
    # Draw domain boundary
    circle = plt.Circle(center, max_radius, fill=False, edgecolor='gray', 
                       linewidth=2, linestyle='--', alpha=0.5)
    ax1.add_patch(circle)
    
    first_frac = fractions[0]
    assignments = results[first_frac]['active_assignments']
    
    for i, (pos_name, coord) in enumerate(position_dict.items()):
        n_at_pos = assignments[pos_name]
        size = 200 + n_at_pos * 50
        ax1.scatter(coord[0], coord[1], s=size, c=[color_cycle[i % len(color_cycle)]], 
                   marker=marker_cycle[i % len(marker_cycle)], edgecolors='black', linewidths=2,
                   label=f'{pos_name}: n={n_at_pos}')
        ax1.annotate(pos_name, (coord[0], coord[1]), textcoords="offset points", 
                    xytext=(10, 10), fontsize=10, fontweight='bold')
    
    ax1.set_xlabel('X (μm)', fontsize=12)
    ax1.set_ylabel('Y (μm)', fontsize=12)
    ax1.set_title(f'Signal Source Locations ({pos_type})', fontsize=14, fontweight='bold')
    ax1.set_aspect('equal')
    ax1.set_xlim(-max_radius*1.1, max_radius*1.1)
    ax1.set_ylim(-max_radius*1.1, max_radius*1.1)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Right: Bar chart of sender counts
    ax2 = axes[1]
    pos_names = list(position_dict.keys())
    counts = [assignments[name] for name in pos_names]
    colors_bar = [color_cycle[i % len(color_cycle)] for i in range(len(pos_names))]
    
    bars = ax2.bar(pos_names, counts, color=colors_bar, edgecolor='black', linewidth=1.5)
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        if count > 0:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, 
                    str(count), ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    ax2.set_xlabel('Position', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Number of Active Senders', fontsize=14, fontweight='bold')
    ax2.set_title('Active Sender Distribution\n(Determines Signal Magnitude)', 
                  fontsize=14, fontweight='bold')
    ax2.set_ylim(0, max(counts) + 2 if max(counts) > 0 else 5)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/active_sender_distribution.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    run_experiment_v5d()

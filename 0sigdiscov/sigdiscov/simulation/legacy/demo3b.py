#!/usr/bin/env python3
"""
Experiment v5e: Long-Range Niche with Stochastic Response Gene Expression
Changes from v5d:
1. **STOCHASTIC RESPONSE**: Receiver gene expression is now stochastic
   - Probability of responding follows Hill function of local concentration
   - Response magnitude drawn from gamma/lognormal when "on"
   - Creates realistic heterogeneity in downstream signaling
2. All v5d features retained (stochastic sender expression, multiple positions)

Response Model:
- P(respond) = P_max * C^n / (Kd_response^n + C^n)  [Hill function]
- If responding: fold_change ~ Gamma(shape, scale) or Lognormal
- Final expression = basal * (1 + stochastic_fold_change * activation)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import os
from core import (calculate_lambda, compute_IND_ring, solve_concentration_field_MM_multipos,
                  get_5_positions, get_random_positions, distribute_active_senders,
                  generate_stochastic_expression_ref, generate_stochastic_response)
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
                           position_dict, active_assignments, frac_str, expressing_mask=None):
    """Plot spatial distribution of cells with active positions marked"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    sender_indices = np.concatenate([active_indices, silent_indices])
    
    other_indices = np.setdiff1d(np.arange(len(all_positions)), 
                                  np.concatenate([sender_indices, receiver_indices]))
    ax.scatter(all_positions[other_indices, 0], all_positions[other_indices, 1],
               s=1, c='lightgray', alpha=0.3, label='Other cells')
    
    ax.scatter(all_positions[receiver_indices, 0], all_positions[receiver_indices, 1],
               s=2, c='blue', alpha=0.5, label='Receivers')
    
    ax.scatter(all_positions[silent_indices, 0], all_positions[silent_indices, 1],
               s=20, c='red', marker='D', alpha=0.6, 
               edgecolors='black', linewidths=0.3,
               label='Silent senders', zorder=9)
    
    n_pos = len(position_dict)
    color_cycle = plt.cm.tab10(np.linspace(0, 1, max(n_pos, 10)))
    marker_cycle = ['*', '^', 'v', '>', '<', 's', 'p', 'h', 'D', 'o']
    
    # Count actually expressing senders per position
    if expressing_mask is not None:
        expressing_indices = active_indices[expressing_mask]
        silent_active_indices = active_indices[~expressing_mask]
        
        # Plot non-expressing active senders (marked differently)
        if len(silent_active_indices) > 0:
            ax.scatter(all_positions[silent_active_indices, 0], 
                      all_positions[silent_active_indices, 1],
                      s=80, c='white', marker='o', edgecolors='red', linewidths=2,
                      label=f'Active-OFF ({len(silent_active_indices)})', zorder=9)
    
    for i, (pos_name, coord) in enumerate(position_dict.items()):
        n_at_pos = active_assignments[pos_name]
        if n_at_pos > 0:
            size = 100 + n_at_pos * 30
            ax.scatter(coord[0], coord[1], s=size, c=[color_cycle[i % len(color_cycle)]], 
                      marker=marker_cycle[i % len(marker_cycle)], edgecolors='black', linewidths=1,
                      label=f'{pos_name}: {n_at_pos} senders', zorder=10)
    
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    
    n_expressing = np.sum(expressing_mask) if expressing_mask is not None else len(active_indices)
    ax.set_title(f'Cell Distribution - {n_expressing}/{len(active_indices)} Active Expressing ({frac_str} Receivers)', 
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/cell_distribution_5pos_{frac_str}.png", dpi=300)
    plt.close()

def plot_concentration_gradient(all_positions, concentrations, active_indices, silent_indices, 
                                position_dict, active_assignments, lambda_val, frac_str, 
                                vmin_conc=1e-3, vmax_conc=1e3, expressing_mask=None):
    """Plot concentration as spatial gradient map with source positions"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    scatter = ax.scatter(all_positions[:, 0], all_positions[:, 1],
                        c=concentrations, s=1, cmap='hot', alpha=0.6,
                        norm=plt.matplotlib.colors.LogNorm(vmin=vmin_conc, vmax=vmax_conc))                        
    
    ax.scatter(all_positions[silent_indices, 0], all_positions[silent_indices, 1], 
               s=50, c='cyan', marker='D', alpha=0.8,
               edgecolors='black', linewidths=0.5, label='Silent senders', zorder=9)
    
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

def plot_expression_distribution(factor_expr, active_indices, expressing_mask, frac_str):
    """Plot histogram of expression levels showing stochastic distribution"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: All cells
    ax1 = axes[0]
    ax1.hist(factor_expr, bins=50, color='gray', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Factor Expression', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('All Cells - Expression Distribution', fontsize=14, fontweight='bold')
    ax1.set_yscale('log')
    
    # Right: Active senders only
    ax2 = axes[1]
    active_expr = factor_expr[active_indices]
    expressing_expr = active_expr[expressing_mask]
    silent_expr = active_expr[~expressing_mask]
    
    if len(silent_expr) > 0:
        ax2.hist(silent_expr, bins=20, color='lightcoral', alpha=0.7, 
                edgecolor='black', label=f'OFF ({len(silent_expr)})')
    if len(expressing_expr) > 0:
        ax2.hist(expressing_expr, bins=20, color='green', alpha=0.7, 
                edgecolor='black', label=f'ON ({len(expressing_expr)})')
    
    ax2.axvline(np.mean(expressing_expr) if len(expressing_expr) > 0 else 0, 
               color='darkgreen', linestyle='--', linewidth=2, label='Mean (ON)')
    ax2.set_xlabel('Factor Expression', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Active Senders - Stochastic Expression', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/expression_distribution_{frac_str}.png", dpi=300)
    plt.close()


def plot_response_distribution(responsive_expr, receiver_indices, responding_mask, 
                               response_probs, concentrations, frac_str, Kd_response):
    """Plot response gene expression showing stochastic concentration-dependent response"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Top-left: Response probability vs concentration
    ax1 = axes[0, 0]
    C = concentrations[receiver_indices]
    
    # Sort by concentration for cleaner plot
    sort_idx = np.argsort(C)
    C_sorted = C[sort_idx]
    probs_sorted = response_probs[sort_idx]
    responded_sorted = responding_mask[sort_idx]
    
    # Bin the data for visualization
    n_bins = 50
    bin_edges = np.logspace(np.log10(max(C.min(), 1e-6)), np.log10(C.max()), n_bins + 1)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    
    # Calculate observed response rate per bin
    observed_rates = []
    expected_probs = []
    for i in range(n_bins):
        mask = (C >= bin_edges[i]) & (C < bin_edges[i+1])
        if np.sum(mask) > 5:  # Only bins with enough samples
            observed_rates.append(np.mean(responding_mask[mask]))
            expected_probs.append(np.mean(response_probs[mask]))
        else:
            observed_rates.append(np.nan)
            expected_probs.append(np.nan)
    
    ax1.scatter(bin_centers, observed_rates, s=50, c='blue', alpha=0.7, label='Observed', zorder=5)
    ax1.plot(bin_centers, expected_probs, 'r-', linewidth=2, label='Expected (Hill)', zorder=4)
    ax1.axvline(Kd_response, color='green', linestyle='--', linewidth=2, label=f'Kd={Kd_response}')
    ax1.set_xscale('log')
    ax1.set_xlabel('Concentration', fontsize=12)
    ax1.set_ylabel('P(respond)', fontsize=12)
    ax1.set_title('Response Probability vs Concentration', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.05, 1.05)
    
    # Top-right: Response expression histogram
    ax2 = axes[0, 1]
    receiver_expr = responsive_expr[receiver_indices]
    responding_expr = receiver_expr[responding_mask]
    non_responding_expr = receiver_expr[~responding_mask]
    
    if len(non_responding_expr) > 0:
        ax2.hist(non_responding_expr, bins=30, color='lightcoral', alpha=0.7, 
                edgecolor='black', label=f'Non-responding ({len(non_responding_expr)})')
    if len(responding_expr) > 0:
        ax2.hist(responding_expr, bins=30, color='green', alpha=0.7, 
                edgecolor='black', label=f'Responding ({len(responding_expr)})')
    
    ax2.set_xlabel('Response Gene Expression', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Response Expression Distribution', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.set_yscale('log')
    
    # Bottom-left: Expression vs concentration scatter
    ax3 = axes[1, 0]
    ax3.scatter(C[~responding_mask], receiver_expr[~responding_mask], 
               s=5, c='lightcoral', alpha=0.3, label='Non-responding')
    ax3.scatter(C[responding_mask], receiver_expr[responding_mask], 
               s=5, c='green', alpha=0.5, label='Responding')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.axvline(Kd_response, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
    ax3.set_xlabel('Concentration', fontsize=12)
    ax3.set_ylabel('Response Expression', fontsize=12)
    ax3.set_title('Expression vs Concentration (Stochastic)', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10, loc='lower right')
    ax3.grid(True, alpha=0.3)
    
    # Bottom-right: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    n_responding = np.sum(responding_mask)
    n_total_recv = len(responding_mask)
    mean_prob = np.mean(response_probs)
    
    # Handle edge case when no cells respond
    if len(responding_expr) > 0:
        resp_mean = responding_expr.mean()
        resp_std = responding_expr.std()
        resp_cv = resp_std / resp_mean if resp_mean > 0 else 0
    else:
        resp_mean = resp_std = resp_cv = 0
    
    summary_text = f"""
    STOCHASTIC RESPONSE SUMMARY
    ===========================
    
    Total Receivers: {n_total_recv:,}
    Responding: {n_responding:,} ({100*n_responding/n_total_recv:.1f}%)
    Non-responding: {n_total_recv - n_responding:,} ({100*(1-n_responding/n_total_recv):.1f}%)
    
    Mean P(respond): {mean_prob:.3f}
    
    Concentration Range:
      Min: {C.min():.2e}
      Max: {C.max():.2e}
      Median: {np.median(C):.2e}
    
    Response Expression (responding only):
      Mean: {resp_mean:.3f}
      Std: {resp_std:.3f}
      CV: {resp_cv:.2f}
    """
    
    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/response_distribution_{frac_str}.png", dpi=300)
    plt.close()


# =============================================================================
# Simulation
# =============================================================================

def run_experiment_v5e():
    print("="*70)
    print("EXPERIMENT v5e: Stochastic Response Gene Expression")
    print("="*70)
    
    # Domain
    domain_size = 5000
    center = np.array([0, 0])
    n_total = 100000
    max_radius = 5000
    domain_area = np.pi * max_radius**2
    
    # Fixed Bio-Parameters
    Kd = 1.0
    k_max = 10.0
   
    # --- CALCULATE REQUIRED D ---
    target_lambda = 1000.0
    ref_fraction = 0.10
    ref_n = (n_total * ref_fraction) / domain_area
    D_calculated = 100.0
    
    print(f"Target Lambda: {target_lambda} um at {ref_fraction*100}% density")
    print(f"Calculated Diffusion Coefficient (D): {D_calculated:.2f} um^2/s")
    
    # Gene Params (matching reference model notation)
    F_BASAL = 0.1             # F_basal: basal factor expression
    F_HIGH = 10               # F_high: high factor expression when "on"
    R_BASAL = 0.1             # R_basal: basal response expression
    FOLD_CHANGE = 10.0         # FC: fold change for response
    BANDWIDTH = 20.0
    
    # Lognormal sigma parameters (separate for expressing vs basal)
    SIGMA_F = 0.1             # σ_f: noise for expressing senders
    SIGMA_F_BASAL = 0.1       # σ_f,b: noise for non-expressing/basal senders
    SIGMA_R = 0.1             # σ_r: noise for responding receivers  
    SIGMA_R_BASAL = 0.1       # σ_r,b: noise for non-responding receivers
    
    # ========== STOCHASTIC SENDER EXPRESSION PARAMETERS ==========
    P_EXPRESS = 0.9           # p_s: Probability that an active sender is "on"
    
    print(f"\nStochastic Sender Expression:")
    print(f"  Q_i = γ_s * F_i")
    print(f"  F_i = S_i * F_high * LN(0,σ_f²) + (1-S_i) * F_basal * LN(0,σ_f,b²)")
    print(f"  P(express) p_s = {P_EXPRESS:.0%}")
    print(f"  F_high = {F_HIGH}, F_basal = {F_BASAL}")
    print(f"  σ_f = {SIGMA_F}, σ_f,b = {SIGMA_F_BASAL}")
    # ============================================================
    
    # ========== STOCHASTIC RESPONSE PARAMETERS (HYBRID MODEL) ==========
    P_RESPOND = 1.0           # p_r: Maximum probability of responding
    RESPONSE_HILL = 1.0       # n: Hill coefficient (1.0 = standard M-M)
    
    # Secretion rate: scales expression -> concentration
    gamma_s = 1.0    # Calibrate so C_max ~ few * Kd
    
    print(f"\nStochastic Response:")
    print(f"  B_j ~ Bernoulli(p_r * C^n / (Kd^n + C^n))")
    print(f"  R_j = B_j * R_basal * (1 + FC * Act) * LN(0,σ_r²)")
    print(f"       + (1-B_j) * R_basal * LN(0,σ_r,b²)")
    print(f"  P_max = {P_RESPOND:.0%}, Hill n = {RESPONSE_HILL}")
    print(f"  R_basal = {R_BASAL}, FC = {FOLD_CHANGE}")
    print(f"  σ_r = {SIGMA_R}, σ_r,b = {SIGMA_R_BASAL}")
    print(f"  Kd = {Kd}, Secretion rate = {gamma_s}")
    # ==================================================================
    
    n_active = 200
    n_silent = 100
    n_senders = n_active + n_silent
    
    # ========== POSITION CONFIGURATION ==========
    N_ACTIVE_POSITIONS = 1           
    USE_RANDOM_POSITIONS = True     
    POSITION_OFFSET = 3000.0         
    MIN_SEPARATION = 0.0           
    
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
    
    # Generate Positions once
    angles = np.random.rand(n_total) * 2 * np.pi
    radii = np.sqrt(np.random.rand(n_total)) * max_radius
    all_positions = np.column_stack([
        center[0] + radii * np.cos(angles),
        center[1] + radii * np.sin(angles)
    ])
    
    # =========================================================================
    # FIXED SENDER SETUP
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
    
    results = {}
    
    for frac in receiver_fractions:
        print(f"\nProcessing {frac*100:.0f}% Receivers...")
        
        non_sender_indices = np.setdiff1d(all_indices, sender_indices)
        n_receivers = int(n_total * frac)
        receiver_indices = np.random.choice(non_sender_indices, n_receivers, replace=False)
        
        # =====================================================================
        # STOCHASTIC FACTOR EXPRESSION (Sender) - Reference Model
        # F_i = S_i * F_high * LN(0,σ_f²) + (1-S_i) * F_basal * LN(0,σ_f,b²)
        # =====================================================================
        factor_expr, expressing_mask = generate_stochastic_expression_ref(
            n_senders=n_total,
            n_active=n_active,
            active_indices=active_indices,
            F_basal=F_BASAL,
            F_high=F_HIGH,
            sigma_f=SIGMA_F,
            sigma_f_basal=SIGMA_F_BASAL,
            p_express=P_EXPRESS
        )
        
        n_expressing = np.sum(expressing_mask)
        print(f"  Stochastic Sender: {n_expressing}/{n_active} active senders expressing")
        
        if n_expressing > 0:
            expr_values = factor_expr[active_indices[expressing_mask]]
            print(f"  Sender expression range: {expr_values.min():.3f} - {expr_values.max():.3f}")
        # =====================================================================
        
        # Diffusion with multiple sources
        n_density = n_receivers / domain_area
        sender_pos_subset = all_positions[sender_indices]
        sender_expr_subset = factor_expr[sender_indices]
        
        concentrations, lambda_val = solve_concentration_field_MM_multipos(
            sender_pos_subset, sender_expr_subset, all_positions,
            n_density, D_calculated, k_max, Kd, p_r=P_RESPOND,
            position_dict=position_dict, active_assignments=active_assignments,
            secretion_rate=gamma_s, active_threshold=0.2
        )
        
        # Print concentration range for debugging
        C_recv = concentrations[receiver_indices]
        print(f"  Concentration at receivers: min={C_recv.min():.2e}, max={C_recv.max():.2e}, median={np.median(C_recv):.2e}")
        
        # =====================================================================
        # STOCHASTIC RESPONSE EXPRESSION (Receiver) - Hybrid Model
        # B_j ~ Bernoulli(p_max * C^n / (Kd^n + C^n))
        # R_j = B_j * R_basal * (1 + FC * Act) * LN(0,σ_r²) 
        #     + (1-B_j) * R_basal * LN(0,σ_r,b²)
        # =====================================================================
        responsive_expr, responding_mask, response_probs = generate_stochastic_response(
            n_total=n_total,
            receiver_indices=receiver_indices,
            concentrations=concentrations,
            basal_responsive=R_BASAL,
            fold_change=FOLD_CHANGE,
            sigma_r=SIGMA_R,
            sigma_r_basal=SIGMA_R_BASAL,
            Kd=Kd,
            p_respond_max=P_RESPOND,
            hill_coef=RESPONSE_HILL
        )
        
        n_responding = np.sum(responding_mask)
        print(f"  Stochastic Response: {n_responding}/{n_receivers} receivers responding ({100*n_responding/n_receivers:.1f}%)")
        
        if n_responding > 0:
            responding_expr = responsive_expr[receiver_indices[responding_mask]]
            print(f"  Response expression range: {responding_expr.min():.3f} - {responding_expr.max():.3f}")
        # =====================================================================
        
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
            'active_assignments': active_assignments.copy(),
            'n_expressing': n_expressing,
            'expressing_mask': expressing_mask.copy(),
            'n_responding': n_responding,
            'responding_mask': responding_mask.copy(),
            'response_probs': response_probs.copy()
        }
        
        print(f"  Lambda: {lambda_val:.0f} um")

        # Generate plots for each fraction
        frac_str = f"{int(frac*100)}pct"
        plot_cell_distribution(all_positions, active_indices, silent_indices, 
                               receiver_indices, position_dict, active_assignments, frac_str,
                               expressing_mask=expressing_mask)
        plot_concentration_gradient(all_positions, concentrations, active_indices, silent_indices, 
                                    position_dict, active_assignments, lambda_val, frac_str, 
                                    vmin_conc=VMIN_CONC, vmax_conc=VMAX_CONC,
                                    expressing_mask=expressing_mask)
        plot_expression_distribution(factor_expr, active_indices, expressing_mask, frac_str)
        
        # NEW: Plot response distribution
        plot_response_distribution(responsive_expr, receiver_indices, responding_mask,
                                   response_probs, concentrations, frac_str, Kd)

    # Visualization
    plot_results_v5e(results, receiver_fractions, center, max_radius, position_dict, 
                     USE_RANDOM_POSITIONS, P_EXPRESS, SIGMA_F, 
                     P_RESPOND, RESPONSE_HILL, Kd, SIGMA_R)
    print(f"\nDone. Results in {OUTPUT_FOLDER}")

def plot_results_v5e(results, fractions, center, max_radius, position_dict, use_random,
                     p_express, sigma_f, p_respond, response_hill, kd, sigma_r):
    
    plt.figure(figsize=(10, 8))
    
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(fractions)))
    
    for f, col in zip(fractions, colors):
        d = results[f]
        dists = np.array([x['distance'] for x in d['ind_curve']])
        vals = np.array([x['I_ND'] for x in d['ind_curve']])
        lambda_v = d['lambda']
        n_expr = d['n_expressing']
        n_resp = d['n_responding']
        
        vals_smooth = gaussian_filter1d(vals, sigma=1.5)
        
        plt.plot(dists, vals_smooth, '-', color=col, linewidth=3, 
                 label=f'{f*100:.0f}% (λ={lambda_v:.0f}μm, on={n_expr}, resp={n_resp})')
        
        plt.axvline(lambda_v, color=col, linestyle=':', linewidth=1.5, alpha=0.6)

    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.xlabel("Distance from Senders (μm)", fontsize=14, fontweight='bold')
    plt.ylabel("I_ND (Ring, bw=100μm)", fontsize=14, fontweight='bold')
    
    pos_type = "Random" if use_random else "Fixed (W,E,N,S,C)"
    plt.title(f"Niche Detection: {len(position_dict)} Positions ({pos_type})\n"
              f"Sender: P(on)={p_express:.0%}, σ_f={sigma_f} | "
              f"Response: P_r={p_respond:.0%}, Hill={response_hill}, Kd={kd}, σ_r={sigma_r}", 
              fontsize=13, fontweight='bold')
    plt.legend(fontsize=8, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-1.1, 1.1)
    plt.xlim(0, 5000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/long_range_niche_stochastic_response.png", dpi=300)
    plt.close()
    
    # Position distribution summary
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax1 = axes[0]
    n_pos = len(position_dict)
    color_cycle = plt.cm.tab10(np.linspace(0, 1, max(n_pos, 10)))
    marker_cycle = ['*', '^', 'v', '>', '<', 's', 'p', 'h', 'D', 'o']
    
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
    
    ax2 = axes[1]
    pos_names = list(position_dict.keys())
    counts = [assignments[name] for name in pos_names]
    colors_bar = [color_cycle[i % len(color_cycle)] for i in range(len(pos_names))]
    
    bars = ax2.bar(pos_names, counts, color=colors_bar, edgecolor='black', linewidth=1.5)
    
    for bar, count in zip(bars, counts):
        if count > 0:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, 
                    str(count), ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    ax2.set_xlabel('Position', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Number of Active Senders', fontsize=14, fontweight='bold')
    ax2.set_title(f'Active Sender Distribution\n(Stochastic: ~{p_express:.0%} expressing)', 
                  fontsize=14, fontweight='bold')
    ax2.set_ylim(0, max(counts) + 2 if max(counts) > 0 else 5)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/active_sender_distribution.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    run_experiment_v5e()

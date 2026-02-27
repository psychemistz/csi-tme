#!/usr/bin/env python3
"""
Experiment v5 with VST Normalization
=====================================
Changes from original:
1. Added VST normalization to mimic real CosMx data
2. Non-expressing cells now have small negative values after VST
3. Raw counts generated first, then transformed
4. Diffusion uses raw expression; I_ND uses VST-normalized expression

VST approaches included:
- Simple: log1p + centering
- Pearson residuals (sctransform-like)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial
from scipy.ndimage import gaussian_filter1d
import pandas as pd
import os
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
# VST Normalization Functions
# =============================================================================

def apply_vst_log1p(raw_expr):
    """
    Simple VST normalization: log1p transform + centering.
    Non-expressing cells will have small negative values.
    
    This mimics the behavior seen in real CosMx VST-normalized data.
    """
    # Log transform (like log1p used in many VST implementations)
    log_expr = np.log1p(raw_expr)
    
    # Center by global mean (key step that creates negative values for low expressors)
    vst_expr = log_expr - np.mean(log_expr)
    
    return vst_expr


def apply_vst_pearson(raw_expr, theta=100, clip_value=10):
    """
    Simulate sctransform-style VST using Pearson residuals.
    
    Parameters:
    -----------
    raw_expr : array
        Raw count data
    theta : float
        Dispersion parameter for negative binomial model
    clip_value : float
        Clip extreme residuals (like sctransform does)
    
    Returns:
    --------
    Pearson residuals (VST-normalized expression)
    """
    # Estimate expected value (simplified: use mean)
    mu = np.mean(raw_expr) + 1e-10  # avoid division by zero
    
    # Pearson residual: (observed - expected) / sqrt(variance)
    # For negative binomial, variance ~ mu + mu^2/theta
    variance = mu + (mu ** 2) / theta
    
    residuals = (raw_expr - mu) / np.sqrt(variance)
    
    # Clip extreme values (like sctransform does)
    residuals = np.clip(residuals, -clip_value, clip_value)
    
    return residuals


def apply_vst_shifted_log(raw_expr, pseudocount=1.0):
    """
    Shifted log normalization that ensures non-expressors get negative values.
    
    This approach:
    1. Adds pseudocount
    2. Log transforms
    3. Shifts so that the pseudocount baseline maps to a small negative value
    """
    log_expr = np.log(raw_expr + pseudocount)
    
    # Shift so that cells with raw_expr=0 have negative values
    # The baseline (raw=0) becomes log(pseudocount) - shift
    baseline = np.log(pseudocount)
    shift = baseline + 0.5  # This makes raw=0 map to -0.5
    
    vst_expr = log_expr - shift
    
    return vst_expr


# =============================================================================
# Core Functions
# =============================================================================

def calculate_lambda(D, n_receivers, k_max, Kd):
    if n_receivers <= 0: return np.inf
    return np.sqrt(D * Kd / (n_receivers * k_max))


def solve_concentration_field_MM(sender_positions, sender_factor_raw, 
                                  cell_positions, n_receivers, 
                                  D, k_max=300.0, Kd=5.0,
                                  active_threshold=1.0):
    """
    Calculate factor concentration field using RAW expression values.
    
    Note: Uses raw expression (not VST) for biophysical diffusion calculation.
    The active_threshold should be set based on raw expression scale.
    """
    lambda_val = calculate_lambda(D, n_receivers, k_max, Kd)
    concentrations = np.zeros(len(cell_positions))
    
    # Identify active senders based on RAW expression
    active_mask = sender_factor_raw > active_threshold
    active_pos = sender_positions[active_mask]
    active_expr = sender_factor_raw[active_mask]
    
    if len(active_pos) > 0:
        # All active senders are at center
        center_pos = active_pos[0]
        total_factor = np.sum(active_expr)
        
        for i, pos in enumerate(cell_positions):
            r = np.linalg.norm(pos - center_pos)
            if r < 1e-3:
                concentrations[i] = total_factor * 100 
            else:
                concentrations[i] = total_factor * np.exp(-r / lambda_val) / np.sqrt(r)
                
    return concentrations, lambda_val


def compute_IND_ring(sender_indices, receiver_indices, 
                     all_positions, factor_expr_vst, responsive_expr_vst,
                     distance, bandwidth=100):
    """
    Compute I_ND using a RING neighborhood.
    
    Note: Uses VST-normalized expression values.
    The z-score normalization in the original is now less critical
    since VST already centers the data, but we keep it for consistency.
    """
    # 1. Global Normalization (on VST data)
    # Since VST data is already centered, this mainly scales by std
    mu_f, sigma_f = np.mean(factor_expr_vst), np.std(factor_expr_vst)
    mu_r, sigma_r = np.mean(responsive_expr_vst), np.std(responsive_expr_vst)
    
    z_s = (factor_expr_vst[sender_indices] - mu_f) / (sigma_f + 1e-10)
    z_r = (responsive_expr_vst[receiver_indices] - mu_r) / (sigma_r + 1e-10)
    
    sender_pos = all_positions[sender_indices]
    receiver_pos = all_positions[receiver_indices]
    
    # 2. Build Weight Matrix (Ring)
    dists = spatial.distance.cdist(sender_pos, receiver_pos)
    
    half_bw = bandwidth / 2.0
    lower = distance - half_bw
    upper = distance + half_bw
    
    W = ((dists > lower) & (dists <= upper)).astype(float)
    
    # Row-normalize
    row_sums = W.sum(axis=1, keepdims=True)
    
    # Check connections
    total_connections = np.sum(W)
    if total_connections == 0:
        return 0.0, 0
        
    row_sums[row_sums == 0] = 1.0 
    W_tilde = W / row_sums
    
    # 3. Compute Metrics
    mean_neighbor_z = np.dot(W_tilde, z_r)
    numerator = np.dot(z_s, mean_neighbor_z)
    norm_z_s = np.linalg.norm(z_s)
    norm_W_z_r = np.linalg.norm(mean_neighbor_z)
    
    if norm_z_s > 1e-10 and norm_W_z_r > 1e-10:
        I_ND = numerator / (norm_z_s * norm_W_z_r)
    else:
        I_ND = 0.0
    
    return I_ND, int(total_connections)


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


def plot_concentration_gradient(all_positions, concentrations, center, lambda_val, frac_str, vmin_conc=1e-3, vmax_conc=1e3):
    """Plot concentration as spatial gradient map with fixed normalization"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Create concentration heatmap with fixed normalization
    scatter = ax.scatter(all_positions[:, 0], all_positions[:, 1],
                        c=concentrations, s=1, cmap='hot', alpha=0.6,
                        norm=plt.matplotlib.colors.LogNorm(vmin=vmin_conc, vmax=vmax_conc))
    
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


def plot_expression_distributions(factor_raw, factor_vst, responsive_raw, responsive_vst,
                                   sender_indices, receiver_indices, frac_str):
    """
    Diagnostic plot showing raw vs VST expression distributions.
    Highlights that non-expressors have negative VST values.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Get cell type masks
    other_indices = np.setdiff1d(np.arange(len(factor_raw)), 
                                  np.concatenate([sender_indices, receiver_indices]))
    
    # --- Factor Expression ---
    # Raw
    ax = axes[0, 0]
    ax.hist(factor_raw[other_indices], bins=50, alpha=0.5, label='Other cells', color='gray', density=True)
    ax.hist(factor_raw[sender_indices], bins=20, alpha=0.7, label='Senders', color='red', density=True)
    ax.set_xlabel('Raw Expression')
    ax.set_ylabel('Density')
    ax.set_title('Factor Gene - Raw Counts')
    ax.legend()
    ax.set_xlim(0, np.percentile(factor_raw, 99))
    
    # VST
    ax = axes[0, 1]
    ax.hist(factor_vst[other_indices], bins=50, alpha=0.5, label='Other cells', color='gray', density=True)
    ax.hist(factor_vst[sender_indices], bins=20, alpha=0.7, label='Senders', color='red', density=True)
    ax.axvline(0, color='black', linestyle='--', linewidth=1, label='Zero')
    ax.set_xlabel('VST Expression')
    ax.set_ylabel('Density')
    ax.set_title('Factor Gene - VST Normalized')
    ax.legend()
    
    # --- Responsive Expression ---
    # Raw
    ax = axes[1, 0]
    ax.hist(responsive_raw[other_indices], bins=50, alpha=0.5, label='Other cells', color='gray', density=True)
    ax.hist(responsive_raw[receiver_indices], bins=50, alpha=0.7, label='Receivers', color='blue', density=True)
    ax.set_xlabel('Raw Expression')
    ax.set_ylabel('Density')
    ax.set_title('Responsive Gene - Raw Counts')
    ax.legend()
    ax.set_xlim(0, np.percentile(responsive_raw, 99))
    
    # VST
    ax = axes[1, 1]
    ax.hist(responsive_vst[other_indices], bins=50, alpha=0.5, label='Other cells', color='gray', density=True)
    ax.hist(responsive_vst[receiver_indices], bins=50, alpha=0.7, label='Receivers', color='blue', density=True)
    ax.axvline(0, color='black', linestyle='--', linewidth=1, label='Zero')
    ax.set_xlabel('VST Expression')
    ax.set_ylabel('Density')
    ax.set_title('Responsive Gene - VST Normalized')
    ax.legend()
    
    plt.suptitle(f'Expression Distributions ({frac_str} Receivers)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/expression_dist_{frac_str}.png", dpi=300)
    plt.close()
    
    # Print summary statistics
    print(f"    Factor VST - Non-senders: mean={np.mean(factor_vst[other_indices]):.3f}, "
          f"min={np.min(factor_vst[other_indices]):.3f}, max={np.max(factor_vst[other_indices]):.3f}")
    print(f"    Factor VST - Senders: mean={np.mean(factor_vst[sender_indices]):.3f}, "
          f"min={np.min(factor_vst[sender_indices]):.3f}, max={np.max(factor_vst[sender_indices]):.3f}")


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
    plt.title("Diffusion Coefficient (D = 100 μm²/s) - VST Normalized", fontsize=16, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.ylim(-1.1, 1.1)
    plt.xlim(0, 5000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/long_range_niche_vst.png", dpi=300)
    plt.close()


# =============================================================================
# Simulation
# =============================================================================

def run_experiment_v5():
    print("="*70)
    print("EXPERIMENT v5 with VST Normalization")
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
    target_lambda = 1535.0
    ref_fraction = 0.10
    ref_n = (n_total * ref_fraction) / domain_area
    D_calculated = 100.0
    
    print(f"Target Lambda: {target_lambda} um at {ref_fraction*100}% density")
    print(f"Calculated Diffusion Coefficient (D): {D_calculated:.2f} um^2/s")
    
    # Gene Params (RAW expression scale)
    # These define the raw count distributions before VST
    BASAL_FACTOR = 0.5       # Low basal expression (will become negative after VST)
    HIGH_FACTOR = 50.0       # High expression for active senders
    BASAL_RESPONSIVE = 0.5   # Low basal for responsive gene
    FOLD_CHANGE = 5.0        # Activation fold change
    BANDWIDTH = 100.0
    sigma_f = 0.8            # Increased noise for more realistic distribution
    sigma_r = 0.8
    
    # Zero-inflation parameters (fraction of cells with zero expression)
    # This creates more realistic CosMx-like data where many cells don't express a gene
    ZERO_INFLATE_FACTOR = 0.7   # 70% of non-sender cells have zero factor expression
    ZERO_INFLATE_RESPONSIVE = 0.5  # 50% of non-receiver cells have zero responsive expression
    
    # Threshold for identifying active senders (in RAW expression scale)
    ACTIVE_THRESHOLD = 5.0
    
    # VST method to use: 'log1p', 'pearson', or 'shifted_log'
    VST_METHOD = 'log1p'
    print(f"VST Method: {VST_METHOD}")
    
    n_active = 20
    n_silent = 0
    n_senders = n_active + n_silent
    
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
    
    results = {}
    
    for frac in receiver_fractions:
        print(f"\nProcessing {frac*100:.0f}% Receivers...")
        
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
        
        # =====================================================================
        # 2. Gene 1: Factor Expression (RAW counts first with zero-inflation)
        # =====================================================================
        factor_raw = BASAL_FACTOR * np.random.lognormal(0, sigma_f, n_total)
        
        # Zero-inflate: set a fraction of non-sender cells to zero
        zero_mask_factor = np.random.rand(n_total) < ZERO_INFLATE_FACTOR
        zero_mask_factor[sender_indices] = False  # Never zero-inflate senders
        factor_raw[zero_mask_factor] = 0.0
        
        # Active senders: high expression
        factor_raw[active_indices] = HIGH_FACTOR * np.random.lognormal(0, sigma_f, n_active)
        factor_raw[silent_indices] = BASAL_FACTOR * np.random.lognormal(0, sigma_f, n_silent)
        
        # Apply VST normalization
        if VST_METHOD == 'log1p':
            factor_vst = apply_vst_log1p(factor_raw)
        elif VST_METHOD == 'pearson':
            factor_vst = apply_vst_pearson(factor_raw)
        elif VST_METHOD == 'shifted_log':
            factor_vst = apply_vst_shifted_log(factor_raw)
        else:
            raise ValueError(f"Unknown VST method: {VST_METHOD}")
        
        # =====================================================================
        # 3. Diffusion (uses RAW expression for biophysics)
        # =====================================================================
        n_density = n_receivers / domain_area
        sender_pos_subset = all_positions[sender_indices]
        sender_expr_subset = factor_raw[sender_indices]  # Use RAW for diffusion
        
        concentrations, lambda_val = solve_concentration_field_MM(
            sender_pos_subset, sender_expr_subset, all_positions,
            n_density, D_calculated, k_max, Kd,
            active_threshold=ACTIVE_THRESHOLD
        )
        
        # =====================================================================
        # 4. Gene 2: Responsive Expression (RAW counts first with zero-inflation)
        # =====================================================================
        responsive_raw = BASAL_RESPONSIVE * np.random.lognormal(0, sigma_r, n_total)
        
        # Zero-inflate: set a fraction of non-receiver cells to zero
        zero_mask_responsive = np.random.rand(n_total) < ZERO_INFLATE_RESPONSIVE
        zero_mask_responsive[receiver_indices] = False  # Never zero-inflate receivers (they respond)
        responsive_raw[zero_mask_responsive] = 0.0
        
        # Receivers respond to concentration
        activation = concentrations[receiver_indices] / (Kd + concentrations[receiver_indices])
        responsive_raw[receiver_indices] = BASAL_RESPONSIVE * (1 + FOLD_CHANGE * activation)
        
        # Apply VST normalization
        if VST_METHOD == 'log1p':
            responsive_vst = apply_vst_log1p(responsive_raw)
        elif VST_METHOD == 'pearson':
            responsive_vst = apply_vst_pearson(responsive_raw)
        elif VST_METHOD == 'shifted_log':
            responsive_vst = apply_vst_shifted_log(responsive_raw)
        
        # =====================================================================
        # 5. Compute I_ND (uses VST-normalized expression)
        # =====================================================================
        ind_curve = []
        for d in test_distances:
            val, n_conn = compute_IND_ring(
                sender_indices, receiver_indices,
                all_positions, factor_vst, responsive_vst,  # VST normalized
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
        plot_concentration_gradient(all_positions, concentrations, center, lambda_val, frac_str, 
                                   vmin_conc=VMIN_CONC, vmax_conc=VMAX_CONC)
        
        # New: Expression distribution diagnostic plot
        plot_expression_distributions(factor_raw, factor_vst, responsive_raw, responsive_vst,
                                      sender_indices, receiver_indices, frac_str)

    # Visualization
    plot_results_v5(results, receiver_fractions, center, max_radius)
    print(f"\nDone. Results in {OUTPUT_FOLDER}")


if __name__ == "__main__":
    run_experiment_v5()

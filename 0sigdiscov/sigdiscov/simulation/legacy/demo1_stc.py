#!/usr/bin/env python3
"""
Stochastic Diffusion-Consumption Model
Based on Oyler-Yaniv et al., 2017, Cell

Key features:
1. Stochastic sending: S_i ~ Bernoulli(p_s)
2. Stochastic receiving: B_j ~ Bernoulli(p_r)
3. Effective receiver density: n_r_eff = n_r * p_r
4. Four sigma parameters: sigma_f, sigma_f_b, sigma_r, sigma_r_b

Factor expression (mixture):
  F_i = S_i * F_high * LogN(0, sigma_f^2) + (1-S_i) * F_basal * LogN(0, sigma_f_b^2)

Response expression (mixture):
  R_j = B_j * R_basal * (1 + FC * Act_j) * LogN(0, sigma_r^2) 
      + (1-B_j) * R_basal * LogN(0, sigma_r_b^2)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial
from scipy.ndimage import gaussian_filter1d
import os
from matplotlib.colors import LinearSegmentedColormap

# Custom colormaps
white_red = LinearSegmentedColormap.from_list("white_red", ["white", "red"])

np.random.seed(42)

OUTPUT_FOLDER = './output'
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Fixed concentration normalization for comparable plots across sigma conditions
VMIN_CONC = 1e-3
VMAX_CONC = 1e3

# =============================================================================
# Core Functions
# =============================================================================

def calculate_lambda(D, n_receivers_eff, k_max, Kd):
    """
    Calculate characteristic length scale.
    lambda_eff = sqrt(D * Kd / (n_r_eff * k_max))
    """
    if n_receivers_eff <= 0: 
        return np.inf
    return np.sqrt(D * Kd / (n_receivers_eff * k_max))


def solve_concentration_field_MM(sender_positions, sender_factor_expr, 
                                  cell_positions, n_receivers_eff, 
                                  D, k_max=300.0, Kd=5.0):
    """
    Calculate factor concentration field using effective receiver density.
    
    C(r) = sum_i (Q_i / sqrt(r_i)) * exp(-r_i / lambda_eff)
    
    Parameters:
    -----------
    n_receivers_eff : float
        Effective receiver density = n_r * p_r (accounts for stochastic receiving)
    """
    lambda_val = calculate_lambda(D, n_receivers_eff, k_max, Kd)
    concentrations = np.zeros(len(cell_positions))
    
    # Identify active senders (those with high expression)
    active_mask = sender_factor_expr > 1.0
    active_pos = sender_positions[active_mask]
    active_expr = sender_factor_expr[active_mask]
    
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
                     all_positions, factor_expr, responsive_expr,
                     distance, bandwidth=100):
    """
    Compute I_ND using a RING neighborhood.
    """
    # 1. Global Normalization
    mu_f, std_f = np.mean(factor_expr), np.std(factor_expr)
    mu_r, std_r = np.mean(responsive_expr), np.std(responsive_expr)
    
    z_s = (factor_expr[sender_indices] - mu_f) / (std_f + 1e-10)
    z_r = (responsive_expr[receiver_indices] - mu_r) / (std_r + 1e-10)
    
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
# Stochastic Expression Functions
# =============================================================================

def generate_factor_expression(n_total, sender_indices, p_s, 
                                F_high, F_basal, sigma_f, sigma_f_b):
    """
    Generate factor expression using mixture model.
    
    F_i = S_i * F_high * LogN(0, sigma_f^2) + (1-S_i) * F_basal * LogN(0, sigma_f_b^2)
    
    where S_i ~ Bernoulli(p_s) for sender cells, S_i = 0 for non-senders
    """
    factor_expr = np.zeros(n_total)
    
    # Non-sender cells: always basal expression
    non_sender_mask = np.ones(n_total, dtype=bool)
    non_sender_mask[sender_indices] = False
    n_non_senders = np.sum(non_sender_mask)
    factor_expr[non_sender_mask] = F_basal * np.random.lognormal(0, sigma_f_b, n_non_senders)
    
    # Sender cells: stochastic S_i ~ Bernoulli(p_s)
    n_senders = len(sender_indices)
    S = np.random.binomial(1, p_s, n_senders)  # Bernoulli draws
    
    # Mixture: S_i * high + (1-S_i) * basal
    high_expr = F_high * np.random.lognormal(0, sigma_f, n_senders)
    basal_expr = F_basal * np.random.lognormal(0, sigma_f_b, n_senders)
    factor_expr[sender_indices] = S * high_expr + (1 - S) * basal_expr
    
    # Return active sender indices (where S_i = 1)
    active_sender_indices = sender_indices[S == 1]
    
    return factor_expr, active_sender_indices


def generate_response_expression(n_total, receiver_indices, p_r, concentrations,
                                  R_basal, FC, Kd, sigma_r, sigma_r_b):
    """
    Generate response expression using mixture model.
    
    R_j = B_j * R_basal * (1 + FC * Act_j) * LogN(0, sigma_r^2) 
        + (1-B_j) * R_basal * LogN(0, sigma_r_b^2)
    
    where B_j ~ Bernoulli(p_r) for receiver cells, B_j = 0 for non-receivers
    Act_j = C_j / (Kd + C_j)
    """
    responsive_expr = np.zeros(n_total)
    
    # Non-receiver cells: always basal expression
    non_receiver_mask = np.ones(n_total, dtype=bool)
    non_receiver_mask[receiver_indices] = False
    n_non_receivers = np.sum(non_receiver_mask)
    responsive_expr[non_receiver_mask] = R_basal * np.random.lognormal(0, sigma_r_b, n_non_receivers)
    
    # Receiver cells: stochastic B_j ~ Bernoulli(p_r)
    n_receivers = len(receiver_indices)
    B = np.random.binomial(1, p_r, n_receivers)  # Bernoulli draws
    
    # Calculate activation for receivers
    C_receivers = concentrations[receiver_indices]
    activation = C_receivers / (Kd + C_receivers)
    
    # Mixture: B_j * activated + (1-B_j) * basal
    activated_expr = R_basal * (1 + FC * activation) * np.random.lognormal(0, sigma_r, n_receivers)
    basal_expr = R_basal * np.random.lognormal(0, sigma_r_b, n_receivers)
    responsive_expr[receiver_indices] = B * activated_expr + (1 - B) * basal_expr
    
    # Return responding receiver indices (where B_j = 1)
    responding_receiver_indices = receiver_indices[B == 1]
    
    return responsive_expr, responding_receiver_indices


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_cell_distribution(all_positions, sender_indices, active_sender_indices,
                           receiver_indices, responding_receiver_indices, 
                           center, frac_str, p_s, p_r):
    """Plot spatial distribution of cells with active/responding status"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Other cells (non-sender, non-receiver)
    other_indices = np.setdiff1d(np.arange(len(all_positions)), 
                                  np.concatenate([sender_indices, receiver_indices]))
    ax.scatter(all_positions[other_indices, 0], all_positions[other_indices, 1],
               s=1, c='lightgray', alpha=0.3, label='Other cells')
    
    # Non-responding receivers
    non_responding = np.setdiff1d(receiver_indices, responding_receiver_indices)
    ax.scatter(all_positions[non_responding, 0], all_positions[non_responding, 1],
               s=2, c='lightblue', alpha=0.3, label=f'Non-responding receivers')
    
    # Responding receivers
    ax.scatter(all_positions[responding_receiver_indices, 0], 
               all_positions[responding_receiver_indices, 1],
               s=2, c='blue', alpha=0.5, label=f'Responding receivers (B=1)')
    
    # Inactive senders
    inactive_senders = np.setdiff1d(sender_indices, active_sender_indices)
    if len(inactive_senders) > 0:
        ax.scatter(all_positions[inactive_senders, 0], all_positions[inactive_senders, 1],
                   s=50, c='pink', marker='*', edgecolors='black', linewidths=0.5, 
                   label='Inactive senders (S=0)', zorder=9)
    
    # Active senders
    ax.scatter(all_positions[active_sender_indices, 0], all_positions[active_sender_indices, 1],
               s=100, c='red', marker='*', edgecolors='black', linewidths=0.5, 
               label='Active senders (S=1)', zorder=10)
    
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Cell Distribution ({frac_str} Receivers, p_s={p_s}, p_r={p_r})', 
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/cell_distribution_{frac_str}.png", dpi=300)
    plt.close()


def plot_concentration_gradient(all_positions, concentrations, center, lambda_val, 
                                 frac_str, vmin_conc=1e-3, vmax_conc=1e3):
    """Plot concentration as spatial gradient map with fixed normalization"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Create concentration heatmap with fixed normalization
    scatter = ax.scatter(all_positions[:, 0], all_positions[:, 1],
                        c=concentrations, s=1, cmap='hot', alpha=0.6,
                        norm=plt.matplotlib.colors.LogNorm(vmin=vmin_conc, vmax=vmax_conc))
    
    # Colorbar removed for side-by-side comparison
    
    # Plot senders at center
    ax.scatter(center[0], center[1], s=200, c='cyan', marker='*', 
               edgecolors='black', linewidths=1, label='Senders', zorder=10)
    
    # Add lambda circle
    circle = plt.Circle(center, lambda_val, fill=False, edgecolor='cyan', 
                       linewidth=2, linestyle='--', label=f'λ_eff = {lambda_val:.0f} μm')
    ax.add_patch(circle)
    
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Concentration Gradient ({frac_str})', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/concentration_gradient_{frac_str}.png", dpi=300)
    plt.close()


def plot_results(results, fractions, p_s, p_r, sigma_params):
    """Plot I_ND curves for all receiver fractions"""
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
                 label=f'{f*100:.0f}% Receivers (λ_eff={lambda_v:.0f}μm)')
        
        # Mark Lambda location on the curve
        plt.axvline(lambda_v, color=col, linestyle=':', linewidth=1.5, alpha=0.6)

    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.xlabel("Distance from Senders (μm)", fontsize=14, fontweight='bold')
    plt.ylabel("I_ND (Ring, bw=100μm)", fontsize=14, fontweight='bold')
    
    sigma_str = f"σ_f={sigma_params['sigma_f']}, σ_f_b={sigma_params['sigma_f_b']}, σ_r={sigma_params['sigma_r']}, σ_r_b={sigma_params['sigma_r_b']}"
    plt.title(f"Stochastic Model (p_s={p_s}, p_r={p_r})\n{sigma_str}", 
              fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.ylim(-1.1, 1.1)
    plt.xlim(0, 5000)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}/ind_curves.png", dpi=300)
    plt.close()


# =============================================================================
# Simulation
# =============================================================================

def run_stochastic_experiment(sigma_f=0.1, sigma_f_b=0.1, sigma_r=0.1, sigma_r_b=0.1,
                               p_s=1.0, p_r=1.0):
    """
    Run stochastic diffusion-consumption experiment.
    
    Parameters:
    -----------
    sigma_f : float
        Log-normal sigma for active sender expression
    sigma_f_b : float
        Log-normal sigma for basal factor expression
    sigma_r : float
        Log-normal sigma for responding receiver expression
    sigma_r_b : float
        Log-normal sigma for basal response expression
    p_s : float
        Probability of sender being active (S_i ~ Bernoulli(p_s))
    p_r : float
        Probability of receiver responding (B_j ~ Bernoulli(p_r))
    """
    print("="*70)
    print("STOCHASTIC DIFFUSION-CONSUMPTION MODEL")
    print("="*70)
    print(f"Parameters:")
    print(f"  p_s = {p_s}, p_r = {p_r}")
    print(f"  sigma_f = {sigma_f}, sigma_f_b = {sigma_f_b}")
    print(f"  sigma_r = {sigma_r}, sigma_r_b = {sigma_r_b}")
    print("="*70)
    
    # Domain
    center = np.array([0, 0])
    n_total = 100000
    max_radius = 5000
    domain_area = np.pi * max_radius**2
    
    # Fixed Bio-Parameters
    k_max = 10.0
    Kd = 30.0
    D = 100.0  # Diffusion coefficient
    
    # Gene Params
    F_BASAL = 0.1
    F_HIGH = 10.0
    R_BASAL = 0.1
    FOLD_CHANGE = 2.0
    BANDWIDTH = 100.0
    
    # Number of potential senders (placed at center)
    n_senders = 20
    
    # Receiver fractions to test
    receiver_fractions = [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
    
    # Distance range for I_ND computation
    test_distances = np.arange(BANDWIDTH/2, 5000, 25)
    
    # Generate positions once (re-used structure)
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
        
        # Fix sender positions to center
        all_positions[sender_indices] = center
        
        # Select receivers from non-senders
        non_sender_indices = np.setdiff1d(all_indices, sender_indices)
        n_receivers = int(n_total * frac)
        receiver_indices = np.random.choice(non_sender_indices, n_receivers, replace=False)
        
        # 2. Factor Expression (mixture model with stochastic S_i)
        factor_expr, active_sender_indices = generate_factor_expression(
            n_total, sender_indices, p_s,
            F_HIGH, F_BASAL, sigma_f, sigma_f_b
        )
        
        # 3. Calculate effective receiver density
        # n_r_eff = n_r * p_r (only responding receivers consume)
        n_density = n_receivers / domain_area
        n_density_eff = n_density * p_r
        
        # 4. Concentration field (using effective density for lambda)
        concentrations, lambda_val = solve_concentration_field_MM(
            all_positions[sender_indices], factor_expr[sender_indices], 
            all_positions, n_density_eff, D, k_max, Kd
        )
        
        # 5. Response Expression (mixture model with stochastic B_j)
        responsive_expr, responding_receiver_indices = generate_response_expression(
            n_total, receiver_indices, p_r, concentrations,
            R_BASAL, FOLD_CHANGE, Kd, sigma_r, sigma_r_b
        )
        
        # 6. Compute I_ND curve
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
            'n_active_senders': len(active_sender_indices),
            'n_responding_receivers': len(responding_receiver_indices)
        }
        
        print(f"  λ_eff: {lambda_val:.0f} μm")
        print(f"  Active senders: {len(active_sender_indices)}/{n_senders}")
        print(f"  Responding receivers: {len(responding_receiver_indices)}/{n_receivers}")
        
        # Generate plots for each fraction
        frac_str = f"{int(frac*100)}pct"
        plot_cell_distribution(all_positions, sender_indices, active_sender_indices,
                               receiver_indices, responding_receiver_indices,
                               center, frac_str, p_s, p_r)
        plot_concentration_gradient(all_positions, concentrations, center, lambda_val, 
                                    frac_str, vmin_conc=VMIN_CONC, vmax_conc=VMAX_CONC)

    # Summary plot
    sigma_params = {
        'sigma_f': sigma_f, 'sigma_f_b': sigma_f_b,
        'sigma_r': sigma_r, 'sigma_r_b': sigma_r_b
    }
    plot_results(results, receiver_fractions, p_s, p_r, sigma_params)
    
    print(f"\nDone. Results in {OUTPUT_FOLDER}")
    return results


if __name__ == "__main__":
    # Example: Run with default parameters
    # You can modify these to test different sigma conditions:
    # (sigma_f, sigma_f_b, sigma_r, sigma_r_b)
    # Condition 1: (1, 1, 1, 1)
    # Condition 2: (0.1, 0.1, 0.1, 0.1)
    # Condition 3: (1, 0.1, 1, 0.1)
    # Condition 4: (0.1, 1, 0.1, 1)
    
    results = run_stochastic_experiment(
        sigma_f=0.1,
        sigma_f_b=0.1,
        sigma_r=0.01,
        sigma_r_b=0.1,
        p_s=1.0,
        p_r=1.0
    )

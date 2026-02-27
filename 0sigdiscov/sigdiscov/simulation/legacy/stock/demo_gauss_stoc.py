#!/usr/bin/env python3
"""
Stochastic Sender Expression Model for I_ND Analysis

Models temporal fluctuations in ligand expression:
1. Transcriptional bursting (fast, Poisson-like)
2. Cell-to-cell variability (slow, correlated)
3. Environmental noise (intermediate)

Tests effects on I_ND detection:
- Static vs dynamic expression
- Different timescales of fluctuation
- Time-averaged vs snapshot measurements
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial, ndimage, stats
from tqdm import tqdm
import matplotlib.gridspec as gridspec

np.random.seed(42)

# =============================================================================
# Stochastic Expression Models
# =============================================================================

def generate_static_expression(n_cells, mean_expr, noise_cv=0.2):
    """
    Static expression: lognormal distribution (current model)
    Represents: stable expression state over measurement window
    """
    return np.random.lognormal(np.log(mean_expr), noise_cv, n_cells)

def generate_bursting_expression(n_cells, mean_expr, burst_freq=0.3, burst_size=5.0):
    """
    Transcriptional bursting model (Poisson bursts)
    
    Models: Stochastic ON/OFF switching of gene expression
    - burst_freq: fraction of cells in "burst" state
    - burst_size: fold-increase during burst
    
    Timescale: minutes to hours (fast fluctuations)
    """
    base_expr = np.random.lognormal(np.log(mean_expr/2), 0.15, n_cells)
    
    # Random cells undergo bursting
    bursting = np.random.rand(n_cells) < burst_freq
    expression = base_expr.copy()
    expression[bursting] *= burst_size
    
    return expression

def generate_correlated_noise_expression(n_cells, positions, mean_expr, 
                                        spatial_scale=200, temporal_scale=0.3):
    """
    Spatially correlated noise (e.g., from local microenvironment)
    
    Models: Local tissue conditions affecting multiple nearby cells
    - spatial_scale: distance over which noise is correlated (μm)
    - temporal_scale: strength of temporal fluctuations
    
    Timescale: hours to days (slow fluctuations)
    """
    base_expr = np.random.lognormal(np.log(mean_expr), 0.15, n_cells)
    
    # Generate spatial noise field
    distances = spatial.distance_matrix(positions, positions)
    # Gaussian correlation: Cov(i,j) = exp(-d²/scale²)
    correlation = np.exp(-(distances / spatial_scale)**2)
    
    # Generate correlated random noise
    try:
        # Ensure positive definite by adding small diagonal
        correlation_stable = correlation + np.eye(n_cells) * 1e-6
        L = np.linalg.cholesky(correlation_stable)
        noise = L @ np.random.randn(n_cells)
    except np.linalg.LinAlgError:
        # Fallback to uncorrelated noise
        noise = np.random.randn(n_cells)
    
    # Apply noise multiplicatively
    expression = base_expr * np.exp(temporal_scale * noise)
    
    return expression

def generate_multistate_expression(n_cells, mean_expr, n_states=3):
    """
    Multi-state expression model (e.g., cell cycle phases)
    
    Models: Discrete expression states (G1, S, G2/M phases)
    - n_states: number of discrete expression levels
    
    Timescale: hours (cell cycle)
    """
    states = np.random.choice(n_states, n_cells)
    state_multipliers = np.linspace(0.5, 2.0, n_states)
    
    expression = np.zeros(n_cells)
    for state_id in range(n_states):
        mask = states == state_id
        n_in_state = np.sum(mask)
        if n_in_state > 0:
            expression[mask] = np.random.lognormal(
                np.log(mean_expr * state_multipliers[state_id]), 
                0.15, 
                n_in_state
            )
    
    return expression

def generate_temporal_trajectory(n_cells, positions, mean_expr, n_timepoints=10,
                                timescale_type='mixed'):
    """
    Generate time series of expression across multiple timepoints
    
    Args:
        timescale_type: 'fast' (bursting), 'slow' (correlated), 'mixed'
    
    Returns:
        trajectory: (n_timepoints, n_cells) array
    """
    trajectory = np.zeros((n_timepoints, n_cells))
    
    if timescale_type == 'fast':
        # Independent snapshots (fast fluctuations)
        for t in range(n_timepoints):
            trajectory[t] = generate_bursting_expression(n_cells, mean_expr)
    
    elif timescale_type == 'slow':
        # Correlated across time (slow drift)
        base = generate_static_expression(n_cells, mean_expr)
        for t in range(n_timepoints):
            # Small random walk
            drift = np.random.randn(n_cells) * 0.1
            trajectory[t] = base * np.exp(drift * t / n_timepoints)
    
    elif timescale_type == 'mixed':
        # Both fast and slow components
        base = generate_correlated_noise_expression(n_cells, positions, mean_expr)
        for t in range(n_timepoints):
            # Slow drift + fast bursting
            bursting = np.random.rand(n_cells) < 0.2
            fast_component = np.ones(n_cells)
            fast_component[bursting] *= 3.0
            trajectory[t] = base * fast_component
    
    return trajectory

# =============================================================================
# Simulation Framework
# =============================================================================

def solve_cytokine_gaussian_decay(sender_positions, sender_expressions,
                                  domain_size=(2000, 2000), 
                                  grid_spacing=20, sigma=100):
    """Simulate cytokine diffusion with Gaussian decay"""
    nx, ny = int(domain_size[0]/grid_spacing), int(domain_size[1]/grid_spacing)
    concentration = np.zeros((ny, nx))
    
    for iy in range(ny):
        for ix in range(nx):
            x, y = ix * grid_spacing, iy * grid_spacing
            grid_pos = np.array([x, y])
            distances = np.linalg.norm(sender_positions - grid_pos, axis=1)
            contributions = sender_expressions * np.exp(-(distances / sigma)**2)
            concentration[iy, ix] = np.sum(contributions)
    
    return concentration

def sample_concentration(positions, concentration_field, grid_spacing, domain_size):
    """Sample concentration with bilinear interpolation"""
    concentrations = np.zeros(len(positions))
    ny, nx = concentration_field.shape
    
    for i, pos in enumerate(positions):
        x_idx = np.clip(pos[0] / grid_spacing, 0, nx - 1.001)
        y_idx = np.clip(pos[1] / grid_spacing, 0, ny - 1.001)
        
        x0, y0 = int(np.floor(x_idx)), int(np.floor(y_idx))
        x1, y1 = min(x0 + 1, nx - 1), min(y0 + 1, ny - 1)
        dx, dy = x_idx - x0, y_idx - y0
        
        c00, c10 = concentration_field[y0, x0], concentration_field[y0, x1]
        c01, c11 = concentration_field[y1, x0], concentration_field[y1, x1]
        
        concentrations[i] = (c00 * (1-dx) * (1-dy) + c10 * dx * (1-dy) +
                           c01 * (1-dx) * dy + c11 * dx * dy)
    
    return concentrations

def compute_IND(sender_positions, sender_expression, 
                receiver_positions, receiver_expression,
                distance, bandwidth=15):
    """Calculate I_ND(d)"""
    z_s = (sender_expression - np.mean(sender_expression)) / (np.std(sender_expression) + 1e-10)
    z_r = (receiver_expression - np.mean(receiver_expression)) / (np.std(receiver_expression) + 1e-10)
    
    distances = spatial.distance_matrix(sender_positions, receiver_positions)
    W = np.zeros((len(sender_positions), len(receiver_positions)))
    W[(distances > distance - bandwidth) & (distances <= distance + bandwidth)] = 1.0
    
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    W_norm = W / row_sums
    
    numerator = np.dot(z_s, np.dot(W_norm, z_r))
    norm_z_s = np.linalg.norm(z_s)
    W_z_r = np.dot(W_norm, z_r)
    norm_W_z_r = np.linalg.norm(W_z_r)
    
    if norm_z_s > 0 and norm_W_z_r > 0:
        I_ND = numerator / (norm_z_s * norm_W_z_r)
    else:
        I_ND = 0
    
    return I_ND

def compute_IND_curve(sender_positions, sender_expr, receiver_positions, receiver_expr,
                     test_distances, bandwidth=15):
    """Compute I_ND for all test distances"""
    results = []
    for dist in test_distances:
        I_ND = compute_IND(sender_positions, sender_expr,
                          receiver_positions, receiver_expr,
                          dist, bandwidth)
        results.append({'distance': dist, 'I_ND': I_ND})
    return results

# =============================================================================
# Main Analysis
# =============================================================================

def run_stochastic_analysis():
    print("=" * 70)
    print("STOCHASTIC SENDER EXPRESSION ANALYSIS")
    print("=" * 70)
    
    domain_size = (2000, 2000)
    grid_spacing = 20
    sigma = 100
    test_distances = np.arange(20, 601, 20)
    
    # ==========================================================================
    # Setup sender clusters
    # ==========================================================================
    print("\n1. Setting up sender clusters...")
    
    n_clusters = 4
    cells_per_cluster = 20
    margin = 200
    
    cluster_centers = [
        [margin, margin],
        [domain_size[0] - margin, margin],
        [margin, domain_size[1] - margin],
        [domain_size[0] - margin, domain_size[1] - margin]
    ]
    
    cluster_expression_means = [1.0, 2.0, 3.5, 5.0]
    
    sender_positions = []
    for center in cluster_centers:
        cluster_pos = np.array(center) + np.random.randn(cells_per_cluster, 2) * 30
        sender_positions.append(cluster_pos)
    sender_positions = np.vstack(sender_positions)
    
    # Setup receivers
    n_receivers_per_cluster = 200
    max_distance = 800
    receiver_positions = []
    
    for center in cluster_centers:
        angles = np.random.rand(n_receivers_per_cluster) * 2 * np.pi
        radii = np.random.rand(n_receivers_per_cluster) * max_distance
        receivers = np.column_stack([
            center[0] + radii * np.cos(angles),
            center[1] + radii * np.sin(angles)
        ])
        receiver_positions.append(receivers)
    
    receiver_positions = np.vstack(receiver_positions)
    receiver_positions[:, 0] = np.clip(receiver_positions[:, 0], 0, domain_size[0])
    receiver_positions[:, 1] = np.clip(receiver_positions[:, 1], 0, domain_size[1])
    
    print(f"   {len(sender_positions)} senders, {len(receiver_positions)} receivers")
    
    # ==========================================================================
    # Test different stochastic expression models
    # ==========================================================================
    print("\n2. Testing stochastic expression models...")
    
    expression_models = {}
    
    # Model 1: Static (baseline)
    print("   Static expression (baseline)...")
    sender_expr_static = []
    for mean in cluster_expression_means:
        sender_expr_static.append(generate_static_expression(cells_per_cluster, mean))
    sender_expr_static = np.concatenate(sender_expr_static)
    expression_models['Static (baseline)'] = sender_expr_static
    
    # Model 2: Transcriptional bursting
    print("   Transcriptional bursting...")
    sender_expr_bursting = []
    for mean in cluster_expression_means:
        sender_expr_bursting.append(
            generate_bursting_expression(cells_per_cluster, mean, burst_freq=0.3, burst_size=5.0)
        )
    sender_expr_bursting = np.concatenate(sender_expr_bursting)
    expression_models['Transcriptional Bursting'] = sender_expr_bursting
    
    # Model 3: Correlated noise (each cluster separately)
    print("   Spatially correlated noise...")
    sender_expr_correlated = []
    for i, (center, mean) in enumerate(zip(cluster_centers, cluster_expression_means)):
        cluster_positions = sender_positions[i*cells_per_cluster:(i+1)*cells_per_cluster]
        sender_expr_correlated.append(
            generate_correlated_noise_expression(
                cells_per_cluster, cluster_positions, mean,
                spatial_scale=100, temporal_scale=0.4
            )
        )
    sender_expr_correlated = np.concatenate(sender_expr_correlated)
    expression_models['Correlated Noise'] = sender_expr_correlated
    
    # Model 4: Multi-state
    print("   Multi-state (cell cycle)...")
    sender_expr_multistate = []
    for mean in cluster_expression_means:
        sender_expr_multistate.append(
            generate_multistate_expression(cells_per_cluster, mean, n_states=3)
        )
    sender_expr_multistate = np.concatenate(sender_expr_multistate)
    expression_models['Multi-state (3 phases)'] = sender_expr_multistate
    
    # ==========================================================================
    # Compute I_ND for each model
    # ==========================================================================
    print("\n3. Computing I_ND curves for each model...")
    
    all_results = {}
    
    for model_name, sender_expr in expression_models.items():
        print(f"   {model_name}...")
        
        # Generate concentration field
        conc_field = solve_cytokine_gaussian_decay(
            sender_positions, sender_expr,
            domain_size=domain_size, grid_spacing=grid_spacing, sigma=sigma
        )
        
        # Sample and generate receiver response
        receiver_conc = sample_concentration(
            receiver_positions, conc_field, grid_spacing, domain_size
        )
        
        Kd = 2.0
        receptor_activation = receiver_conc / (Kd + receiver_conc)
        receiver_expr = (1 + 10 * receptor_activation) * \
                       np.random.lognormal(0, 0.15, len(receiver_positions))
        
        # Compute I_ND curve
        results = compute_IND_curve(
            sender_positions, sender_expr,
            receiver_positions, receiver_expr,
            test_distances
        )
        
        all_results[model_name] = {
            'results': results,
            'sender_expr': sender_expr,
            'conc_field': conc_field,
            'receiver_expr': receiver_expr
        }
    
    # ==========================================================================
    # Temporal analysis: Multiple snapshots
    # ==========================================================================
    print("\n4. Temporal analysis: multiple snapshots...")
    
    n_timepoints = 10
    timescale_types = ['fast', 'slow', 'mixed']
    temporal_results = {}
    
    for ts_type in timescale_types:
        print(f"   Timescale: {ts_type}...")
        
        # Generate temporal trajectories for each cluster
        trajectories = []
        for i, mean in enumerate(cluster_expression_means):
            cluster_positions = sender_positions[i*cells_per_cluster:(i+1)*cells_per_cluster]
            traj = generate_temporal_trajectory(
                cells_per_cluster, cluster_positions, mean,
                n_timepoints=n_timepoints, timescale_type=ts_type
            )
            trajectories.append(traj)
        
        # Combine all clusters
        full_trajectory = np.hstack(trajectories)  # (n_timepoints, n_senders)
        
        # Compute I_ND at each timepoint
        I_ND_over_time = []
        for t in range(n_timepoints):
            sender_expr_t = full_trajectory[t]
            
            # Generate concentration and response
            conc_field_t = solve_cytokine_gaussian_decay(
                sender_positions, sender_expr_t,
                domain_size=domain_size, grid_spacing=grid_spacing, sigma=sigma
            )
            
            receiver_conc_t = sample_concentration(
                receiver_positions, conc_field_t, grid_spacing, domain_size
            )
            
            receptor_activation_t = receiver_conc_t / (Kd + receiver_conc_t)
            receiver_expr_t = (1 + 10 * receptor_activation_t) * \
                            np.random.lognormal(0, 0.15, len(receiver_positions))
            
            # Compute I_ND at peak distance (100 μm)
            I_ND_100 = compute_IND(
                sender_positions, sender_expr_t,
                receiver_positions, receiver_expr_t,
                distance=100, bandwidth=15
            )
            
            I_ND_over_time.append(I_ND_100)
        
        temporal_results[ts_type] = {
            'trajectory': full_trajectory,
            'I_ND_over_time': I_ND_over_time,
            'mean_I_ND': np.mean(I_ND_over_time),
            'std_I_ND': np.std(I_ND_over_time)
        }
    
    # ==========================================================================
    # Visualization
    # ==========================================================================
    print("\n5. Creating visualizations...")
    
    fig = plt.figure(figsize=(20, 12))
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.35)
    
    # ========== Panel 1: Expression distributions ==========
    ax1 = fig.add_subplot(gs[0, 0])
    
    for model_name, data in all_results.items():
        sender_expr = data['sender_expr']
        ax1.hist(sender_expr, bins=30, alpha=0.4, label=model_name, density=True)
    
    ax1.set_xlabel('Sender Expression', fontsize=10)
    ax1.set_ylabel('Probability Density', fontsize=10)
    ax1.set_title('(A) Expression Distributions\nDifferent Stochastic Models', 
                 fontsize=11, fontweight='bold')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # ========== Panel 2: Coefficient of Variation ==========
    ax2 = fig.add_subplot(gs[0, 1])
    
    model_names = []
    cvs = []
    colors_cv = []
    
    for model_name, data in all_results.items():
        sender_expr = data['sender_expr']
        cv = np.std(sender_expr) / np.mean(sender_expr)
        model_names.append(model_name.split('(')[0].strip())
        cvs.append(cv)
        colors_cv.append('blue' if 'Static' in model_name else 'orange')
    
    bars = ax2.bar(range(len(model_names)), cvs, color=colors_cv, alpha=0.7, edgecolor='black')
    ax2.set_xticks(range(len(model_names)))
    ax2.set_xticklabels(model_names, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Coefficient of Variation', fontsize=10)
    ax2.set_title('(B) Expression Noise Level\nCV = σ/μ', 
                 fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add values on bars
    for bar, cv in zip(bars, cvs):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{cv:.2f}', ha='center', va='bottom', fontsize=9)
    
    # ========== Panel 3: Spatial distribution of expression ==========
    ax3 = fig.add_subplot(gs[0, 2])
    
    # Show bursting model as example
    sender_expr_burst = all_results['Transcriptional Bursting']['sender_expr']
    scatter = ax3.scatter(sender_positions[:, 0], sender_positions[:, 1],
                         c=sender_expr_burst, s=150, marker='*',
                         cmap='Reds', edgecolors='black', linewidths=1,
                         vmin=0, vmax=np.percentile(sender_expr_burst, 95))
    
    ax3.set_xlim(0, domain_size[0])
    ax3.set_ylim(0, domain_size[1])
    ax3.set_xlabel('X (μm)', fontsize=10)
    ax3.set_ylabel('Y (μm)', fontsize=10)
    ax3.set_title('(C) Spatial Pattern\nTranscriptional Bursting Example',
                 fontsize=11, fontweight='bold')
    ax3.set_aspect('equal')
    plt.colorbar(scatter, ax=ax3, label='Expression', pad=0.02)
    
    # ========== Panel 4: Example concentration fields ==========
    ax4 = fig.add_subplot(gs[0, 3])
    
    # Show static vs bursting concentration fields
    conc_static = all_results['Static (baseline)']['conc_field']
    conc_burst = all_results['Transcriptional Bursting']['conc_field']
    
    # Compute difference
    conc_diff = conc_burst - conc_static
    
    im4 = ax4.imshow(conc_diff, cmap='RdBu_r', origin='lower',
                    extent=[0, domain_size[0], 0, domain_size[1]],
                    vmin=-np.percentile(np.abs(conc_diff), 95),
                    vmax=np.percentile(np.abs(conc_diff), 95))
    
    ax4.set_xlabel('X (μm)', fontsize=10)
    ax4.set_ylabel('Y (μm)', fontsize=10)
    ax4.set_title('(D) Concentration Difference\nBursting - Static',
                 fontsize=11, fontweight='bold')
    plt.colorbar(im4, ax=ax4, label='ΔConcentration', pad=0.02)
    
    # ========== Panel 5-8: I_ND curves for each model ==========
    colors_models = ['blue', 'red', 'green', 'orange']
    
    ax5 = fig.add_subplot(gs[1, :])
    
    for (model_name, data), color in zip(all_results.items(), colors_models):
        results = data['results']
        distances = [r['distance'] for r in results]
        I_ND_vals = [r['I_ND'] for r in results]
        
        ax5.plot(distances, I_ND_vals, '-o', linewidth=2.5, markersize=5,
                label=model_name, color=color, alpha=0.8)
    
    ax5.axvline(sigma, color='black', linestyle=':', linewidth=2,
               alpha=0.5, label=f'σ = {sigma} μm')
    ax5.axhline(0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
    
    ax5.set_xlabel('Distance (μm)', fontsize=11)
    ax5.set_ylabel('I_ND(d)', fontsize=11)
    ax5.set_title('(E) I_ND Curves: Effect of Stochastic Expression Models',
                 fontsize=12, fontweight='bold')
    ax5.legend(fontsize=10, loc='upper right', ncol=2)
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim(0, 600)
    
    # ========== Panel 9: Peak I_ND comparison ==========
    ax6 = fig.add_subplot(gs[2, 0])
    
    peak_I_NDs = []
    peak_distances = []
    for model_name, data in all_results.items():
        results = data['results']
        I_vals = [r['I_ND'] for r in results]
        dists = [r['distance'] for r in results]
        peak_idx = np.argmax(I_vals)
        peak_I_NDs.append(I_vals[peak_idx])
        peak_distances.append(dists[peak_idx])
    
    x_pos = np.arange(len(model_names))
    bars = ax6.bar(x_pos, peak_I_NDs, color=colors_models, alpha=0.7, edgecolor='black')
    
    # Add peak distance as text
    for i, (bar, dist) in enumerate(zip(bars, peak_distances)):
        height = bar.get_height()
        ax6.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{dist}μm', ha='center', va='bottom', fontsize=8)
    
    ax6.set_xticks(x_pos)
    ax6.set_xticklabels(model_names, rotation=45, ha='right', fontsize=9)
    ax6.set_ylabel('Peak I_ND', fontsize=10)
    ax6.set_title('(F) Peak I_ND Comparison\n(numbers = peak distance)',
                 fontsize=11, fontweight='bold')
    ax6.grid(True, alpha=0.3, axis='y')
    
    # ========== Panel 10: Temporal trajectories ==========
    ax7 = fig.add_subplot(gs[2, 1])
    
    # Show expression trajectories for one representative cell from each cluster
    timepoints = np.arange(n_timepoints)
    
    for ts_type, data in temporal_results.items():
        trajectory = data['trajectory']
        # Average over cells in high-expressing cluster (cluster 4)
        cluster_4_trajectory = trajectory[:, 60:80].mean(axis=1)
        ax7.plot(timepoints, cluster_4_trajectory, '-o', linewidth=2,
                markersize=6, label=f'{ts_type.capitalize()} timescale', alpha=0.7)
    
    ax7.set_xlabel('Time Point', fontsize=10)
    ax7.set_ylabel('Mean Expression (Cluster 4)', fontsize=10)
    ax7.set_title('(G) Temporal Expression Dynamics\nHigh-expressing cluster',
                 fontsize=11, fontweight='bold')
    ax7.legend(fontsize=9)
    ax7.grid(True, alpha=0.3)
    
    # ========== Panel 11: I_ND temporal variability ==========
    ax8 = fig.add_subplot(gs[2, 2])
    
    for ts_type, data in temporal_results.items():
        I_ND_over_time = data['I_ND_over_time']
        mean_I_ND = data['mean_I_ND']
        std_I_ND = data['std_I_ND']
        
        ax8.plot(timepoints, I_ND_over_time, '-s', linewidth=2,
                markersize=6, label=f'{ts_type.capitalize()}: μ={mean_I_ND:.3f}, σ={std_I_ND:.3f}',
                alpha=0.7)
    
    # Reference: static baseline
    static_I_ND = all_results['Static (baseline)']['results'][5]['I_ND']  # at 100 μm
    ax8.axhline(static_I_ND, color='black', linestyle='--', linewidth=2,
               label=f'Static baseline: {static_I_ND:.3f}', alpha=0.5)
    
    ax8.set_xlabel('Time Point', fontsize=10)
    ax8.set_ylabel('I_ND at d=100μm', fontsize=10)
    ax8.set_title('(H) Temporal I_ND Variability\nMeasured at optimal distance',
                 fontsize=11, fontweight='bold')
    ax8.legend(fontsize=8, loc='best')
    ax8.grid(True, alpha=0.3)
    
    # ========== Panel 12: Summary statistics ==========
    ax9 = fig.add_subplot(gs[2, 3])
    ax9.axis('off')
    
    summary_text = "STOCHASTIC EXPRESSION SUMMARY\n" + "="*40 + "\n\n"
    
    summary_text += "EXPRESSION MODELS:\n"
    for model_name, data in all_results.items():
        expr = data['sender_expr']
        peak_I = max([r['I_ND'] for r in data['results']])
        summary_text += f"  {model_name}:\n"
        summary_text += f"    Mean: {np.mean(expr):.2f}\n"
        summary_text += f"    CV: {np.std(expr)/np.mean(expr):.2f}\n"
        summary_text += f"    Peak I_ND: {peak_I:.3f}\n\n"
    
    summary_text += "TEMPORAL DYNAMICS:\n"
    for ts_type, data in temporal_results.items():
        summary_text += f"  {ts_type.capitalize()} timescale:\n"
        summary_text += f"    Mean I_ND: {data['mean_I_ND']:.3f}\n"
        summary_text += f"    Std I_ND: {data['std_I_ND']:.3f}\n"
        summary_text += f"    CV: {data['std_I_ND']/data['mean_I_ND']:.2f}\n\n"
    
    summary_text += "KEY FINDINGS:\n"
    summary_text += "  • Bursting ↑ variance → ↑ CV\n"
    summary_text += "  • Correlated noise → spatial\n"
    summary_text += "    structure preserved\n"
    summary_text += "  • Fast fluctuations → high\n"
    summary_text += "    temporal variability in I_ND\n"
    summary_text += "  • Slow drift → stable I_ND\n"
    
    ax9.text(0.05, 0.95, summary_text, transform=ax9.transAxes,
            fontsize=8, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
    
    fig.suptitle('Stochastic Sender Expression Effects on I_ND Detection',
                fontsize=14, fontweight='bold', y=0.98)
    
    plt.savefig('./IND_stochastic_expression.png',
                dpi=200, bbox_inches='tight')
    print("✓ Saved: IND_stochastic_expression.png")
    
    # ==========================================================================
    # Print Summary
    # ==========================================================================
    print("\n" + "=" * 70)
    print("STOCHASTIC EXPRESSION ANALYSIS COMPLETE")
    print("=" * 70)
    
    print("\nExpression Model Comparison:")
    for model_name, data in all_results.items():
        expr = data['sender_expr']
        results = data['results']
        peak_I = max([r['I_ND'] for r in results])
        peak_d = [r['distance'] for r in results][np.argmax([r['I_ND'] for r in results])]
        print(f"\n{model_name}:")
        print(f"  CV: {np.std(expr)/np.mean(expr):.3f}")
        print(f"  Peak I_ND: {peak_I:.3f} at {peak_d} μm")
    
    print("\nTemporal Dynamics:")
    for ts_type, data in temporal_results.items():
        print(f"\n{ts_type.capitalize()} timescale:")
        print(f"  Mean I_ND: {data['mean_I_ND']:.3f} ± {data['std_I_ND']:.3f}")
        print(f"  Temporal CV: {data['std_I_ND']/data['mean_I_ND']:.2f}")
    
    print("=" * 70)
    
    return all_results, temporal_results

if __name__ == "__main__":
    results_static, results_temporal = run_stochastic_analysis()
    plt.show()

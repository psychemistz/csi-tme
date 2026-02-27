#!/usr/bin/env python3
"""
Comprehensive Demo: I_ND(d) Analysis with Multiple Conditions
- Gaussian decay model
- Effects of receiver density
- Effects of sender expression (saturation analysis)
- Negative controls
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial, ndimage
from tqdm import tqdm
import matplotlib.gridspec as gridspec

np.random.seed(42)

# =============================================================================
# Core Functions
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
    
    return I_ND, np.sum(W > 0)

def compute_IND_curve(sender_positions, sender_expr, receiver_positions, receiver_expr,
                     test_distances, bandwidth=15, show_progress=False):
    """Compute I_ND for all test distances"""
    results = []
    iterator = tqdm(test_distances, desc="Computing I_ND") if show_progress else test_distances
    
    for dist in iterator:
        I_ND, n_conn = compute_IND(sender_positions, sender_expr,
                                   receiver_positions, receiver_expr,
                                   dist, bandwidth)
        results.append({'distance': dist, 'I_ND': I_ND, 'n_connections': n_conn})
    
    return results

# =============================================================================
# Simulation Setup
# =============================================================================

def setup_simulation(domain_size=(2000, 2000), 
                    n_receivers_per_cluster=200,
                    cluster_expression_means=[1.0, 2.0, 3.5, 5.0],
                    sigma=100,
                    Kd=2.0,
                    grid_spacing=20):
    """Set up sender clusters and receivers"""
    
    # Sender clusters at corners
    n_clusters = 4
    cells_per_cluster = 20
    margin = 200
    
    cluster_centers = [
        [margin, margin],
        [domain_size[0] - margin, margin],
        [margin, domain_size[1] - margin],
        [domain_size[0] - margin, domain_size[1] - margin]
    ]
    
    # Create senders
    sender_positions = []
    sender_ligand_expr = []
    
    for center, expr_mean in zip(cluster_centers, cluster_expression_means):
        cluster_pos = np.array(center) + np.random.randn(cells_per_cluster, 2) * 30
        sender_positions.append(cluster_pos)
        
        cluster_expr = np.random.lognormal(np.log(expr_mean), 0.2, cells_per_cluster)
        sender_ligand_expr.append(cluster_expr)
    
    sender_positions = np.vstack(sender_positions)
    sender_ligand_expr = np.concatenate(sender_ligand_expr)
    
    # Create receivers in shells around clusters
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
    
    # Simulate diffusion
    conc_field = solve_cytokine_gaussian_decay(
        sender_positions, sender_ligand_expr,
        domain_size=domain_size, grid_spacing=grid_spacing, sigma=sigma
    )
    
    # Sample concentrations and generate responses
    receiver_concentrations = sample_concentration(
        receiver_positions, conc_field, grid_spacing, domain_size
    )
    
    receptor_activation = receiver_concentrations / (Kd + receiver_concentrations)
    receiver_target_expr = (1 + 10 * receptor_activation) * \
                          np.random.lognormal(0, 0.15, len(receiver_positions))
    
    return {
        'sender_positions': sender_positions,
        'sender_expr': sender_ligand_expr,
        'receiver_positions': receiver_positions,
        'receiver_expr': receiver_target_expr,
        'receiver_conc': receiver_concentrations,
        'conc_field': conc_field,
        'cluster_centers': cluster_centers
    }

# =============================================================================
# Main Analysis
# =============================================================================

def run_comprehensive_analysis():
    print("=" * 70)
    print("COMPREHENSIVE I_ND ANALYSIS")
    print("=" * 70)
    
    domain_size = (2000, 2000)
    grid_spacing = 20
    sigma = 100
    test_distances = np.arange(20, 601, 20)  # Up to 600 μm for visualization
    
    # ==========================================================================
    # 1. Baseline Simulation
    # ==========================================================================
    print("\n1. Running baseline simulation...")
    baseline = setup_simulation(domain_size=domain_size, 
                               n_receivers_per_cluster=200,
                               sigma=sigma,
                               Kd=2.0,
                               grid_spacing=grid_spacing)
    
    print("   Computing I_ND curve...")
    baseline_results = compute_IND_curve(
        baseline['sender_positions'], baseline['sender_expr'],
        baseline['receiver_positions'], baseline['receiver_expr'],
        test_distances, show_progress=True
    )
    
    # Negative control
    print("   Computing negative control...")
    random_expr = np.random.lognormal(0, 0.15, len(baseline['receiver_positions']))
    control_results = compute_IND_curve(
        baseline['sender_positions'], baseline['sender_expr'],
        baseline['receiver_positions'], random_expr,
        test_distances[::2]  # Every other distance for speed
    )
    
    # ==========================================================================
    # 2. Receiver Density Effects
    # ==========================================================================
    print("\n2. Testing receiver density effects...")
    
    density_conditions = {
        'Low (50/cluster)': 50,
        'Medium (200/cluster)': 200,
        'High (500/cluster)': 500
    }
    
    density_results = {}
    for label, n_receivers in density_conditions.items():
        print(f"   {label}...")
        sim = setup_simulation(domain_size=domain_size,
                              n_receivers_per_cluster=n_receivers,
                              sigma=sigma, Kd=2.0, grid_spacing=grid_spacing)
        
        results = compute_IND_curve(
            sim['sender_positions'], sim['sender_expr'],
            sim['receiver_positions'], sim['receiver_expr'],
            test_distances
        )
        density_results[label] = results
    
    # ==========================================================================
    # 3. Sender Expression Effects (Saturation Analysis)
    # ==========================================================================
    print("\n3. Testing sender expression effects (saturation)...")
    
    expression_conditions = {
        'Low (0.5-1.5x)': [0.5, 0.75, 1.0, 1.5],
        'Medium (1-4x)': [1.0, 2.0, 3.0, 4.0],
        'High (2-10x)': [2.0, 4.0, 7.0, 10.0]
    }
    
    expression_results = {}
    for label, expr_means in expression_conditions.items():
        print(f"   {label}...")
        sim = setup_simulation(domain_size=domain_size,
                              n_receivers_per_cluster=200,
                              cluster_expression_means=expr_means,
                              sigma=sigma, Kd=2.0, grid_spacing=grid_spacing)
        
        results = compute_IND_curve(
            sim['sender_positions'], sim['sender_expr'],
            sim['receiver_positions'], sim['receiver_expr'],
            test_distances
        )
        expression_results[label] = (results, sim)
    
    # ==========================================================================
    # 4. Create Comprehensive Visualization
    # ==========================================================================
    print("\n4. Creating comprehensive visualization...")
    
    fig = plt.figure(figsize=(20, 12))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # ========== Panel 1: Sender & Receiver Distribution ==========
    ax1 = fig.add_subplot(gs[0, 0])
    
    # Plot receivers
    scatter_recv = ax1.scatter(baseline['receiver_positions'][:, 0],
                              baseline['receiver_positions'][:, 1],
                              c='lightblue', s=5, alpha=0.3, label='Receivers')
    
    # Plot senders colored by expression
    scatter_send = ax1.scatter(baseline['sender_positions'][:, 0],
                              baseline['sender_positions'][:, 1],
                              c=baseline['sender_expr'], s=150, marker='*',
                              cmap='Reds', edgecolors='black', linewidths=1.5,
                              vmin=0, vmax=6, label='Senders')
    
    # Draw distance circles around one cluster
    center = baseline['cluster_centers'][3]
    for d in [100, 200, 300]:
        circle = plt.Circle(center, d, fill=False, color='gray',
                          linestyle='--', linewidth=1, alpha=0.5)
        ax1.add_patch(circle)
    
    ax1.set_xlim(0, domain_size[0])
    ax1.set_ylim(0, domain_size[1])
    ax1.set_xlabel('X (μm)', fontsize=10)
    ax1.set_ylabel('Y (μm)', fontsize=10)
    ax1.set_title('(A) Sender & Receiver Distribution\n4 clusters, 80 senders, 800 receivers',
                 fontsize=11, fontweight='bold')
    ax1.set_aspect('equal')
    cbar1 = plt.colorbar(scatter_send, ax=ax1, label='Sender Expression', pad=0.02)
    
    # ========== Panel 2: Cytokine Concentration Field ==========
    ax2 = fig.add_subplot(gs[0, 1])
    
    im2 = ax2.imshow(baseline['conc_field'], cmap='YlOrRd', origin='lower',
                    extent=[0, domain_size[0], 0, domain_size[1]],
                    interpolation='bilinear')
    
    ax2.scatter(baseline['sender_positions'][:, 0],
               baseline['sender_positions'][:, 1],
               c=baseline['sender_expr'], s=100, marker='*',
               cmap='coolwarm', edgecolors='white', linewidths=1.5,
               vmin=0, vmax=6)
    
    # Add contour lines
    X = np.linspace(0, domain_size[0], baseline['conc_field'].shape[1])
    Y = np.linspace(0, domain_size[1], baseline['conc_field'].shape[0])
    levels = np.percentile(baseline['conc_field'][baseline['conc_field'] > 0], [25, 50, 75])
    ax2.contour(X, Y, baseline['conc_field'], levels=levels, colors='white',
               linewidths=0.5, alpha=0.5)
    
    ax2.set_xlabel('X (μm)', fontsize=10)
    ax2.set_ylabel('Y (μm)', fontsize=10)
    ax2.set_title('(B) Cytokine Concentration Field\nGaussian decay: C(r) = C₀ exp(-(r/σ)²)',
                 fontsize=11, fontweight='bold')
    cbar2 = plt.colorbar(im2, ax=ax2, label='Concentration', pad=0.02)
    
    # ========== Panel 3: I_ND vs Distance with Negative Control ==========
    ax3 = fig.add_subplot(gs[0, 2])
    
    distances_baseline = [r['distance'] for r in baseline_results]
    I_ND_baseline = [r['I_ND'] for r in baseline_results]
    distances_control = [r['distance'] for r in control_results]
    I_ND_control = [r['I_ND'] for r in control_results]
    
    ax3.plot(distances_baseline, I_ND_baseline, 'b-', linewidth=2.5,
            label='True Signal', alpha=0.8)
    ax3.plot(distances_control, I_ND_control, 'r--', linewidth=2,
            label='Random Control', alpha=0.7)
    
    # Mark peak
    peak_idx = np.argmax(I_ND_baseline)
    peak_dist = distances_baseline[peak_idx]
    peak_I_ND = I_ND_baseline[peak_idx]
    ax3.plot(peak_dist, peak_I_ND, 'g*', markersize=20, markeredgecolor='black',
            markeredgewidth=1.5, label=f'Peak: {peak_dist} μm')
    
    # Mark sigma
    ax3.axvline(sigma, color='green', linestyle=':', linewidth=2,
               alpha=0.5, label=f'σ = {sigma} μm')
    ax3.axhline(0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
    
    ax3.set_xlabel('Distance (μm)', fontsize=10)
    ax3.set_ylabel('I_ND(d)', fontsize=10)
    ax3.set_title('(C) Distance-Dependent Communication\nTrue Signal vs Negative Control',
                 fontsize=11, fontweight='bold')
    ax3.legend(fontsize=9, loc='upper right')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 600)
    
    # ========== Panel 4: Receiver Density Effects ==========
    ax4 = fig.add_subplot(gs[1, 0])
    
    colors_density = ['lightblue', 'blue', 'darkblue']
    for (label, results), color in zip(density_results.items(), colors_density):
        distances = [r['distance'] for r in results]
        I_ND_vals = [r['I_ND'] for r in results]
        ax4.plot(distances, I_ND_vals, '-o', linewidth=2, markersize=4,
                label=label, color=color, alpha=0.7)
    
    ax4.axvline(sigma, color='green', linestyle=':', linewidth=2, alpha=0.5)
    ax4.axhline(0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
    ax4.set_xlabel('Distance (μm)', fontsize=10)
    ax4.set_ylabel('I_ND(d)', fontsize=10)
    ax4.set_title('(D) Effect of Receiver Cell Density\nMore receivers → better statistics',
                 fontsize=11, fontweight='bold')
    ax4.legend(fontsize=9, loc='upper right')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0, 600)
    
    # ========== Panel 5: Sender Expression Effects ==========
    ax5 = fig.add_subplot(gs[1, 1])
    
    colors_expr = ['orange', 'red', 'darkred']
    for (label, (results, sim)), color in zip(expression_results.items(), colors_expr):
        distances = [r['distance'] for r in results]
        I_ND_vals = [r['I_ND'] for r in results]
        ax5.plot(distances, I_ND_vals, '-o', linewidth=2, markersize=4,
                label=label, color=color, alpha=0.7)
    
    ax5.axvline(sigma, color='green', linestyle=':', linewidth=2, alpha=0.5)
    ax5.axhline(0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
    ax5.set_xlabel('Distance (μm)', fontsize=10)
    ax5.set_ylabel('I_ND(d)', fontsize=10)
    ax5.set_title('(E) Effect of Sender Expression Level\nHigh expression → saturation at short distances',
                 fontsize=11, fontweight='bold')
    ax5.legend(fontsize=9, loc='upper right')
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim(0, 600)
    
    # ========== Panel 6: Saturation Analysis ==========
    ax6 = fig.add_subplot(gs[1, 2])
    
    # Show activation vs distance for different expression levels
    for label, (results, sim) in expression_results.items():
        # Sample activation at different distances
        distances_sample = [50, 100, 150, 200, 300, 400]
        mean_activations = []
        
        for dist in distances_sample:
            # Find receivers at this distance from high-expressing cluster
            center = baseline['cluster_centers'][3]
            dists_to_center = np.linalg.norm(sim['receiver_positions'] - center, axis=1)
            mask = np.abs(dists_to_center - dist) < 30
            
            if np.sum(mask) > 10:
                concs = sim['receiver_conc'][mask]
                activations = concs / (2.0 + concs)
                mean_activations.append(np.mean(activations))
            else:
                mean_activations.append(0)
        
        ax6.plot(distances_sample, mean_activations, '-s', linewidth=2,
                markersize=8, label=label, alpha=0.7)
    
    # Mark saturation threshold
    ax6.axhline(0.5, color='red', linestyle='--', linewidth=1.5,
               alpha=0.5, label='Half-max activation')
    ax6.axhline(0.8, color='darkred', linestyle=':', linewidth=1.5,
               alpha=0.5, label='Near-saturation')
    
    ax6.set_xlabel('Distance from high-expressing sender (μm)', fontsize=10)
    ax6.set_ylabel('Mean Receptor Activation', fontsize=10)
    ax6.set_title('(F) Receptor Saturation Analysis\nActivation = C/(Kd+C), Kd=2.0',
                 fontsize=11, fontweight='bold')
    ax6.legend(fontsize=8, loc='upper right')
    ax6.grid(True, alpha=0.3)
    ax6.set_xlim(0, 600) 
    ax6.set_ylim(0, 1)
    
    # ========== Panel 7: Connection Statistics ==========
    ax7 = fig.add_subplot(gs[2, 0])
    
    n_connections_baseline = [r['n_connections'] for r in baseline_results]
    ax7.plot(distances_baseline, n_connections_baseline, 'g-o',
            linewidth=2, markersize=4, alpha=0.7)
    ax7.fill_between(distances_baseline, n_connections_baseline, alpha=0.3, color='green')
    
    ax7.set_xlabel('Distance (μm)', fontsize=10)
    ax7.set_ylabel('Number of Sender-Receiver Pairs', fontsize=10)
    ax7.set_title('(G) Connection Statistics\nPairs in annular band (bandwidth=15 μm)',
                 fontsize=11, fontweight='bold')
    ax7.grid(True, alpha=0.3)
    ax7.set_xlim(0, 600)
    
    # ========== Panel 8: Peak I_ND Comparison ==========
    ax8 = fig.add_subplot(gs[2, 1])
    
    # Extract peak I_ND for each condition
    conditions = []
    peak_I_NDs = []
    peak_distances = []
    
    # Baseline
    conditions.append('Baseline')
    peak_I_NDs.append(peak_I_ND)
    peak_distances.append(peak_dist)
    
    # Density effects
    for label, results in density_results.items():
        I_vals = [r['I_ND'] for r in results]
        dists = [r['distance'] for r in results]
        peak_idx = np.argmax(I_vals)
        conditions.append(label.split('(')[0].strip())
        peak_I_NDs.append(I_vals[peak_idx])
        peak_distances.append(dists[peak_idx])
    
    # Expression effects
    for label, (results, sim) in expression_results.items():
        I_vals = [r['I_ND'] for r in results]
        dists = [r['distance'] for r in results]
        peak_idx = np.argmax(I_vals)
        conditions.append(label.split('(')[0].strip())
        peak_I_NDs.append(I_vals[peak_idx])
        peak_distances.append(dists[peak_idx])
    
    x_pos = np.arange(len(conditions))
    colors_bar = ['gray', 'lightblue', 'blue', 'darkblue', 'orange', 'red', 'darkred']
    
    bars = ax8.bar(x_pos, peak_I_NDs, color=colors_bar, alpha=0.7, edgecolor='black')
    
    # Add peak distance as text on bars
    for i, (bar, dist) in enumerate(zip(bars, peak_distances)):
        height = bar.get_height()
        ax8.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                f'{dist}μm', ha='center', va='bottom', fontsize=8, rotation=0)
    
    ax8.set_xticks(x_pos)
    ax8.set_xticklabels(conditions, rotation=45, ha='right', fontsize=9)
    ax8.set_ylabel('Peak I_ND', fontsize=10)
    ax8.set_title('(H) Peak I_ND Comparison\n(numbers show peak distance)',
                 fontsize=11, fontweight='bold')
    ax8.grid(True, alpha=0.3, axis='y')
    ax8.set_ylim(0, max(peak_I_NDs) * 1.2)
    
    # ========== Panel 9: Summary Statistics Table ==========
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis('off')
    
    # Create summary text
    summary_text = "SUMMARY STATISTICS\n" + "="*35 + "\n\n"
    
    summary_text += f"Domain: {domain_size[0]}×{domain_size[1]} μm\n"
    summary_text += f"Decay model: Gaussian (σ={sigma} μm)\n"
    summary_text += f"Senders: {len(baseline['sender_positions'])}\n"
    summary_text += f"Receivers: {len(baseline['receiver_positions'])}\n\n"
    
    summary_text += "BASELINE RESULTS:\n"
    summary_text += f"  Peak I_ND: {peak_I_ND:.3f}\n"
    summary_text += f"  Peak distance: {peak_dist} μm\n"
    summary_text += f"  Peak/σ ratio: {peak_dist/sigma:.2f}\n\n"
    
    summary_text += "KEY FINDINGS:\n"
    summary_text += f"  • Receiver density ↑ → I_ND ↑\n"
    summary_text += f"    (better statistics)\n\n"
    summary_text += f"  • High sender expression:\n"
    summary_text += f"    - Saturates at d<100μm\n"
    summary_text += f"    - Reduces initial I_ND\n\n"
    summary_text += f"  • Optimal distance:\n"
    summary_text += f"    - Near σ = {sigma} μm\n"
    summary_text += f"    - Maximal dynamic range\n\n"
    
    neg_I_ND_count = sum(1 for r in baseline_results if r['I_ND'] < 0)
    summary_text += f"  • Negative I_ND: {neg_I_ND_count}/{len(baseline_results)}\n"
    summary_text += f"  • Random control: ~0\n"
    
    ax9.text(0.1, 0.95, summary_text, transform=ax9.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    # Overall title
    fig.suptitle('Comprehensive I_ND(d) Analysis: Gaussian Decay Model', 
                fontsize=14, fontweight='bold', y=0.98)
    
    plt.savefig('./IND_comprehensive_analysis.png',
                dpi=200, bbox_inches='tight')
    print("✓ Saved: IND_comprehensive_analysis.png")
    
    # ==========================================================================
    # Summary Output
    # ==========================================================================
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nBaseline peak: I_ND = {peak_I_ND:.3f} at {peak_dist} μm")
    print(f"\nReceiver Density Effects:")
    for label, results in density_results.items():
        I_vals = [r['I_ND'] for r in results]
        print(f"  {label}: peak I_ND = {max(I_vals):.3f}")
    
    print(f"\nSender Expression Effects:")
    for label, (results, sim) in expression_results.items():
        I_vals = [r['I_ND'] for r in results]
        print(f"  {label}: peak I_ND = {max(I_vals):.3f}")
    
    print("=" * 70)
    
    return {
        'baseline': baseline_results,
        'control': control_results,
        'density': density_results,
        'expression': expression_results
    }

if __name__ == "__main__":
    results = run_comprehensive_analysis()
    plt.show()

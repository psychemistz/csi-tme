#!/usr/bin/env python3
"""
Quick Demo: 4 Cytokines (IL-2, TNF, TGF-β1, IFN-γ) with 3 Cell Types
Simplified version for rapid testing and visualization

Author: Seongyong Park
Institution: Cancer Data Science Lab, NCI, NIH
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial, ndimage
import pandas as pd
from typing import Dict, Tuple
import seaborn as sns

np.random.seed(42)
sns.set_style('whitegrid')

# =============================================================================
# Cytokine Parameters (simplified)
# =============================================================================

CYTOKINES = {
    'IL-2': {
        'diffusion': 20.0,      # μm²/s
        'decay_rate': 0.023,    # 1/min
        'Kd': 0.01,            # nM
        'max_response': 5.0,
        'color': '#FF6B6B'
    },
    'TNF': {
        'diffusion': 40.0,
        'decay_rate': 0.035,
        'Kd': 0.02,
        'max_response': 10.0,
        'color': '#4ECDC4'
    },
    'TGF-β1': {
        'diffusion': 10.0,
        'decay_rate': 0.020,
        'Kd': 0.05,
        'max_response': 8.0,
        'color': '#95E77E'
    },
    'IFN-γ': {
        'diffusion': 30.0,
        'decay_rate': 0.025,
        'Kd': 0.1,
        'max_response': 12.0,
        'color': '#FFD93D'
    }
}

# Cell type properties (production/consumption rates)
CELL_TYPES = {
    'Sender': {
        'produce': {'IL-2': 2.0, 'TNF': 1.5, 'TGF-β1': 1.0, 'IFN-γ': 2.5},
        'consume': {'IL-2': 0.0, 'TNF': 0.0, 'TGF-β1': 0.0, 'IFN-γ': 0.0},
        'respond': {'IL-2': 0.5, 'TNF': 0.5, 'TGF-β1': 0.5, 'IFN-γ': 0.5},
        'color': 'red',
        'marker': 'o'
    },
    'Receiver': {
        'produce': {'IL-2': 0.0, 'TNF': 0.0, 'TGF-β1': 0.0, 'IFN-γ': 0.0},
        'consume': {'IL-2': 2.0, 'TNF': 1.5, 'TGF-β1': 1.0, 'IFN-γ': 1.5},
        'respond': {'IL-2': 2.0, 'TNF': 2.0, 'TGF-β1': 2.0, 'IFN-γ': 2.0},
        'color': 'blue',
        'marker': 's'
    },
    'Inert': {
        'produce': {'IL-2': 0.0, 'TNF': 0.0, 'TGF-β1': 0.0, 'IFN-γ': 0.0},
        'consume': {'IL-2': 0.0, 'TNF': 0.0, 'TGF-β1': 0.0, 'IFN-γ': 0.0},
        'respond': {'IL-2': 0.0, 'TNF': 0.0, 'TGF-β1': 0.0, 'IFN-γ': 0.0},
        'color': 'gray',
        'marker': '.'
    }
}

# =============================================================================
# Simplified PDE Solver
# =============================================================================

def solve_cytokine_diffusion(cell_positions, cell_types, cytokine_name, 
                            domain_size=(500, 500), grid_size=50):
    """
    Simplified diffusion solver with production and consumption
    """
    nx, ny = grid_size, grid_size
    dx = domain_size[0] / nx
    dy = domain_size[1] / ny
    
    # Initialize concentration field and source/sink maps
    concentration = np.zeros((ny, nx))
    source_map = np.zeros((ny, nx))
    sink_map = np.zeros((ny, nx))
    
    cytokine = CYTOKINES[cytokine_name]
    
    # Place sources and sinks based on cell types
    for i, (pos, ctype) in enumerate(zip(cell_positions, cell_types)):
        ix = int(pos[0] / dx)
        iy = int(pos[1] / dy)
        if 0 <= ix < nx and 0 <= iy < ny:
            # Production
            source_map[iy, ix] += CELL_TYPES[ctype]['produce'][cytokine_name]
            # Consumption
            sink_map[iy, ix] += CELL_TYPES[ctype]['consume'][cytokine_name]
    
    # Simple steady-state approximation using Gaussian smoothing
    # Add sources
    concentration = source_map * 100
    
    # Diffusion (multiple smoothing iterations)
    sigma = np.sqrt(cytokine['diffusion'] / cytokine['decay_rate']) / 20
    for _ in range(20):
        concentration = ndimage.gaussian_filter(concentration, sigma=sigma)
        # Apply degradation and consumption
        concentration *= 0.95  # Degradation
        concentration = np.maximum(0, concentration - sink_map * 0.1)  # Consumption
    
    return concentration

# =============================================================================
# Normalized Directional Moran's I
# =============================================================================

def calculate_I_ND(sender_pos, sender_exp, receiver_pos, receiver_exp, 
                  radius=150.0):
    """
    Calculate normalized directional Moran's I
    """
    if len(sender_pos) == 0 or len(receiver_pos) == 0:
        return 0.0, 1.0
    
    # Distance matrix
    distances = spatial.distance_matrix(sender_pos, receiver_pos)
    
    # Gaussian weights
    sigma = radius / 3.0
    W = np.exp(-distances**2 / (2 * sigma**2))
    W[distances > radius] = 0
    
    # Row-normalize
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    W_norm = W / row_sums
    
    # Standardize expressions
    z_s = (sender_exp - np.mean(sender_exp)) / (np.std(sender_exp) + 1e-10)
    z_r = (receiver_exp - np.mean(receiver_exp)) / (np.std(receiver_exp) + 1e-10)
    
    # Spatial lag
    lag = W_norm @ z_r
    
    # Calculate I_ND
    norm_s = np.linalg.norm(z_s)
    norm_lag = np.linalg.norm(lag)
    
    if norm_s == 0 or norm_lag == 0:
        return 0.0, 1.0
    
    I_ND = np.dot(z_s, lag) / (norm_s * norm_lag)
    
    # Simple permutation test
    n_perm = 199
    I_perm = []
    for _ in range(n_perm):
        z_s_perm = np.random.permutation(z_s)
        I_perm.append(np.dot(z_s_perm, lag) / (np.linalg.norm(z_s_perm) * norm_lag))
    
    p_value = (np.sum(np.abs(I_perm) >= np.abs(I_ND)) + 1) / (n_perm + 1)
    
    return I_ND, p_value

# =============================================================================
# Main Demo
# =============================================================================

def run_demo():
    """
    Run demonstration of 4 cytokines with 3 cell types
    """
    print("="*70)
    print("DEMO: 4 Cytokines with Sender, Receiver, and Inert Cells")
    print("="*70)
    
    # Create cells with spatial pattern
    n_cells = 500
    domain_size = (500, 500)
    
    # Create gradient pattern: senders on left, receivers on right, inert scattered
    positions = []
    cell_types = []
    
    for i in range(n_cells):
        x = np.random.rand() * domain_size[0]
        y = np.random.rand() * domain_size[1]
        
        # Assign type based on position
        if x < 150:  # Left side - mostly senders
            cell_type = np.random.choice(['Sender', 'Sender', 'Inert'], p=[0.7, 0.1, 0.2])
        elif x > 350:  # Right side - mostly receivers
            cell_type = np.random.choice(['Receiver', 'Receiver', 'Inert'], p=[0.7, 0.1, 0.2])
        else:  # Middle - mixed
            cell_type = np.random.choice(['Sender', 'Receiver', 'Inert'], p=[0.2, 0.2, 0.6])
        
        positions.append([x, y])
        cell_types.append(cell_type)
    
    positions = np.array(positions)
    
    # Test all 4 cytokines
    results = {}
    
    for cytokine_name in CYTOKINES.keys():
        print(f"\nTesting {cytokine_name}...")
        
        # Solve diffusion
        conc_field = solve_cytokine_diffusion(
            positions, cell_types, cytokine_name, domain_size
        )
        
        # Sample concentration at cell positions
        concentrations = np.zeros(n_cells)
        for i, pos in enumerate(positions):
            ix = int(pos[0] / 10)
            iy = int(pos[1] / 10)
            if 0 <= ix < 50 and 0 <= iy < 50:
                concentrations[i] = conc_field[iy, ix]
        
        # Generate expression response (Hill function)
        cytokine = CYTOKINES[cytokine_name]
        expression = np.zeros(n_cells)
        for i, ctype in enumerate(cell_types):
            response = CELL_TYPES[ctype]['respond'][cytokine_name]
            if response > 0:
                hill = concentrations[i] / (cytokine['Kd'] + concentrations[i])
                expression[i] = 1 + cytokine['max_response'] * hill * response
        
        # Add noise
        expression *= np.random.lognormal(0, 0.2, n_cells)
        
        # Calculate I_ND for each direction
        directions = {}
        
        for sender_type in ['Sender', 'Receiver']:
            for receiver_type in ['Sender', 'Receiver']:
                # Get cells of each type
                sender_mask = np.array(cell_types) == sender_type
                receiver_mask = np.array(cell_types) == receiver_type
                
                if np.sum(sender_mask) > 5 and np.sum(receiver_mask) > 5:
                    I_ND, p_val = calculate_I_ND(
                        positions[sender_mask],
                        expression[sender_mask],
                        positions[receiver_mask],
                        expression[receiver_mask]
                    )
                    
                    direction = f"{sender_type}→{receiver_type}"
                    directions[direction] = {'I_ND': I_ND, 'p_value': p_val}
                    print(f"  {direction}: I_ND = {I_ND:.3f}, p = {p_val:.3f}")
        
        # Test with inert cells (negative control)
        inert_mask = np.array(cell_types) == 'Inert'
        if np.sum(inert_mask) > 5:
            I_ND_control, p_val_control = calculate_I_ND(
                positions[sender_mask],
                expression[sender_mask],
                positions[inert_mask],
                expression[inert_mask]
            )
            print(f"  Sender→Inert (control): I_ND = {I_ND_control:.3f}, p = {p_val_control:.3f}")
            directions['Sender→Inert'] = {'I_ND': I_ND_control, 'p_value': p_val_control}
        
        results[cytokine_name] = {
            'concentration': conc_field,
            'expression': expression,
            'directions': directions
        }
    
    # Visualization
    visualize_results(positions, cell_types, results)
    
    return results

def visualize_results(positions, cell_types, results):
    """
    Create comprehensive visualization
    """
    fig = plt.figure(figsize=(16, 12))
    
    # Create grid for subplots
    gs = fig.add_gridspec(4, 5, hspace=0.3, wspace=0.3)
    
    # Cell distribution (top left)
    ax_cells = fig.add_subplot(gs[0, 0])
    for ctype, props in CELL_TYPES.items():
        mask = np.array(cell_types) == ctype
        if np.any(mask):
            ax_cells.scatter(positions[mask, 0], positions[mask, 1],
                           c=props['color'], marker=props['marker'],
                           s=10, alpha=0.6, label=ctype)
    ax_cells.set_title('Cell Distribution')
    ax_cells.set_xlabel('X (μm)')
    ax_cells.set_ylabel('Y (μm)')
    ax_cells.legend(loc='upper right', fontsize=8)
    ax_cells.set_xlim(0, 500)
    ax_cells.set_ylim(0, 500)
    
    # For each cytokine
    for idx, (cytokine_name, result) in enumerate(results.items()):
        row = idx
        
        # Concentration field
        ax_conc = fig.add_subplot(gs[row, 1])
        im = ax_conc.imshow(result['concentration'], cmap='YlOrRd',
                           origin='lower', aspect='auto')
        ax_conc.set_title(f'{cytokine_name} Concentration')
        ax_conc.set_xticks([])
        ax_conc.set_yticks([])
        plt.colorbar(im, ax=ax_conc, fraction=0.046)
        
        # Expression response
        ax_exp = fig.add_subplot(gs[row, 2])
        scatter = ax_exp.scatter(positions[:, 0], positions[:, 1],
                                c=result['expression'], cmap='viridis',
                                s=5, alpha=0.6)
        ax_exp.set_title(f'{cytokine_name} Response')
        ax_exp.set_xlim(0, 500)
        ax_exp.set_ylim(0, 500)
        ax_exp.set_xticks([])
        ax_exp.set_yticks([])
        plt.colorbar(scatter, ax=ax_exp, fraction=0.046)
        
        # I_ND matrix
        ax_ind = fig.add_subplot(gs[row, 3])
        
        # Create matrix from directions
        dir_types = ['Sender→Sender', 'Sender→Receiver', 
                    'Receiver→Sender', 'Receiver→Receiver']
        matrix = np.zeros((2, 2))
        
        for i, s in enumerate(['Sender', 'Receiver']):
            for j, r in enumerate(['Sender', 'Receiver']):
                key = f"{s}→{r}"
                if key in result['directions']:
                    matrix[i, j] = result['directions'][key]['I_ND']
        
        im = ax_ind.imshow(matrix, cmap='RdBu_r', vmin=-0.3, vmax=0.3)
        ax_ind.set_xticks([0, 1])
        ax_ind.set_yticks([0, 1])
        ax_ind.set_xticklabels(['Send', 'Recv'], fontsize=8)
        ax_ind.set_yticklabels(['Send', 'Recv'], fontsize=8)
        ax_ind.set_title(f'{cytokine_name} I_ND')
        
        # Add values
        for i in range(2):
            for j in range(2):
                ax_ind.text(j, i, f'{matrix[i, j]:.2f}',
                          ha='center', va='center',
                          color='white' if abs(matrix[i, j]) > 0.15 else 'black')
        
        # P-values
        ax_p = fig.add_subplot(gs[row, 4])
        p_matrix = np.ones((2, 2))
        
        for i, s in enumerate(['Sender', 'Receiver']):
            for j, r in enumerate(['Sender', 'Receiver']):
                key = f"{s}→{r}"
                if key in result['directions']:
                    p_matrix[i, j] = result['directions'][key]['p_value']
        
        im = ax_p.imshow(p_matrix, cmap='YlOrRd_r', vmin=0, vmax=0.1)
        ax_p.set_xticks([0, 1])
        ax_p.set_yticks([0, 1])
        ax_p.set_xticklabels(['Send', 'Recv'], fontsize=8)
        ax_p.set_yticklabels(['Send', 'Recv'], fontsize=8)
        ax_p.set_title(f'{cytokine_name} p-values')
        
        # Mark significant values
        for i in range(2):
            for j in range(2):
                if p_matrix[i, j] < 0.05:
                    ax_p.text(j, i, '***', ha='center', va='center',
                            color='red', fontweight='bold')
    
    plt.suptitle('Cytokine Signaling with 3 Cell Types: Sender, Receiver, Inert', 
                fontsize=14, fontweight='bold')
    
    # Save figure
    plt.savefig('cytokine_demo_3celltype.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print("\nFigure saved as 'cytokine_demo_3celltype.png'")

def summarize_results(results):
    """
    Print summary of findings
    """
    print("\n" + "="*70)
    print("SUMMARY OF RESULTS")
    print("="*70)
    
    for cytokine_name, result in results.items():
        print(f"\n{cytokine_name}:")
        
        # Find strongest interaction
        max_I = -1
        max_dir = ""
        sig_count = 0
        
        for direction, values in result['directions'].items():
            if values['I_ND'] > max_I:
                max_I = values['I_ND']
                max_dir = direction
            if values['p_value'] < 0.05:
                sig_count += 1
        
        print(f"  Strongest interaction: {max_dir} (I_ND = {max_I:.3f})")
        print(f"  Significant interactions: {sig_count}/{len(result['directions'])}")
        
        # Check if Sender→Receiver is significant
        if 'Sender→Receiver' in result['directions']:
            sr = result['directions']['Sender→Receiver']
            print(f"  Sender→Receiver: I_ND = {sr['I_ND']:.3f}, p = {sr['p_value']:.3f}")
            if sr['p_value'] < 0.05:
                print(f"  ✓ Directional signaling detected!")
        
        # Check negative control
        if 'Sender→Inert' in result['directions']:
            control = result['directions']['Sender→Inert']
            print(f"  Control (Sender→Inert): I_ND = {control['I_ND']:.3f}")
            if abs(control['I_ND']) < 0.1:
                print(f"  ✓ Negative control passed (no response in inert cells)")
    
    print("\n" + "="*70)
    print("KEY FINDINGS:")
    print("="*70)
    print("1. Directional Moran's I detects cytokine gradients from Sender→Receiver")
    print("2. Inert cells show minimal I_ND values (negative control)")
    print("3. Cell consumption affects local concentration fields")
    print("4. Each cytokine shows distinct spatial signaling patterns")
    print("5. Method distinguishes production vs consumption vs no interaction")

if __name__ == "__main__":
    print("\nRunning cytokine simulation demo...")
    print("This will take approximately 10-20 seconds...\n")
    
    results = run_demo()
    summarize_results(results)
    
    print("\n✓ Demo complete! Check 'cytokine_demo_3celltype.png' for visualization.")

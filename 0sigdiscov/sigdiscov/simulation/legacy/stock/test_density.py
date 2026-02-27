#!/usr/bin/env python3
"""
Gaussian Decay Spatial Transcriptomics Simulation
With configurable receiver vs inert density ratios

Author: Seongyong Park
Institution: Cancer Data Science Lab, NCI, NIH
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial, stats
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from typing import Tuple, Optional, Dict, List
import warnings
warnings.filterwarnings('ignore')

sns.set_style('whitegrid')
np.random.seed(42)

# =============================================================================
# Configuration Class
# =============================================================================

class SpatialConfig:
    """Configuration for spatial analysis"""
    
    def __init__(self, sigma_relationship: str = 'standard', sigma_cutoff_multiplier: float = 3.0):
        self.sigma_relationship = sigma_relationship
        self.sigma_cutoff_multiplier = sigma_cutoff_multiplier
        self.MIN_CELLS_THRESHOLD = 10
        self.WEIGHT_THRESHOLD = 1e-6
        
    def get_sigma(self, radius: float) -> float:
        if self.sigma_relationship == 'tight':
            return radius / 4
        elif self.sigma_relationship == 'standard':
            return radius / 3
        else:
            return radius / 2
    
    def get_max_signaling_distance(self, lambda_val: float) -> float:
        """Get theoretical cutoff distance based on sigma"""
        sigma = self.get_sigma(lambda_val)
        return sigma * self.sigma_cutoff_multiplier

# =============================================================================
# Tissue Generation - WITH DENSITY CONTROL
# =============================================================================

class TissueGenerator:
    @staticmethod
    def create_tissue_with_sender_cluster(
        n_senders: int = 500,
        n_non_senders: int = 3500,  # Total receivers + inert
        receiver_fraction: float = 0.5,  # Fraction of non-senders that are receivers
        domain_size: float = 5000,
        cluster_radius: float = 600
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """
        Create tissue with configurable receiver/inert density
        
        Parameters:
        -----------
        receiver_fraction : float
            Fraction of non-sender cells that are receivers
            0.9 = 90% receivers, 10% inert
            0.5 = 50% receivers, 50% inert
            0.1 = 10% receivers, 90% inert
        """
        positions = []
        cell_types = []
        center = np.array([domain_size/2, domain_size/2])
        
        # Clustered senders
        for _ in range(n_senders):
            r = np.random.exponential(cluster_radius/3)
            r = min(r, cluster_radius)
            theta = np.random.uniform(0, 2*np.pi)
            
            x = center[0] + r * np.cos(theta)
            y = center[1] + r * np.sin(theta)
            x = np.clip(x, 0, domain_size)
            y = np.clip(y, 0, domain_size)
            
            positions.append([x, y])
            cell_types.append('Sender')
        
        # Calculate numbers
        n_receivers = int(n_non_senders * receiver_fraction)
        n_inert = n_non_senders - n_receivers
        
        # Uniformly distributed receivers
        for _ in range(n_receivers):
            positions.append([
                np.random.uniform(0, domain_size),
                np.random.uniform(0, domain_size)
            ])
            cell_types.append('Receiver')
        
        # Uniformly distributed inert cells
        for _ in range(n_inert):
            positions.append([
                np.random.uniform(0, domain_size),
                np.random.uniform(0, domain_size)
            ])
            cell_types.append('Inert')
        
        return np.array(positions), cell_types, center

# =============================================================================
# Expression Generation
# =============================================================================

class ExpressionGenerator:
    """Generate expression with ~20% zero counts overall"""
    
    def __init__(self, expression_params: Optional[Dict] = None):
        self.params = expression_params or {
            'sender_expressing_fraction': 0.90,
            'receiver_baseline_fraction': 0.85,
            'technical_dropout': 0.05,
            'sender_high_mean': 30,
            'sender_high_dispersion': 0.5,
            'receiver_baseline_mean': 5,
            'receiver_induced_mean': 25,
            'inert_expression_rate': 0.85
        }
    
    @staticmethod
    def generate_zinb_counts(mean: float, size: float, zero_prob: float, n: int = 1) -> np.ndarray:
        is_zero = np.random.rand(n) < zero_prob
        counts = np.zeros(n)
        non_zero_mask = ~is_zero
        
        if np.any(non_zero_mask):
            n_nonzero = np.sum(non_zero_mask)
            p = size / (size + mean)
            counts[non_zero_mask] = np.random.negative_binomial(size, p, n_nonzero)
        
        return counts if n > 1 else counts[0]
    
    def generate_expression_with_gaussian_decay(
        self,
        positions: np.ndarray,
        cell_types: List[str],
        lambda_val: float,
        config: SpatialConfig
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        n_cells = len(cell_types)
        expression = np.zeros(n_cells)
        is_expressing = np.zeros(n_cells, dtype=bool)
        
        sender_mask = np.array([ct == 'Sender' for ct in cell_types])
        receiver_mask = np.array([ct == 'Receiver' for ct in cell_types])
        inert_mask = np.array([ct == 'Inert' for ct in cell_types])
        
        sender_indices = np.where(sender_mask)[0]
        n_expressing = int(len(sender_indices) * self.params['sender_expressing_fraction'])
        expressing_sender_indices = np.random.choice(sender_indices, n_expressing, replace=False)
        is_expressing[expressing_sender_indices] = True
        
        concentration = self._calculate_gaussian_signal_field(
            positions, expressing_sender_indices, lambda_val, config
        )
        
        expression = self._generate_cell_expression(
            positions, cell_types, concentration,
            sender_mask, receiver_mask, inert_mask, is_expressing
        )
        
        expression_log = np.log1p(expression)
        sigma = config.get_sigma(lambda_val)
        
        return expression, expression_log, is_expressing, concentration, sigma
    
    def _calculate_gaussian_signal_field(
        self,
        positions: np.ndarray,
        expressing_sender_indices: np.ndarray,
        lambda_val: float,
        config: SpatialConfig
    ) -> np.ndarray:
        n_cells = len(positions)
        concentration = np.zeros(n_cells)
        
        if len(expressing_sender_indices) == 0:
            return concentration
        
        sigma = config.get_sigma(lambda_val)
        gaussian_factor = -1.0 / (2 * sigma * sigma)
        sender_positions = positions[expressing_sender_indices]
        
        for i in range(n_cells):
            diff = positions[i] - sender_positions
            dist_sq = np.sum(diff * diff, axis=1)
            signals = np.exp(dist_sq * gaussian_factor)
            signals[signals < config.WEIGHT_THRESHOLD] = 0
            concentration[i] = np.sum(signals)
        
        if np.max(concentration) > 0:
            concentration = concentration / np.max(concentration)
        
        return concentration
    
    def _generate_cell_expression(
        self,
        positions: np.ndarray,
        cell_types: List[str],
        concentration: np.ndarray,
        sender_mask: np.ndarray,
        receiver_mask: np.ndarray,
        inert_mask: np.ndarray,
        is_expressing: np.ndarray
    ) -> np.ndarray:
        n_cells = len(positions)
        expression = np.zeros(n_cells)
        
        for i in range(n_cells):
            if sender_mask[i]:
                if is_expressing[i]:
                    mean_exp = np.random.gamma(5, 6)
                    expression[i] = self.generate_zinb_counts(
                        mean_exp, self.params['sender_high_dispersion'],
                        self.params['technical_dropout']
                    )
                else:
                    if np.random.rand() < 0.5:
                        expression[i] = self.generate_zinb_counts(3, 2.0, 0.2)
            
            elif receiver_mask[i]:
                has_baseline = np.random.rand() < self.params['receiver_baseline_fraction']
                
                if has_baseline:
                    baseline = self.generate_zinb_counts(
                        self.params['receiver_baseline_mean'], 2.0,
                        self.params['technical_dropout']
                    )
                else:
                    baseline = 0
                
                response_prob = concentration[i] * 0.8
                
                if np.random.rand() < response_prob:
                    induced_mean = concentration[i] * self.params['receiver_induced_mean']
                    induced = self.generate_zinb_counts(
                        induced_mean, 1.0, self.params['technical_dropout']
                    )
                    expression[i] = baseline + induced
                else:
                    expression[i] = baseline
            
            else:
                if np.random.rand() < self.params['inert_expression_rate']:
                    expression[i] = self.generate_zinb_counts(4, 2.0, 0.1)
        
        return expression

# =============================================================================
# Weight Matrix Builder
# =============================================================================

class WeightMatrixBuilder:
    """Build spatial weight matrices - NO distance cutoff applied"""
    
    def __init__(self, config: SpatialConfig):
        self.config = config
    
    def build_gaussian_weight_matrix(
        self,
        sender_coords: np.ndarray,
        receiver_coords: np.ndarray,
        radius: float,
        inner_radius: Optional[float] = None,
        lambda_val: Optional[float] = None
    ) -> Tuple[np.ndarray, float]:
        if lambda_val is not None:
            sigma = self.config.get_sigma(lambda_val)
        else:
            sigma = self.config.get_sigma(radius)
            
        gaussian_factor = -1.0 / (2 * sigma * sigma)
        
        distances = spatial.distance_matrix(sender_coords, receiver_coords)
        distances_sq = distances ** 2
        
        mask = distances <= radius
        if inner_radius is not None:
            mask = mask & (distances > inner_radius)
        
        W = np.exp(distances_sq * gaussian_factor) * mask
        W[W < self.config.WEIGHT_THRESHOLD] = 0
        
        total_weight = np.sum(W)
        
        row_sums = W.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        W_normalized = W / row_sums
        
        return W_normalized, total_weight

# =============================================================================
# Moran's I Calculator
# =============================================================================

class MoransICalculator:
    """Calculate I_ND for ALL distance ranges"""
    
    def __init__(self, config: SpatialConfig):
        self.config = config
        self.weight_builder = WeightMatrixBuilder(config)
    
    def calculate_I_ND(
        self,
        positions: np.ndarray,
        cell_types: List[str],
        expression: np.ndarray,
        max_radius: float = 4000,
        step: float = 50,
        lambda_val: Optional[float] = None
    ) -> Tuple[List[float], List[float]]:
        sender_mask = np.array([ct == 'Sender' for ct in cell_types])
        receiver_mask = np.array([ct == 'Receiver' for ct in cell_types])
        
        sender_pos = positions[sender_mask]
        sender_exp = expression[sender_mask]
        receiver_pos = positions[receiver_mask]
        receiver_exp = expression[receiver_mask]
        
        global_mean = np.mean(expression)
        global_std = np.std(expression)
        
        if global_std == 0:
            return [], []
        
        z_s = (sender_exp - global_mean) / (global_std + 1e-10)
        z_r = (receiver_exp - global_mean) / (global_std + 1e-10)
        
        radii = list(range(int(step), int(max_radius) + 1, int(step)))
        I_ND_values = []
        
        for radius in radii:
            inner_radius = max(0, radius - step)
            
            W_norm, total_weight = self.weight_builder.build_gaussian_weight_matrix(
                sender_pos, receiver_pos, radius, inner_radius, lambda_val=lambda_val
            )
            
            # FIXED: Set to 0.0 instead of nan
            if total_weight < 5:
                I_ND_values.append(0.0)
                continue
            
            lag = W_norm @ z_r
            
            norm_s = np.linalg.norm(z_s)
            norm_lag = np.linalg.norm(lag)
            
            if norm_s > 0 and norm_lag > 0:
                I_ND = np.dot(z_s, lag) / (norm_s * norm_lag)
            else:
                I_ND = 0.0
            
            I_ND_values.append(I_ND)
        
        return radii, I_ND_values

# =============================================================================
# Visualization - DENSITY COMPARISON
# =============================================================================

def plot_density_comparison(
    all_results: Dict,
    density_scenarios: List[Tuple[float, str]],
    lambda_val: int,
    config: SpatialConfig,
    step_size: float,
    save_path: Optional[str] = None
):
    """Plot comparison across different receiver densities for a single lambda"""
    
    fig = plt.figure(figsize=(20, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    colors = ['#e74c3c', '#f39c12', '#2ecc71']  # Red, Orange, Green
    
    # 1. I_ND comparison across densities
    ax1 = fig.add_subplot(gs[0, :2])
    
    for i, (frac, label) in enumerate(density_scenarios):
        results = all_results[frac][lambda_val]
        df = results['df']
        cutoff = config.get_max_signaling_distance(lambda_val)
        
        ax1.plot(df['radius'], df['I_ND_smooth'], 
                color=colors[i], linewidth=2.5,
                label=label, alpha=0.8)
        ax1.axvline(x=cutoff, color=colors[i], linestyle=':', alpha=0.3)
    
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax1.fill_between([0, 4000], -0.05, 0.05, alpha=0.1, color='gray')
    ax1.set_xlabel('Distance (μm)', fontsize=12)
    ax1.set_ylabel('I_ND', fontsize=12)
    ax1.set_title(f'I_ND vs Receiver Density (λ={lambda_val}μm)', 
                 fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 4000)
    
    # 2. Cell composition
    ax2 = fig.add_subplot(gs[0, 2])
    
    comp_data = []
    for frac, label in density_scenarios:
        cell_types = all_results[frac][lambda_val]['cell_types']
        n_senders = sum(1 for ct in cell_types if ct == 'Sender')
        n_receivers = sum(1 for ct in cell_types if ct == 'Receiver')
        n_inert = sum(1 for ct in cell_types if ct == 'Inert')
        total = len(cell_types)
        
        comp_data.append([n_senders/total*100, n_receivers/total*100, n_inert/total*100])
    
    comp_data = np.array(comp_data)
    x = np.arange(len(density_scenarios))
    width = 0.6
    
    ax2.bar(x, comp_data[:, 0], width, label='Senders', color='red', alpha=0.7)
    ax2.bar(x, comp_data[:, 1], width, bottom=comp_data[:, 0], 
           label='Receivers', color='blue', alpha=0.7)
    ax2.bar(x, comp_data[:, 2], width, bottom=comp_data[:, 0] + comp_data[:, 1],
           label='Inert', color='gray', alpha=0.7)
    
    ax2.set_ylabel('% of Total Cells', fontsize=11)
    ax2.set_title('Cell Type Composition', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels([label.split('(')[0].strip() for _, label in density_scenarios])
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 3-5. Spatial distributions for each density
    for idx, (frac, label) in enumerate(density_scenarios):
        ax = fig.add_subplot(gs[1, idx])
        
        positions = all_results[frac][lambda_val]['positions']
        cell_types = all_results[frac][lambda_val]['cell_types']
        counts = all_results[frac][lambda_val]['counts']
        
        sender_mask = np.array([ct == 'Sender' for ct in cell_types])
        receiver_mask = np.array([ct == 'Receiver' for ct in cell_types])
        inert_mask = np.array([ct == 'Inert' for ct in cell_types])
        
        # Plot senders
        sender_pos = positions[sender_mask]
        sender_counts = counts[sender_mask]
        ax.scatter(sender_pos[:, 0], sender_pos[:, 1],
                  c=sender_counts, s=15, alpha=0.8, cmap='Reds', 
                  vmin=0, vmax=40, edgecolor='black', linewidth=0.5)
        
        # Sample receivers
        receiver_pos = positions[receiver_mask]
        n_receivers = len(receiver_pos)
        if n_receivers > 0:
            sample_size = min(500, n_receivers)
            sample_idx = np.random.choice(n_receivers, sample_size, replace=False)
            ax.scatter(receiver_pos[sample_idx, 0], receiver_pos[sample_idx, 1],
                      c='blue', s=5, alpha=0.3, label=f'Receivers (n={n_receivers})')
        
        # Sample inert
        inert_pos = positions[inert_mask]
        n_inert = len(inert_pos)
        if n_inert > 0:
            sample_size = min(500, n_inert)
            sample_idx = np.random.choice(n_inert, sample_size, replace=False)
            ax.scatter(inert_pos[sample_idx, 0], inert_pos[sample_idx, 1],
                      c='gray', s=3, alpha=0.2, label=f'Inert (n={n_inert})')
        
        cutoff = config.get_max_signaling_distance(lambda_val)
        circle = plt.Circle((2500, 2500), cutoff, color='red', 
                           fill=False, linestyle='--', linewidth=2)
        ax.add_patch(circle)
        
        ax.set_xlim(0, 5000)
        ax.set_ylim(0, 5000)
        ax.set_aspect('equal')
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
    
    plt.suptitle(f'Receiver Density Impact on Spatial Signaling (λ={lambda_val}μm, step={step_size}μm)',
                fontsize=15, fontweight='bold')
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    plt.tight_layout()
    plt.show()

def plot_comprehensive_summary(
    all_results: Dict,
    density_scenarios: List[Tuple[float, str]],
    lambda_values: List[int],
    config: SpatialConfig,
    save_path: Optional[str] = None
):
    """Summary comparison across all densities and lambdas"""
    
    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    colors_density = ['#e74c3c', '#f39c12', '#2ecc71']
    
    # Plot for each lambda
    for lambda_idx, lambda_val in enumerate(lambda_values[:3]):  # First 3 lambdas
        ax = fig.add_subplot(gs[lambda_idx, 0])
        
        for density_idx, (frac, label) in enumerate(density_scenarios):
            df = all_results[frac][lambda_val]['df']
            ax.plot(df['radius'], df['I_ND_smooth'],
                   color=colors_density[density_idx], linewidth=2,
                   label=label, alpha=0.8)
        
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Distance (μm)', fontsize=10)
        ax.set_ylabel('I_ND', fontsize=10)
        ax.set_title(f'λ={lambda_val}μm', fontsize=11, fontweight='bold')
        if lambda_idx == 0:
            ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 2000)
    
    # Summary metrics heatmap
    ax_heat = fig.add_subplot(gs[:, 1:])
    
    # Collect peak I_ND values
    peak_data = np.zeros((len(lambda_values), len(density_scenarios)))
    
    for i, lambda_val in enumerate(lambda_values):
        for j, (frac, label) in enumerate(density_scenarios):
            df = all_results[frac][lambda_val]['df']
            peak_data[i, j] = df['I_ND_smooth'].max()
    
    im = ax_heat.imshow(peak_data, cmap='YlOrRd', aspect='auto', vmin=0, vmax=peak_data.max())
    
    ax_heat.set_xticks(np.arange(len(density_scenarios)))
    ax_heat.set_yticks(np.arange(len(lambda_values)))
    ax_heat.set_xticklabels([label.split('(')[0].strip() for _, label in density_scenarios])
    ax_heat.set_yticklabels([f'λ={lv}μm' for lv in lambda_values])
    
    # Add text annotations
    for i in range(len(lambda_values)):
        for j in range(len(density_scenarios)):
            text = ax_heat.text(j, i, f'{peak_data[i, j]:.3f}',
                              ha="center", va="center", color="black", fontsize=10)
    
    ax_heat.set_title('Peak I_ND Comparison\n(Receiver Density × Lambda)', 
                     fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax_heat, label='Peak I_ND')
    
    plt.suptitle('Comprehensive Density Analysis', fontsize=16, fontweight='bold')
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    plt.tight_layout()
    plt.show()

# =============================================================================
# Main - WITH DENSITY SCENARIOS
# =============================================================================

def run_density_simulation(
    lambda_values: List[float] = [200, 500, 1000, 2000],
    density_scenarios: List[Tuple[float, str]] = [
        (0.9, 'High Receiver (90% vs 10%)'),
        (0.5, 'Balanced (50% vs 50%)'),
        (0.1, 'Low Receiver (10% vs 90%)')
    ],
    expression_fraction: float = 0.9,
    sigma_relationship: str = 'standard',
    sigma_cutoff_multiplier: float = 3.0,
    step_size: float = 20,
    max_radius: float = 4000,
    visualize: bool = True,
    save_figures: bool = True
) -> Dict:
    """
    Run simulation across multiple receiver density scenarios
    
    Parameters:
    -----------
    density_scenarios : List[Tuple[float, str]]
        List of (receiver_fraction, label) tuples
    """
    
    print("="*80)
    print("DENSITY-BASED SIMULATION")
    print(f"Scenarios: {len(density_scenarios)}")
    for frac, label in density_scenarios:
        print(f"  - {label}")
    print(f"Lambda values: {lambda_values}")
    print(f"Step size: {step_size} μm")
    print("="*80)
    
    config = SpatialConfig(sigma_relationship, sigma_cutoff_multiplier)
    tissue_gen = TissueGenerator()
    expr_gen = ExpressionGenerator()
    moran_calc = MoransICalculator(config)
    
    expr_gen.params['sender_expressing_fraction'] = expression_fraction
    
    # Store results for all density scenarios
    all_results = {}
    
    for receiver_frac, scenario_label in density_scenarios:
        print(f"\n{'='*80}")
        print(f"SCENARIO: {scenario_label}")
        print(f"{'='*80}")
        
        # Generate tissue with specific density
        positions, cell_types, center = tissue_gen.create_tissue_with_sender_cluster(
            n_senders=500,
            n_non_senders=3500,
            receiver_fraction=receiver_frac
        )
        
        # Count cells
        n_senders = sum(1 for ct in cell_types if ct == 'Sender')
        n_receivers = sum(1 for ct in cell_types if ct == 'Receiver')
        n_inert = sum(1 for ct in cell_types if ct == 'Inert')
        
        print(f"Cell composition:")
        print(f"  Senders: {n_senders}")
        print(f"  Receivers: {n_receivers} ({n_receivers/(n_senders+n_receivers+n_inert)*100:.1f}%)")
        print(f"  Inert: {n_inert} ({n_inert/(n_senders+n_receivers+n_inert)*100:.1f}%)")
        
        results_for_density = {}
        
        for lambda_val in lambda_values:
            print(f"\n  Processing λ = {lambda_val} μm...")
            
            counts, expr_log, is_expressing, concentration, sigma = \
                expr_gen.generate_expression_with_gaussian_decay(
                    positions, cell_types, lambda_val, config
                )
            
            overall_zeros = np.sum(counts == 0) / len(counts) * 100
            print(f"    Zero counts: {overall_zeros:.1f}%")
            
            radii, I_ND_values = moran_calc.calculate_I_ND(
                positions, cell_types, expr_log,
                max_radius=max_radius,
                step=step_size,
                lambda_val=lambda_val
            )
            
            df = pd.DataFrame({'radius': radii, 'I_ND': I_ND_values})
            window = max(3, int(100 / step_size))
            df['I_ND_smooth'] = df['I_ND'].rolling(window, center=True, min_periods=1).mean()
            
            results_for_density[lambda_val] = {
                'df': df, 'counts': counts, 'expression_log': expr_log,
                'concentration': concentration, 'sigma': sigma,
                'positions': positions, 'cell_types': cell_types
            }
            
            peak_I_ND = df['I_ND_smooth'].max()
            peak_location = df.loc[df['I_ND_smooth'].idxmax(), 'radius']
            print(f"    Peak I_ND: {peak_I_ND:.3f} at {peak_location} μm")
        
        all_results[receiver_frac] = results_for_density
    
    # Visualization
    if visualize:
        # Plot for each lambda
        for lambda_val in lambda_values:
            save_path = f'density_comparison_lambda{lambda_val}.png' if save_figures else None
            plot_density_comparison(all_results, density_scenarios, lambda_val, 
                                   config, step_size, save_path)
        
        # Comprehensive summary
        save_path = 'density_comprehensive_summary.png' if save_figures else None
        plot_comprehensive_summary(all_results, density_scenarios, lambda_values, 
                                  config, save_path)
    
    print("\n" + "="*80)
    print("✓ Density simulation complete")
    print("="*80)
    
    return all_results

if __name__ == "__main__":
    results = run_density_simulation(
        lambda_values=[200, 500, 1000, 2000],
        density_scenarios=[
            (0.9, 'High Receiver (90% vs 10%)'),
            (0.5, 'Balanced (50% vs 50%)'),
            (0.1, 'Low Receiver (10% vs 90%)')
        ],
        expression_fraction=0.9,
        step_size=20,
        visualize=True,
        save_figures=True
    )

#!/usr/bin/env python3
"""
Gaussian Decay Spatial Transcriptomics Simulation
Calculate I_ND across ALL distance ranges with configurable step size

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
# Tissue Generation
# =============================================================================

class TissueGenerator:
    @staticmethod
    def create_tissue_with_sender_cluster(
        n_senders: int = 500,
        n_receivers: int = 2500, 
        n_inert: int = 1000,
        domain_size: float = 5000,
        cluster_radius: float = 600
    ) -> Tuple[np.ndarray, List[str], np.ndarray]:
        positions = []
        cell_types = []
        center = np.array([domain_size/2, domain_size/2])
        
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
        
        for _ in range(n_receivers):
            positions.append([
                np.random.uniform(0, domain_size),
                np.random.uniform(0, domain_size)
            ])
            cell_types.append('Receiver')
        
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
        """
        Calculate I_ND for ALL distances up to max_radius
        
        Parameters:
        -----------
        step : float
            Distance step size in μm (default: 50)
            Each annular ring will have width = step
        """
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
# Visualization
# =============================================================================

def plot_comprehensive_analysis(
    results: Dict, 
    lambda_values: List[float], 
    config: SpatialConfig, 
    step_size: float,
    save_path: Optional[str] = None
):
    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(3, 4, figure=fig, hspace=0.3, wspace=0.25)
    colors = plt.cm.tab10(np.linspace(0, 0.8, len(lambda_values)))
    
    # 1. I_ND curves - ALL RANGES
    ax1 = fig.add_subplot(gs[0, :2])
    for i, lambda_val in enumerate(lambda_values):
        df = results[lambda_val]['df']
        cutoff = config.get_max_signaling_distance(lambda_val)
        
        ax1.plot(df['radius'], df['I_ND_smooth'], 
                color=colors[i], linewidth=2.5,
                label=f'λ={lambda_val}μm', alpha=0.8)
        
        ax1.axvline(x=cutoff, color=colors[i], linestyle=':', alpha=0.3, linewidth=1.5)
    
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax1.fill_between([0, 4000], -0.05, 0.05, alpha=0.1, color='gray')
    ax1.set_xlabel('Distance (μm)', fontsize=11)
    ax1.set_ylabel('I_ND', fontsize=11)
    ax1.set_title(f'I_ND Across ALL Distances (step={step_size}μm, dashed=cutoff)', 
                 fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 4000)
    ax1.set_ylim(-0.5, 1.0)
    
    # 2. Decay functions
    ax2 = fig.add_subplot(gs[0, 2:])
    distances = np.linspace(0, 2000, 500)
    for i, lambda_val in enumerate(lambda_values[:3]):
        sigma = config.get_sigma(lambda_val)
        cutoff = config.get_max_signaling_distance(lambda_val)
        
        gauss_decay = np.exp(-distances**2 / (2 * sigma**2))
        ax2.plot(distances, gauss_decay, color=colors[i],
                linewidth=2, label=f'λ={lambda_val}, σ={sigma:.0f}')
        
        ax2.axvline(x=cutoff, color=colors[i], linestyle=':', alpha=0.5)
        ax2.axvspan(cutoff, 2000, alpha=0.05, color=colors[i])
    
    ax2.set_xlabel('Distance (μm)', fontsize=11)
    ax2.set_ylabel('Weight', fontsize=11)
    ax2.set_title('Gaussian Decay (shaded = beyond cutoff)', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-4, 1)
    
    # 3. Zero count statistics
    ax3 = fig.add_subplot(gs[1, 0])
    zero_data = []
    labels = []
    
    for lambda_val in lambda_values:
        counts = results[lambda_val]['counts']
        overall_zero = np.sum(counts == 0) / len(counts) * 100
        zero_data.append(overall_zero)
        labels.append(f'λ={lambda_val}')
    
    ax3.bar(range(len(lambda_values)), zero_data, color='coral', alpha=0.7)
    ax3.axhline(y=20, color='red', linestyle='--', label='Target 20%', linewidth=2)
    ax3.set_xlabel('Lambda (μm)', fontsize=11)
    ax3.set_ylabel('% Zero Counts', fontsize=11)
    ax3.set_title('Overall Zero Counts', fontsize=11, fontweight='bold')
    ax3.set_xticks(range(len(lambda_values)))
    ax3.set_xticklabels(labels)
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Count distribution
    ax4 = fig.add_subplot(gs[1, 1])
    if 500 in results:
        counts = results[500]['counts']
        cell_types = results[500]['cell_types']
        sender_mask = np.array([ct == 'Sender' for ct in cell_types])
        receiver_mask = np.array([ct == 'Receiver' for ct in cell_types])
        
        bins = np.arange(0, min(np.max(counts), 60) + 2) - 0.5
        ax4.hist(counts[sender_mask], bins=bins, alpha=0.6, 
                color='red', label='Senders', density=True)
        ax4.hist(counts[receiver_mask], bins=bins, alpha=0.6,
                color='blue', label='Receivers', density=True)
        ax4.set_xlabel('Raw Counts', fontsize=11)
        ax4.set_ylabel('Density', fontsize=11)
        ax4.set_title('Count Distribution (λ=500)', fontsize=11, fontweight='bold')
        ax4.legend(fontsize=10)
        ax4.set_xlim(-0.5, 50)
    
    # 5. Spatial visualization
    ax5 = fig.add_subplot(gs[1, 2:])
    if 500 in results:
        positions = results[500]['positions']
        cell_types = results[500]['cell_types']
        counts = results[500]['counts']
        sender_mask = np.array([ct == 'Sender' for ct in cell_types])
        receiver_mask = np.array([ct == 'Receiver' for ct in cell_types])
        
        cutoff_500 = config.get_max_signaling_distance(500)
        circle = plt.Circle((2500, 2500), cutoff_500, color='red', 
                           fill=False, linestyle='--', linewidth=2, 
                           label=f'Cutoff={cutoff_500:.0f}μm')
        ax5.add_patch(circle)
        
        sender_counts = counts[sender_mask]
        sender_pos = positions[sender_mask]
        ax5.scatter(sender_pos[:, 0], sender_pos[:, 1],
                   c=sender_counts, s=20, alpha=0.8, cmap='Reds', 
                   vmin=0, vmax=40, edgecolor='black', linewidth=0.5,
                   label='Senders')
        
        receiver_counts = counts[receiver_mask]
        receiver_pos = positions[receiver_mask]
        ax5.scatter(receiver_pos[::5, 0], receiver_pos[::5, 1],
                   c=receiver_counts[::5], s=8, alpha=0.5, cmap='Blues', 
                   vmin=0, vmax=30, label='Receivers')
        
        ax5.set_xlim(0, 5000)
        ax5.set_ylim(0, 5000)
        ax5.set_aspect('equal')
        ax5.set_title('Spatial Distribution (λ=500)', fontsize=11, fontweight='bold')
        ax5.legend(fontsize=9, loc='upper right')
    
    # 6. Summary table
    ax6 = fig.add_subplot(gs[2, :])
    summary_data = []
    
    for lambda_val in lambda_values:
        df = results[lambda_val]['df']
        sigma = results[lambda_val]['sigma']
        cutoff = config.get_max_signaling_distance(lambda_val)
        counts = results[lambda_val]['counts']
        
        overall_zero = np.sum(counts == 0) / len(counts) * 100
        peak = df['I_ND_smooth'].max()
        peak_loc = df.loc[df['I_ND_smooth'].idxmax(), 'radius']
        
        within = df[df['radius'] <= cutoff]['I_ND_smooth']
        beyond = df[df['radius'] > cutoff]['I_ND_smooth']
        
        below_zero = df[df['I_ND_smooth'] < 0.05]
        if len(below_zero) > 0:
            decay_dist = below_zero.iloc[0]['radius']
        else:
            decay_dist = '>4000'
        
        summary_data.append([
            f'{lambda_val}',
            f'{sigma:.0f}',
            f'{cutoff:.0f}',
            f'{overall_zero:.1f}%',
            f'{peak:.3f}',
            f'{peak_loc}',
            f'{within.mean():.3f}' if len(within) > 0 else 'N/A',
            f'{beyond.mean():.4f}' if len(beyond) > 0 else 'N/A',
            f'{decay_dist}'
        ])
    
    table = ax6.table(cellText=summary_data,
                     colLabels=['λ', 'σ', 'Cutoff', 'Zero%', 'Peak', 'Peak@', 
                               'I_ND≤cut', 'I_ND>cut', 'Decay@'],
                     cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)
    
    for i in range(len(summary_data) + 1):
        for j in range(9):
            cell = table[(i, j)]
            if i == 0:
                cell.set_facecolor('#E8E8E8')
                cell.set_text_props(weight='bold')
            else:
                cell.set_facecolor('#F5F5F5' if i % 2 == 0 else 'white')
    
    ax6.axis('off')
    ax6.set_title(f'Summary Metrics (step={step_size}μm)', fontsize=12, fontweight='bold')
    
    plt.suptitle(f'I_ND Calculated Across ALL Distance Ranges (step={step_size}μm)',
                fontsize=16, fontweight='bold')
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    plt.tight_layout()
    plt.show()

# =============================================================================
# Main
# =============================================================================

def run_simulation(
    lambda_values: List[float] = [200, 500, 1000, 2000],
    expression_fraction: float = 0.9,
    sigma_relationship: str = 'standard',
    sigma_cutoff_multiplier: float = 3.0,
    step_size: float = 50,
    max_radius: float = 4000,
    visualize: bool = True,
    save_figures: bool = True
) -> Dict:
    """
    Run simulation calculating I_ND at ALL distances
    
    Parameters:
    -----------
    step_size : float
        Distance step size in μm (default: 50)
        Ring width for each calculation = step_size
        Calculates at: step_size, 2*step_size, 3*step_size, ...
    max_radius : float
        Maximum distance to calculate (default: 4000)
    """
    
    print("="*80)
    print("CALCULATE I_ND ACROSS ALL DISTANCE RANGES")
    print(f"Step size: {step_size} μm (ring width = {step_size} μm)")
    print(f"Max radius: {max_radius} μm")
    print("="*80)
    
    config = SpatialConfig(sigma_relationship, sigma_cutoff_multiplier)
    tissue_gen = TissueGenerator()
    expr_gen = ExpressionGenerator()
    moran_calc = MoransICalculator(config)
    
    positions, cell_types, center = tissue_gen.create_tissue_with_sender_cluster()
    results = {}
    expr_gen.params['sender_expressing_fraction'] = expression_fraction
    
    for lambda_val in lambda_values:
        print(f"\nProcessing λ = {lambda_val} μm...")
        sigma = config.get_sigma(lambda_val)
        cutoff = config.get_max_signaling_distance(lambda_val)
        print(f"  σ = {sigma:.0f} μm, Theoretical cutoff = {cutoff:.0f} μm")
        
        counts, expr_log, is_expressing, concentration, sigma = \
            expr_gen.generate_expression_with_gaussian_decay(
                positions, cell_types, lambda_val, config
            )
        
        overall_zeros = np.sum(counts == 0) / len(counts) * 100
        print(f"  Zero counts: {overall_zeros:.1f}%")
        
        # Calculate with specified step size
        radii, I_ND_values = moran_calc.calculate_I_ND(
            positions, cell_types, expr_log, 
            max_radius=max_radius,
            step=step_size,
            lambda_val=lambda_val
        )
        
        df = pd.DataFrame({'radius': radii, 'I_ND': I_ND_values})
        
        # Smooth based on step size (adjust window)
        window = max(3, int(100 / step_size))  # ~100μm smoothing window
        df['I_ND_smooth'] = df['I_ND'].rolling(window, center=True, min_periods=1).mean()
        
        results[lambda_val] = {
            'df': df, 'counts': counts, 'expression_log': expr_log,
            'concentration': concentration, 'sigma': sigma,
            'positions': positions, 'cell_types': cell_types
        }
        
        peak_I_ND = df['I_ND_smooth'].max()
        peak_location = df.loc[df['I_ND_smooth'].idxmax(), 'radius']
        print(f"  Peak I_ND: {peak_I_ND:.3f} at {peak_location} μm")
        
        short = df[df['radius'] <= 500]['I_ND_smooth'].mean()
        medium = df[(df['radius'] > 500) & (df['radius'] <= 2000)]['I_ND_smooth'].mean()
        long_range = df[df['radius'] > 2000]['I_ND_smooth'].mean()
        
        print(f"  I_ND 0-500μm: {short:.3f}")
        print(f"  I_ND 500-2000μm: {medium:.3f}")
        print(f"  I_ND >2000μm: {long_range:.4f}")
    
    if visualize:
        save_path = f'gaussian_simulation_step{int(step_size)}.png' if save_figures else None
        plot_comprehensive_analysis(results, lambda_values, config, step_size, save_path)
    
    print("\n" + "="*80)
    print(f"✓ I_ND calculated with {step_size}μm steps across all distances")
    print("="*80)
    
    return results

if __name__ == "__main__":
    results = run_simulation(
        lambda_values=[200, 500, 1000, 2000, 3000],
        expression_fraction=1.0,
        sigma_relationship='standard',
        sigma_cutoff_multiplier=3.0,
        step_size=20,  # CONFIGURABLE: 20 μm bins with 20 μm ring width
        max_radius=5000,  # CONFIGURABLE: maximum distance
        visualize=True,
        save_figures=True
    )

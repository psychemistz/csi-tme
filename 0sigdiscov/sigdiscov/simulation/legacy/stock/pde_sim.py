#!/usr/bin/env python3
"""
Comprehensive Simulation Framework for Validating Normalized Directional Moran's I
Detection of Cytokine-Mediated Cell-to-Cell Communication

Author: Seongyong Park
Institution: CDSL, NCI, NIH
Date: 17th. Nov.2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats, spatial, ndimage
from scipy.sparse import csr_matrix
from scipy.optimize import minimize_scalar
from typing import Tuple, Dict, List, Optional
import warnings
from dataclasses import dataclass
import multiprocessing as mp
from functools import partial

# Simple progress indicator if tqdm not available
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, desc=None, **kwargs):
        if desc:
            print(f"{desc}...")
        return iterable

warnings.filterwarnings('ignore')

# Set random seed for reproducibility
np.random.seed(42)

# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class CytokineParams:
    """Parameters for cytokine diffusion and signaling"""
    name: str
    diffusion_coefficient: float  # μm²/s
    degradation_rate: float      # 1/min
    secretion_rate: float        # molecules/cell/hour
    Kd: float                    # nM (dissociation constant)
    hill_coefficient: float      # Hill coefficient for dose-response
    max_fold_change: float       # Maximum expression fold change
    
# Realistic cytokine parameters from literature
CYTOKINE_LIBRARY = {
    'IL-2': CytokineParams(
        name='IL-2',
        diffusion_coefficient=20.0,
        degradation_rate=0.023,  # t½ ≈ 30 min
        secretion_rate=1000.0,
        Kd=0.01,  # 10 pM
        hill_coefficient=1.5,
        max_fold_change=5.0
    ),
    'IL-6': CytokineParams(
        name='IL-6',
        diffusion_coefficient=30.0,
        degradation_rate=0.014,  # t½ ≈ 50 min
        secretion_rate=2000.0,
        Kd=0.05,  # 50 pM
        hill_coefficient=1.8,
        max_fold_change=8.0
    ),
    'TNF': CytokineParams(
        name='TNF',
        diffusion_coefficient=40.0,
        degradation_rate=0.035,  # t½ ≈ 20 min
        secretion_rate=500.0,
        Kd=0.02,  # 20 pM
        hill_coefficient=2.0,
        max_fold_change=10.0
    ),
    'CXCL12': CytokineParams(
        name='CXCL12',
        diffusion_coefficient=100.0,
        degradation_rate=0.012,  # t½ ≈ 60 min
        secretion_rate=3000.0,
        Kd=5.0,  # 5 nM
        hill_coefficient=1.2,
        max_fold_change=4.0
    )
}

# =============================================================================
# Core PDE Solver for Cytokine Diffusion
# =============================================================================

class CytokineDiffusionSimulator:
    """
    Solves reaction-diffusion PDE for cytokine concentration in tissue
    ∂C/∂t = D∇²C - ηC + S - U
    """
    
    def __init__(self, 
                 domain_size: Tuple[float, float] = (1000, 1000),  # μm
                 grid_spacing: float = 10.0,  # μm
                 dt: float = 0.1):  # seconds
        
        self.domain_size = domain_size
        self.grid_spacing = grid_spacing
        self.dt = dt
        
        # Create computational grid
        self.nx = int(domain_size[0] / grid_spacing)
        self.ny = int(domain_size[1] / grid_spacing)
        self.x = np.linspace(0, domain_size[0], self.nx)
        self.y = np.linspace(0, domain_size[1], self.ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
    def solve_steady_state(self,
                          producer_positions: np.ndarray,
                          consumer_positions: np.ndarray,
                          cytokine_params: CytokineParams,
                          consumer_density: Optional[np.ndarray] = None,
                          max_time: float = 3600.0,  # seconds (1 hour)
                          tolerance: float = 1e-6) -> np.ndarray:
        """
        Solve diffusion equation to steady state
        
        Returns:
            concentration: 2D array of cytokine concentration (nM)
        """
        
        # Initialize concentration field
        C = np.zeros((self.ny, self.nx))
        
        # Convert parameters to appropriate units
        D = cytokine_params.diffusion_coefficient  # μm²/s
        eta = cytokine_params.degradation_rate / 60.0  # convert to 1/s
        S = cytokine_params.secretion_rate / 3600.0  # molecules/cell/s
        
        # Stability criterion for explicit scheme
        max_dt = 0.25 * self.grid_spacing**2 / D
        if self.dt > max_dt:
            self.dt = max_dt * 0.9
            
        # Create source term (producers)
        source = np.zeros((self.ny, self.nx))
        for pos in producer_positions:
            ix = int(pos[0] / self.grid_spacing)
            iy = int(pos[1] / self.grid_spacing)
            if 0 <= ix < self.nx and 0 <= iy < self.ny:
                source[iy, ix] += S
                
        # Create sink term (consumers) - simplified uptake
        sink = np.zeros((self.ny, self.nx))
        if consumer_density is not None:
            sink = consumer_density * 1e-4  # uptake rate
        
        # Time evolution
        t = 0
        converged = False
        
        while t < max_time and not converged:
            C_old = C.copy()
            
            # Finite difference for Laplacian (5-point stencil)
            laplacian = ndimage.laplace(C) / self.grid_spacing**2
            
            # Update concentration
            dC_dt = D * laplacian - eta * C + source - sink * C
            C = C + self.dt * dC_dt
            
            # Neumann boundary conditions (no flux)
            C[0, :] = C[1, :]
            C[-1, :] = C[-2, :]
            C[:, 0] = C[:, 1]
            C[:, -1] = C[:, -2]
            
            # Check convergence
            if np.max(np.abs(C - C_old)) < tolerance:
                converged = True
                
            t += self.dt
            
        # Convert to nM (assuming appropriate scaling)
        C = C * 1e3  # arbitrary scaling for demonstration
        
        return C
    
    def sample_at_positions(self, 
                           concentration: np.ndarray,
                           positions: np.ndarray) -> np.ndarray:
        """Sample concentration field at cell positions"""
        
        values = np.zeros(len(positions))
        for i, pos in enumerate(positions):
            ix = int(pos[0] / self.grid_spacing)
            iy = int(pos[1] / self.grid_spacing)
            if 0 <= ix < self.nx and 0 <= iy < self.ny:
                values[i] = concentration[iy, ix]
        return values

# =============================================================================
# Cell Placement and Tissue Architecture
# =============================================================================

def generate_cell_positions(n_cells: int,
                           pattern: str = 'random',
                           domain_size: Tuple[float, float] = (1000, 1000),
                           **kwargs) -> np.ndarray:
    """
    Generate cell positions with different spatial patterns
    
    Patterns:
        - 'random': Uniform random distribution
        - 'clustered': Thomas cluster process
        - 'gradient': Density gradient from edge
        - 'two_regions': Two distinct regions
    """
    
    if pattern == 'random':
        positions = np.random.rand(n_cells, 2) * domain_size
        
    elif pattern == 'clustered':
        # Thomas cluster process
        n_clusters = kwargs.get('n_clusters', 10)
        cluster_std = kwargs.get('cluster_std', 50)
        
        positions = []
        cells_per_cluster = n_cells // n_clusters
        
        for _ in range(n_clusters):
            center = np.random.rand(2) * domain_size
            cluster_positions = np.random.randn(cells_per_cluster, 2) * cluster_std + center
            positions.append(cluster_positions)
            
        positions = np.vstack(positions)
        # Clip to domain
        positions = np.clip(positions, 0, domain_size)
        
    elif pattern == 'gradient':
        # Exponential gradient from left edge
        positions = []
        for _ in range(n_cells):
            x = np.random.exponential(scale=200)  # decay length 200 μm
            y = np.random.rand() * domain_size[1]
            positions.append([x, y])
        positions = np.array(positions)
        positions = np.clip(positions, 0, domain_size)
        
    elif pattern == 'two_regions':
        # Two distinct regions separated by boundary
        boundary = domain_size[0] / 2
        n_left = n_cells // 2
        n_right = n_cells - n_left
        
        left_positions = np.random.rand(n_left, 2) * [boundary, domain_size[1]]
        right_positions = np.random.rand(n_right, 2) * [boundary, domain_size[1]]
        right_positions[:, 0] += boundary
        
        positions = np.vstack([left_positions, right_positions])
        
    else:
        raise ValueError(f"Unknown pattern: {pattern}")
        
    return positions[:n_cells]  # Ensure exact number

# =============================================================================
# Gene Expression Response Models
# =============================================================================

def hill_function(concentration: np.ndarray,
                 Kd: float,
                 n: float = 1.5) -> np.ndarray:
    """Hill equation for dose-response"""
    return concentration**n / (Kd**n + concentration**n)

def generate_expression_response(concentrations: np.ndarray,
                                cytokine_params: CytokineParams,
                                baseline_expression: float = 1.0,
                                noise_cv: float = 0.3) -> np.ndarray:
    """
    Generate gene expression response to cytokine concentration
    """
    
    # Calculate receptor activation
    activation = hill_function(concentrations, 
                              cytokine_params.Kd, 
                              cytokine_params.hill_coefficient)
    
    # Generate expression with fold change
    expression = baseline_expression * (1 + cytokine_params.max_fold_change * activation)
    
    # Add biological noise (log-normal)
    if noise_cv > 0:
        noise = np.random.lognormal(mean=0, sigma=noise_cv, size=len(expression))
        expression = expression * noise
        
    return expression

# =============================================================================
# Normalized Directional Moran's I Implementation
# =============================================================================

class DirectionalMoranI:
    """
    Calculate normalized directional Moran's I statistic
    """
    
    def __init__(self, 
                 positions: np.ndarray,
                 bandwidth: float = 150.0,  # μm
                 n_directions: int = 8):
        
        self.positions = positions
        self.n_cells = len(positions)
        self.bandwidth = bandwidth
        self.n_directions = n_directions
        
        # Precompute distance matrix and angles
        self.distances = spatial.distance_matrix(positions, positions)
        self.angles = self._compute_angles()
        
        # Create spatial weights with distance decay
        self.weights = np.exp(-self.distances / bandwidth)
        np.fill_diagonal(self.weights, 0)
        
        # Normalize weights (row standardization)
        row_sums = self.weights.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        self.weights = self.weights / row_sums
        
    def _compute_angles(self) -> np.ndarray:
        """Compute angles between all cell pairs"""
        angles = np.zeros((self.n_cells, self.n_cells))
        
        for i in range(self.n_cells):
            for j in range(self.n_cells):
                if i != j:
                    dx = self.positions[j, 0] - self.positions[i, 0]
                    dy = self.positions[j, 1] - self.positions[i, 1]
                    angles[i, j] = np.arctan2(dy, dx) * 180 / np.pi
                    if angles[i, j] < 0:
                        angles[i, j] += 360
                        
        return angles
    
    def calculate_directional_I(self, 
                               expression: np.ndarray,
                               direction: float,
                               angular_width: float = 45.0) -> float:
        """
        Calculate Moran's I for a specific direction
        
        Args:
            expression: Gene expression values
            direction: Direction in degrees (0-360)
            angular_width: Width of angular bin in degrees
        """
        
        # Create directional weights
        dir_weights = self.weights.copy()
        
        # Mask weights outside angular range
        for i in range(self.n_cells):
            for j in range(self.n_cells):
                if i != j:
                    angle_diff = np.abs(self.angles[i, j] - direction)
                    # Handle circular difference
                    if angle_diff > 180:
                        angle_diff = 360 - angle_diff
                    if angle_diff > angular_width / 2:
                        dir_weights[i, j] = 0
                        
        # Calculate Moran's I with directional weights
        z = expression - np.mean(expression)
        
        numerator = np.sum(dir_weights * np.outer(z, z))
        denominator = np.sum(z**2)
        
        n = self.n_cells
        W = np.sum(dir_weights)
        
        if W > 0 and denominator > 0:
            I = (n / W) * (numerator / denominator)
        else:
            I = 0
            
        return I
    
    def calculate_all_directions(self, 
                                expression: np.ndarray) -> Dict[float, float]:
        """Calculate Moran's I for all directions"""
        
        directions = np.linspace(0, 360, self.n_directions, endpoint=False)
        I_values = {}
        
        for direction in directions:
            I_values[direction] = self.calculate_directional_I(expression, direction)
            
        return I_values
    
    def permutation_test(self,
                        expression: np.ndarray,
                        direction: float,
                        n_permutations: int = 999) -> Tuple[float, float]:
        """
        Permutation test for significance of directional Moran's I
        
        Returns:
            (I_observed, p_value)
        """
        
        I_observed = self.calculate_directional_I(expression, direction)
        
        # Generate null distribution
        I_null = []
        for _ in range(n_permutations):
            expr_perm = np.random.permutation(expression)
            I_perm = self.calculate_directional_I(expr_perm, direction)
            I_null.append(I_perm)
            
        I_null = np.array(I_null)
        
        # Calculate p-value (two-tailed)
        p_value = np.sum(np.abs(I_null) >= np.abs(I_observed)) / n_permutations
        
        return I_observed, p_value
    
    def find_peak_direction(self, 
                           expression: np.ndarray) -> Tuple[float, float]:
        """
        Find direction with maximum Moran's I
        
        Returns:
            (peak_direction, peak_I_value)
        """
        
        I_values = self.calculate_all_directions(expression)
        peak_dir = max(I_values, key=I_values.get)
        
        # Refine with finer search
        def neg_I(direction):
            return -self.calculate_directional_I(expression, direction % 360)
        
        result = minimize_scalar(neg_I, 
                               bounds=(peak_dir - 22.5, peak_dir + 22.5),
                               method='bounded')
        
        refined_dir = result.x % 360
        refined_I = -result.fun
        
        return refined_dir, refined_I

# =============================================================================
# Validation Scenarios
# =============================================================================

class ValidationScenario:
    """Base class for validation scenarios"""
    
    def __init__(self, 
                 n_cells: int = 1000,
                 domain_size: Tuple[float, float] = (1000, 1000)):
        
        self.n_cells = n_cells
        self.domain_size = domain_size
        self.simulator = CytokineDiffusionSimulator(domain_size)
        
    def run(self) -> Dict:
        """Run validation scenario and return results"""
        raise NotImplementedError

class DirectionalGradientValidation(ValidationScenario):
    """Test 1: Unidirectional cytokine gradient detection"""
    
    def run(self, cytokine: str = 'IL-2', noise_cv: float = 0.3) -> Dict:
        
        cytokine_params = CYTOKINE_LIBRARY[cytokine]
        
        # Place producers on left edge
        n_producers = int(self.n_cells * 0.1)
        producer_positions = np.column_stack([
            np.random.rand(n_producers) * 100,  # x: 0-100 μm
            np.random.rand(n_producers) * self.domain_size[1]  # y: full range
        ])
        
        # Place other cells throughout domain
        all_positions = generate_cell_positions(self.n_cells, 'random', self.domain_size)
        
        # Simulate cytokine diffusion
        concentration_field = self.simulator.solve_steady_state(
            producer_positions,
            all_positions,
            cytokine_params
        )
        
        # Sample concentration at cell positions
        concentrations = self.simulator.sample_at_positions(concentration_field, all_positions)
        
        # Generate expression response
        expression = generate_expression_response(concentrations, cytokine_params, noise_cv=noise_cv)
        
        # Calculate directional Moran's I
        moran_calc = DirectionalMoranI(all_positions)
        
        # Test different directions
        I_0 = moran_calc.calculate_directional_I(expression, 0)      # Along gradient
        I_90 = moran_calc.calculate_directional_I(expression, 90)    # Perpendicular
        I_180 = moran_calc.calculate_directional_I(expression, 180)  # Against gradient
        
        # Find peak direction
        peak_dir, peak_I = moran_calc.find_peak_direction(expression)
        
        # Significance test
        _, p_value = moran_calc.permutation_test(expression, peak_dir, n_permutations=199)
        
        return {
            'scenario': 'Directional Gradient',
            'cytokine': cytokine,
            'true_direction': 0,
            'I_along_gradient': I_0,
            'I_perpendicular': I_90,
            'I_against_gradient': I_180,
            'peak_direction': peak_dir,
            'peak_I': peak_I,
            'angular_error': min(abs(peak_dir - 0), 360 - abs(peak_dir - 0)),
            'p_value': p_value,
            'gradient_detected': I_0 > I_90 * 1.5 and I_0 > I_180 * 2.0,
            'positions': all_positions,
            'expression': expression,
            'concentration_field': concentration_field
        }

class DistanceDependentValidation(ValidationScenario):
    """Test 2: Distance-dependent signaling range"""
    
    def run(self, cytokine: str = 'IL-6', noise_cv: float = 0.3) -> Dict:
        
        cytokine_params = CYTOKINE_LIBRARY[cytokine]
        
        # Place producers in center
        n_producers = int(self.n_cells * 0.05)
        center = np.array(self.domain_size) / 2
        producer_positions = np.random.randn(n_producers, 2) * 50 + center
        
        # Place receivers at different distances
        distance_ranges = [(50, 100), (100, 150), (150, 200), (200, 300), (300, 400)]
        I_by_distance = []
        
        for min_dist, max_dist in distance_ranges:
            # Generate cells in annulus
            n_cells_ring = 200
            angles = np.random.rand(n_cells_ring) * 2 * np.pi
            radii = np.random.uniform(min_dist, max_dist, n_cells_ring)
            
            positions = np.column_stack([
                center[0] + radii * np.cos(angles),
                center[1] + radii * np.sin(angles)
            ])
            
            # Simulate diffusion
            concentration_field = self.simulator.solve_steady_state(
                producer_positions,
                positions,
                cytokine_params
            )
            
            # Sample and generate expression
            concentrations = self.simulator.sample_at_positions(concentration_field, positions)
            expression = generate_expression_response(concentrations, cytokine_params, noise_cv=noise_cv)
            
            # Calculate Moran's I
            moran_calc = DirectionalMoranI(positions, bandwidth=(max_dist - min_dist) / 2)
            I_values = moran_calc.calculate_all_directions(expression)
            mean_I = np.mean(list(I_values.values()))
            
            I_by_distance.append({
                'distance': (min_dist + max_dist) / 2,
                'mean_I': mean_I,
                'max_I': max(I_values.values())
            })
            
        # Calculate characteristic length scale
        D = cytokine_params.diffusion_coefficient
        eta = cytokine_params.degradation_rate / 60  # to 1/s
        characteristic_length = np.sqrt(D / eta)
        
        return {
            'scenario': 'Distance Dependent',
            'cytokine': cytokine,
            'I_by_distance': I_by_distance,
            'characteristic_length': characteristic_length,
            'diffusion_coefficient': D,
            'degradation_rate': eta
        }

class NegativeControlValidation(ValidationScenario):
    """Test 3: Negative controls - what should NOT be detected"""
    
    def run(self, noise_cv: float = 0.3) -> Dict:
        
        results = {}
        
        # Control 1: Random expression (no spatial pattern)
        positions = generate_cell_positions(self.n_cells, 'random', self.domain_size)
        expression = np.random.lognormal(0, noise_cv, self.n_cells)
        
        moran_calc = DirectionalMoranI(positions)
        I_values_random = moran_calc.calculate_all_directions(expression)
        cv_random = np.std(list(I_values_random.values())) / np.mean(list(I_values_random.values()))
        
        results['random_expression'] = {
            'I_values': I_values_random,
            'CV': cv_random,
            'isotropic': cv_random < 0.2  # Should be isotropic
        }
        
        # Control 2: Isotropic clustering (no directional signal)
        positions = generate_cell_positions(self.n_cells, 'clustered', self.domain_size,
                                          n_clusters=10, cluster_std=50)
        
        # Expression correlated with local density but not directional
        tree = spatial.cKDTree(positions)
        local_density = np.array([len(tree.query_ball_point(pos, 100)) for pos in positions])
        expression = local_density + np.random.normal(0, noise_cv, self.n_cells)
        
        moran_calc = DirectionalMoranI(positions)
        I_values_cluster = moran_calc.calculate_all_directions(expression)
        cv_cluster = np.std(list(I_values_cluster.values())) / np.mean(list(I_values_cluster.values()))
        
        results['clustered_isotropic'] = {
            'I_values': I_values_cluster,
            'CV': cv_cluster,
            'isotropic': cv_cluster < 0.2
        }
        
        # Control 3: Contact-dependent signaling (very short range)
        positions = generate_cell_positions(self.n_cells, 'random', self.domain_size)
        
        # Only cells within 20 μm show elevated expression
        contact_threshold = 20
        expression = np.ones(self.n_cells)
        
        for i in range(self.n_cells):
            neighbors = np.where(moran_calc.distances[i] < contact_threshold)[0]
            if len(neighbors) > 2:
                expression[i] = 5.0
                
        expression = expression * np.random.lognormal(0, noise_cv, self.n_cells)
        
        moran_calc_contact = DirectionalMoranI(positions, bandwidth=30)  # Short bandwidth
        I_values_contact = moran_calc_contact.calculate_all_directions(expression)
        
        results['contact_dependent'] = {
            'I_values': I_values_contact,
            'mean_I': np.mean(list(I_values_contact.values())),
            'weak_signal': np.mean(list(I_values_contact.values())) < 0.1
        }
        
        return {
            'scenario': 'Negative Controls',
            'results': results
        }

class ChemokineGradientValidation(ValidationScenario):
    """Validation: CXCL12 chemokine gradient directing migration"""
    
    def run(self, noise_cv: float = 0.3) -> Dict:
        
        cytokine_params = CYTOKINE_LIBRARY['CXCL12']
        
        # Stromal cells at tissue boundary produce CXCL12
        n_stromal = int(self.n_cells * 0.1)
        stromal_positions = np.column_stack([
            np.ones(n_stromal) * self.domain_size[0] * 0.9,  # Right edge
            np.random.rand(n_stromal) * self.domain_size[1]
        ])
        
        # CXCR4+ lymphocytes distributed throughout
        lymphocyte_positions = generate_cell_positions(
            int(self.n_cells * 0.6), 'random', self.domain_size
        )
        
        # Other cells
        other_positions = generate_cell_positions(
            self.n_cells - n_stromal - len(lymphocyte_positions),
            'random', self.domain_size
        )
        
        all_positions = np.vstack([stromal_positions, lymphocyte_positions, other_positions])
        
        # Cell type labels
        cell_types = np.concatenate([
            np.ones(n_stromal) * 0,  # Stromal
            np.ones(len(lymphocyte_positions)) * 1,  # CXCR4+ lymphocytes
            np.ones(len(other_positions)) * 2  # Others
        ])
        
        # Simulate CXCL12 gradient
        concentration_field = self.simulator.solve_steady_state(
            stromal_positions,
            all_positions,
            cytokine_params
        )
        
        # Sample concentrations
        concentrations = self.simulator.sample_at_positions(concentration_field, all_positions)
        
        # Generate expression response (only in CXCR4+ cells)
        expression = np.ones(len(all_positions))
        lymphocyte_mask = cell_types == 1
        expression[lymphocyte_mask] = generate_expression_response(
            concentrations[lymphocyte_mask],
            cytokine_params,
            noise_cv=noise_cv
        )
        
        # Calculate directional Moran's I
        moran_calc = DirectionalMoranI(all_positions)
        
        # Expected direction: from tissue edge (right) inward (left) = 180°
        expected_direction = 180
        
        I_values = moran_calc.calculate_all_directions(expression)
        peak_dir, peak_I = moran_calc.find_peak_direction(expression)
        
        # Test significance
        _, p_value = moran_calc.permutation_test(expression, peak_dir, n_permutations=199)
        
        # Angular error
        angular_error = min(abs(peak_dir - expected_direction), 
                          360 - abs(peak_dir - expected_direction))
        
        return {
            'scenario': 'CXCL12 Chemokine Gradient',
            'expected_direction': expected_direction,
            'detected_direction': peak_dir,
            'angular_error': angular_error,
            'peak_I': peak_I,
            'p_value': p_value,
            'direction_correct': angular_error < 30,
            'I_values': I_values,
            'positions': all_positions,
            'cell_types': cell_types,
            'expression': expression,
            'concentration_field': concentration_field
        }

# =============================================================================
# Statistical Power Analysis
# =============================================================================

def power_analysis(n_simulations: int = 100,
                  fold_changes: List[float] = [1.5, 2.0, 3.0, 5.0],
                  noise_levels: List[float] = [0.2, 0.3, 0.5, 0.8],
                  n_cells: int = 1000,
                  alpha: float = 0.05) -> pd.DataFrame:
    """
    Comprehensive power analysis for cytokine detection
    """
    
    results = []
    
    for fold_change in fold_changes:
        for noise_cv in noise_levels:
            
            detected = 0
            angular_errors = []
            
            for sim in tqdm(range(n_simulations), 
                          desc=f"FC={fold_change}, Noise={noise_cv}"):
                
                # Create custom cytokine with specified fold change
                custom_cytokine = CytokineParams(
                    name='Custom',
                    diffusion_coefficient=30.0,
                    degradation_rate=0.02,
                    secretion_rate=1000.0,
                    Kd=0.05,
                    hill_coefficient=1.5,
                    max_fold_change=fold_change
                )
                
                # Generate gradient scenario
                validator = DirectionalGradientValidation(n_cells=n_cells)
                
                # Override cytokine params
                CYTOKINE_LIBRARY['Custom'] = custom_cytokine
                result = validator.run('Custom', noise_cv)
                
                # Check detection
                if result['p_value'] < alpha:
                    detected += 1
                    angular_errors.append(result['angular_error'])
                    
            power = detected / n_simulations
            mean_angular_error = np.mean(angular_errors) if angular_errors else np.nan
            
            results.append({
                'fold_change': fold_change,
                'noise_cv': noise_cv,
                'power': power,
                'n_detected': detected,
                'mean_angular_error': mean_angular_error,
                'n_simulations': n_simulations
            })
            
    return pd.DataFrame(results)

# =============================================================================
# Visualization Functions
# =============================================================================

def visualize_cytokine_field(concentration_field: np.ndarray,
                            positions: np.ndarray,
                            expression: np.ndarray,
                            title: str = "Cytokine Diffusion and Response"):
    """Visualize cytokine concentration field and cellular response"""
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Cytokine concentration field
    im1 = axes[0].imshow(concentration_field, cmap='YlOrRd', origin='lower')
    axes[0].set_title('Cytokine Concentration Field')
    axes[0].set_xlabel('X (μm)')
    axes[0].set_ylabel('Y (μm)')
    plt.colorbar(im1, ax=axes[0], label='Concentration (nM)')
    
    # Cell positions colored by expression
    scatter = axes[1].scatter(positions[:, 0], positions[:, 1],
                             c=expression, cmap='viridis', s=20)
    axes[1].set_title('Cell Expression Response')
    axes[1].set_xlabel('X (μm)')
    axes[1].set_ylabel('Y (μm)')
    plt.colorbar(scatter, ax=axes[1], label='Expression (AU)')
    
    # Directional Moran's I polar plot
    moran_calc = DirectionalMoranI(positions)
    I_values = moran_calc.calculate_all_directions(expression)
    
    directions = np.array(list(I_values.keys())) * np.pi / 180
    I_vals = list(I_values.values())
    
    ax3 = plt.subplot(133, projection='polar')
    ax3.plot(directions, I_vals, 'b-o')
    ax3.fill_between(directions, 0, I_vals, alpha=0.3)
    ax3.set_title("Directional Moran's I")
    ax3.set_theta_zero_location('E')
    
    plt.suptitle(title)
    plt.tight_layout()
    
    return fig

def plot_power_curves(power_df: pd.DataFrame):
    """Plot power analysis results"""
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Power vs fold change for different noise levels
    for noise in power_df['noise_cv'].unique():
        subset = power_df[power_df['noise_cv'] == noise]
        axes[0].plot(subset['fold_change'], subset['power'], 
                    marker='o', label=f'CV={noise}')
    
    axes[0].axhline(y=0.8, color='r', linestyle='--', label='80% Power')
    axes[0].set_xlabel('Fold Change')
    axes[0].set_ylabel('Statistical Power')
    axes[0].set_title('Power vs Effect Size')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Angular error vs noise
    for fc in power_df['fold_change'].unique():
        subset = power_df[power_df['fold_change'] == fc]
        axes[1].plot(subset['noise_cv'], subset['mean_angular_error'],
                    marker='s', label=f'FC={fc}')
    
    axes[1].axhline(y=30, color='r', linestyle='--', label='30° threshold')
    axes[1].set_xlabel('Noise (CV)')
    axes[1].set_ylabel('Mean Angular Error (degrees)')
    axes[1].set_title('Directional Accuracy vs Noise')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

# =============================================================================
# Main Validation Pipeline
# =============================================================================

def run_complete_validation():
    """Run complete validation suite"""
    
    print("=" * 60)
    print("NORMALIZED DIRECTIONAL MORAN'S I VALIDATION")
    print("Cytokine-Mediated Cell Communication Detection")
    print("=" * 60)
    
    results = {}
    
    # Test 1: Directional Gradient Detection
    print("\n1. Testing Directional Gradient Detection...")
    gradient_validator = DirectionalGradientValidation()
    gradient_result = gradient_validator.run('IL-2', noise_cv=0.3)
    results['gradient'] = gradient_result
    
    print(f"   True direction: {gradient_result['true_direction']}°")
    print(f"   Detected direction: {gradient_result['peak_direction']:.1f}°")
    print(f"   Angular error: {gradient_result['angular_error']:.1f}°")
    print(f"   I along gradient: {gradient_result['I_along_gradient']:.3f}")
    print(f"   I perpendicular: {gradient_result['I_perpendicular']:.3f}")
    print(f"   p-value: {gradient_result['p_value']:.4f}")
    print(f"   ✓ Gradient detected: {gradient_result['gradient_detected']}")
    
    # Test 2: Distance-Dependent Signaling
    print("\n2. Testing Distance-Dependent Signaling Range...")
    distance_validator = DistanceDependentValidation()
    distance_result = distance_validator.run('IL-6', noise_cv=0.3)
    results['distance'] = distance_result
    
    print(f"   Characteristic length: {distance_result['characteristic_length']:.1f} μm")
    for item in distance_result['I_by_distance']:
        print(f"   Distance {item['distance']:.0f} μm: I = {item['mean_I']:.3f}")
    
    # Test 3: Negative Controls
    print("\n3. Running Negative Controls...")
    control_validator = NegativeControlValidation()
    control_result = control_validator.run(noise_cv=0.3)
    results['controls'] = control_result
    
    for control_name, control_data in control_result['results'].items():
        if 'isotropic' in control_data:
            print(f"   {control_name}: Isotropic = {control_data['isotropic']}")
        elif 'weak_signal' in control_data:
            print(f"   {control_name}: Weak signal = {control_data['weak_signal']}")
    
    # Test 4: CXCL12 Chemokine Gradient
    print("\n4. Testing CXCL12 Chemokine Gradient...")
    chemokine_validator = ChemokineGradientValidation()
    chemokine_result = chemokine_validator.run(noise_cv=0.3)
    results['chemokine'] = chemokine_result
    
    print(f"   Expected direction: {chemokine_result['expected_direction']}°")
    print(f"   Detected direction: {chemokine_result['detected_direction']:.1f}°")
    print(f"   Angular error: {chemokine_result['angular_error']:.1f}°")
    print(f"   Peak I: {chemokine_result['peak_I']:.3f}")
    print(f"   p-value: {chemokine_result['p_value']:.4f}")
    print(f"   ✓ Direction correct: {chemokine_result['direction_correct']}")
    
    # Visualization
    print("\n5. Generating Visualizations...")
    
    # Visualize gradient scenario
    if 'concentration_field' in gradient_result:
        fig1 = visualize_cytokine_field(
            gradient_result['concentration_field'],
            gradient_result['positions'],
            gradient_result['expression'],
            "IL-2 Gradient Detection"
        )
        plt.savefig('gradient_validation.png', dpi=300, bbox_inches='tight')
    
    # Visualize chemokine scenario  
    if 'concentration_field' in chemokine_result:
        fig2 = visualize_cytokine_field(
            chemokine_result['concentration_field'],
            chemokine_result['positions'],
            chemokine_result['expression'],
            "CXCL12 Chemokine Gradient"
        )
        plt.savefig('chemokine_validation.png', dpi=300, bbox_inches='tight')
    
    # Power Analysis (reduced for demo)
    print("\n6. Running Power Analysis (this may take a few minutes)...")
    power_df = power_analysis(
        n_simulations=50,  # Reduced for demo
        fold_changes=[1.5, 2.0, 3.0],
        noise_levels=[0.2, 0.3, 0.5],
        n_cells=500  # Reduced for speed
    )
    
    print("\nPower Analysis Results:")
    print(power_df.to_string())
    
    # Plot power curves
    fig3 = plot_power_curves(power_df)
    plt.savefig('power_analysis.png', dpi=300, bbox_inches='tight')
    
    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    
    success_criteria = {
        'Gradient Detection': gradient_result['gradient_detected'],
        'Angular Accuracy': gradient_result['angular_error'] < 30,
        'Statistical Significance': gradient_result['p_value'] < 0.05,
        'Negative Controls Pass': control_result['results']['random_expression']['isotropic'],
        'Chemokine Direction': chemokine_result['direction_correct'],
        'Power > 80% at FC=2': any(
            (row['fold_change'] == 2.0 and row['noise_cv'] <= 0.3 and row['power'] >= 0.8)
            for _, row in power_df.iterrows()
        ) if len(power_df) > 0 else False
    }
    
    for criterion, passed in success_criteria.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{criterion:30s}: {status}")
    
    overall_success = all(success_criteria.values())
    print(f"\nOverall Validation: {'✓ SUCCESS' if overall_success else '✗ NEEDS IMPROVEMENT'}")
    
    return results, power_df

# =============================================================================
# Entry Point
# =============================================================================

if __name__ == "__main__":
    # Run complete validation
    results, power_df = run_complete_validation()
    
    # Save results
    import pickle
    with open('validation_results.pkl', 'wb') as f:
        pickle.dump({'results': results, 'power_analysis': power_df}, f)
    
    print("\nResults saved to validation_results.pkl")
    print("Figures saved as PNG files")
    
    # Show plots
    plt.show()

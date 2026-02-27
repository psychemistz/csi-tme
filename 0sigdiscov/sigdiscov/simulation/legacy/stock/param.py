#!/usr/bin/env python3
"""
Cytokine Parameter Reference with Literature Values
Comprehensive database for IL-2, TNF, TGF-β1, IFN-γ

Author: Seongyong Park
Institution: Cancer Data Science Lab, NCI, NIH
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List

# =============================================================================
# Detailed Cytokine Parameters with Literature References
# =============================================================================

@dataclass
class DetailedCytokineParams:
    """Extended cytokine parameters with biological context"""
    # Basic properties
    name: str
    molecular_weight: float  # kDa
    
    # Diffusion parameters
    diffusion_coefficient: float  # μm²/s
    diffusion_reference: str
    
    # Kinetics
    degradation_rate: float  # 1/s
    half_life_minutes: float  # minutes
    kinetics_reference: str
    
    # Production
    secretion_rate_per_cell: float  # molecules/cell/hour
    max_concentration_tissue: float  # nM
    production_reference: str
    
    # Receptor binding
    Kd: float  # nM (dissociation constant)
    receptor_density: float  # receptors/cell
    binding_reference: str
    
    # Signaling
    hill_coefficient: float
    EC50: float  # nM (half-maximal response)
    max_fold_change: float
    signaling_range: float  # μm (characteristic length)
    
    # Cell type specific production rates (relative)
    production_by_cell_type: Dict[str, float] = field(default_factory=dict)
    
    # Cell type specific consumption rates (relative)
    consumption_by_cell_type: Dict[str, float] = field(default_factory=dict)
    
    # Downstream genes
    primary_targets: List[str] = field(default_factory=list)
    secondary_targets: List[str] = field(default_factory=list)
    
    # Biological function
    primary_function: str = ""
    cell_types_affected: List[str] = field(default_factory=list)

# =============================================================================
# IL-2 (Interleukin-2)
# =============================================================================

IL2_PARAMS = DetailedCytokineParams(
    name="IL-2",
    molecular_weight=15.5,
    
    # Diffusion
    diffusion_coefficient=20.0,  # μm²/s
    diffusion_reference="Thurley et al., 2015, Cell Reports",
    
    # Kinetics
    degradation_rate=3.85e-4,  # 1/s (t½ = 30 min)
    half_life_minutes=30,
    kinetics_reference="Fallon & Lauffenburger, 2000, Cytokine",
    
    # Production
    secretion_rate_per_cell=1000,  # molecules/cell/hour
    max_concentration_tissue=1.0,  # nM
    production_reference="Busse et al., 2010, PNAS",
    
    # Receptor binding
    Kd=0.01,  # 10 pM (high affinity IL-2R)
    receptor_density=5000,  # receptors/cell
    binding_reference="Wang & Smith, 1987, J Exp Med",
    
    # Signaling
    hill_coefficient=1.5,
    EC50=0.02,  # nM
    max_fold_change=5.0,
    signaling_range=230,  # μm (√(D/η))
    
    # Cell type specific production
    production_by_cell_type={
        'CD4_T': 2.0,      # Major producer
        'CD8_T': 1.5,      # High producer
        'Treg': 0.5,       # Low producer
        'NK': 0.3,         # Minor producer
        'Macrophage': 0.0,
        'Fibroblast': 0.0,
        'B_cell': 0.1
    },
    
    # Cell type specific consumption
    consumption_by_cell_type={
        'CD4_T': 1.0,      # Auto/paracrine
        'CD8_T': 2.0,      # Major consumer
        'Treg': 3.0,       # Highest consumer (survival)
        'NK': 1.5,         # High consumer
        'Macrophage': 0.5,
        'Fibroblast': 0.0,
        'B_cell': 0.5
    },
    
    # Targets
    primary_targets=['IL2RA', 'IL2RB', 'IL2RG', 'STAT5A', 'STAT5B'],
    secondary_targets=['IFNG', 'GZMB', 'PRF1', 'TNF', 'BCL2'],
    
    # Function
    primary_function="T cell proliferation and survival",
    cell_types_affected=['CD4_T', 'CD8_T', 'Treg', 'NK']
)

# =============================================================================
# TNF-α (Tumor Necrosis Factor alpha)
# =============================================================================

TNF_PARAMS = DetailedCytokineParams(
    name="TNF-α",
    molecular_weight=17.3,  # Monomer, 51.9 for trimer
    
    # Diffusion
    diffusion_coefficient=40.0,  # μm²/s
    diffusion_reference="Chevalier et al., 2021, J Cell Sci",
    
    # Kinetics
    degradation_rate=5.78e-4,  # 1/s (t½ = 20 min)
    half_life_minutes=20,
    kinetics_reference="Oliver et al., 1993, Blood",
    
    # Production
    secretion_rate_per_cell=500,  # molecules/cell/hour
    max_concentration_tissue=0.5,  # nM
    production_reference="Xue et al., 2015, Immunity",
    
    # Receptor binding
    Kd=0.02,  # 20 pM (TNFR1)
    receptor_density=3000,  # receptors/cell
    binding_reference="Grell et al., 1995, Cell",
    
    # Signaling
    hill_coefficient=2.0,
    EC50=0.05,  # nM
    max_fold_change=10.0,
    signaling_range=260,  # μm
    
    # Cell type specific production
    production_by_cell_type={
        'Macrophage': 3.0,  # Major producer
        'DC': 2.0,          # High producer
        'CD4_T': 0.5,
        'CD8_T': 0.5,
        'NK': 1.0,
        'Neutrophil': 2.5,
        'Fibroblast': 0.2
    },
    
    # Cell type specific consumption
    consumption_by_cell_type={
        'Macrophage': 1.0,
        'DC': 1.0,
        'CD4_T': 1.5,
        'CD8_T': 1.5,
        'Fibroblast': 2.0,  # Major responder
        'Endothelial': 2.5,
        'Epithelial': 2.0
    },
    
    # Targets
    primary_targets=['NFKB1', 'RELA', 'JUN', 'FOS', 'TNFAIP3'],
    secondary_targets=['IL6', 'IL1B', 'CXCL8', 'ICAM1', 'VCAM1'],
    
    # Function
    primary_function="Inflammation and apoptosis",
    cell_types_affected=['Macrophage', 'Fibroblast', 'Endothelial', 'T_cells']
)

# =============================================================================
# TGF-β1 (Transforming Growth Factor beta 1)
# =============================================================================

TGFB1_PARAMS = DetailedCytokineParams(
    name="TGF-β1",
    molecular_weight=25.0,  # Dimer: 25 kDa
    
    # Diffusion (slower due to ECM binding)
    diffusion_coefficient=10.0,  # μm²/s
    diffusion_reference="Batlle & Massagué, 2019, Nature Reviews Cancer",
    
    # Kinetics
    degradation_rate=2.31e-4,  # 1/s (t½ = 50 min)
    half_life_minutes=50,
    kinetics_reference="Wakefield et al., 1990, J Cell Biol",
    
    # Production
    secretion_rate_per_cell=200,  # molecules/cell/hour
    max_concentration_tissue=2.0,  # nM
    production_reference="Flavell et al., 2010, Nature Reviews Immunology",
    
    # Receptor binding
    Kd=0.05,  # 50 pM (TGFβR2)
    receptor_density=10000,  # receptors/cell
    binding_reference="Massagué, 2012, Nature Reviews Mol Cell Biol",
    
    # Signaling
    hill_coefficient=1.8,
    EC50=0.1,  # nM
    max_fold_change=8.0,
    signaling_range=210,  # μm
    
    # Cell type specific production
    production_by_cell_type={
        'Treg': 2.5,        # Major producer
        'Macrophage': 1.5,
        'DC': 1.0,
        'Fibroblast': 2.0,  # Major producer
        'CAF': 3.0,         # Cancer-associated fibroblasts
        'Platelet': 2.0,
        'Epithelial': 0.5
    },
    
    # Cell type specific consumption
    consumption_by_cell_type={
        'CD4_T': 1.5,       # Suppressed by TGF-β
        'CD8_T': 1.5,       # Suppressed by TGF-β
        'Treg': 1.0,
        'Fibroblast': 2.0,  # Activation to myofibroblast
        'CAF': 2.5,
        'Epithelial': 2.0,  # EMT
        'Endothelial': 1.5
    },
    
    # Targets
    primary_targets=['SMAD2', 'SMAD3', 'SMAD4', 'SMAD7', 'TGFBR1'],
    secondary_targets=['SERPINE1', 'ACTA2', 'COL1A1', 'FN1', 'CTGF'],
    
    # Function
    primary_function="Immunosuppression and fibrosis",
    cell_types_affected=['T_cells', 'Fibroblasts', 'Epithelial', 'Endothelial']
)

# =============================================================================
# IFN-γ (Interferon gamma)
# =============================================================================

IFNG_PARAMS = DetailedCytokineParams(
    name="IFN-γ",
    molecular_weight=34.0,  # Homodimer
    
    # Diffusion
    diffusion_coefficient=30.0,  # μm²/s
    diffusion_reference="Thibaut et al., 2020, Cytokine",
    
    # Kinetics
    degradation_rate=2.89e-4,  # 1/s (t½ = 40 min)
    half_life_minutes=40,
    kinetics_reference="Pernis et al., 1995, Immunol Today",
    
    # Production
    secretion_rate_per_cell=800,  # molecules/cell/hour
    max_concentration_tissue=1.5,  # nM
    production_reference="Schoenborn & Wilson, 2007, Advances in Immunology",
    
    # Receptor binding
    Kd=0.1,  # 100 pM (IFNGR1/2)
    receptor_density=8000,  # receptors/cell
    binding_reference="Marsters et al., 1995, J Biol Chem",
    
    # Signaling
    hill_coefficient=1.6,
    EC50=0.2,  # nM
    max_fold_change=12.0,
    signaling_range=320,  # μm
    
    # Cell type specific production
    production_by_cell_type={
        'CD4_T': 2.0,       # Th1 cells
        'CD8_T': 3.0,       # Major producer
        'NK': 2.5,          # Major producer
        'NKT': 2.0,
        'Macrophage': 0.2,
        'DC': 0.1,
        'B_cell': 0.0
    },
    
    # Cell type specific consumption
    consumption_by_cell_type={
        'Macrophage': 3.0,  # Major responder (M1 activation)
        'DC': 2.5,          # Major responder
        'CD4_T': 0.5,
        'CD8_T': 0.5,
        'B_cell': 1.5,
        'Epithelial': 2.0,
        'Endothelial': 1.5
    },
    
    # Targets
    primary_targets=['STAT1', 'IRF1', 'JAK1', 'JAK2', 'SOCS1'],
    secondary_targets=['HLA-DRA', 'HLA-DRB1', 'CD274', 'CXCL9', 'CXCL10'],
    
    # Function
    primary_function="Macrophage activation and antigen presentation",
    cell_types_affected=['Macrophage', 'DC', 'B_cell', 'Endothelial']
)

# =============================================================================
# Cell Type Interaction Matrix
# =============================================================================

def create_interaction_matrix():
    """
    Create comprehensive interaction matrix for all cytokines and cell types
    """
    cytokines = ['IL-2', 'TNF-α', 'TGF-β1', 'IFN-γ']
    cell_types = ['CD4_T', 'CD8_T', 'Treg', 'NK', 'Macrophage', 
                  'DC', 'Fibroblast', 'B_cell']
    
    # Production matrix
    production_matrix = pd.DataFrame(index=cell_types, columns=cytokines)
    
    production_matrix.loc['CD4_T'] = [2.0, 0.5, 0.2, 2.0]
    production_matrix.loc['CD8_T'] = [1.5, 0.5, 0.1, 3.0]
    production_matrix.loc['Treg'] = [0.5, 0.0, 2.5, 0.0]
    production_matrix.loc['NK'] = [0.3, 1.0, 0.0, 2.5]
    production_matrix.loc['Macrophage'] = [0.0, 3.0, 1.5, 0.2]
    production_matrix.loc['DC'] = [0.0, 2.0, 1.0, 0.1]
    production_matrix.loc['Fibroblast'] = [0.0, 0.2, 2.0, 0.0]
    production_matrix.loc['B_cell'] = [0.1, 0.0, 0.0, 0.0]
    
    # Consumption matrix
    consumption_matrix = pd.DataFrame(index=cell_types, columns=cytokines)
    
    consumption_matrix.loc['CD4_T'] = [1.0, 1.5, 1.5, 0.5]
    consumption_matrix.loc['CD8_T'] = [2.0, 1.5, 1.5, 0.5]
    consumption_matrix.loc['Treg'] = [3.0, 0.5, 1.0, 0.2]
    consumption_matrix.loc['NK'] = [1.5, 0.5, 0.5, 1.0]
    consumption_matrix.loc['Macrophage'] = [0.5, 1.0, 1.0, 3.0]
    consumption_matrix.loc['DC'] = [0.3, 1.0, 0.5, 2.5]
    consumption_matrix.loc['Fibroblast'] = [0.0, 2.0, 2.0, 0.5]
    consumption_matrix.loc['B_cell'] = [0.5, 0.5, 0.5, 1.5]
    
    return production_matrix, consumption_matrix

# =============================================================================
# Spatial Signaling Properties
# =============================================================================

def calculate_signaling_properties():
    """
    Calculate key spatial signaling properties for each cytokine
    """
    properties = []
    
    for params in [IL2_PARAMS, TNF_PARAMS, TGFB1_PARAMS, IFNG_PARAMS]:
        # Characteristic length scale
        lambda_decay = np.sqrt(params.diffusion_coefficient / params.degradation_rate)
        
        # Time to reach steady state
        t_steady = 3 / params.degradation_rate  # ~3 time constants
        
        # Maximum signaling distance (99% decay)
        max_distance = lambda_decay * np.log(100)  # 99% reduction
        
        # Peclet number (advection vs diffusion, assuming v ~ 0.1 μm/s)
        v_typical = 0.1  # μm/s (typical interstitial flow)
        Pe = v_typical * lambda_decay / params.diffusion_coefficient
        
        properties.append({
            'Cytokine': params.name,
            'λ (μm)': lambda_decay,
            't_steady (min)': t_steady / 60,
            'Max range (μm)': max_distance,
            'Peclet number': Pe,
            'Diffusion-dominated': Pe < 1
        })
    
    return pd.DataFrame(properties)

# =============================================================================
# Export Functions
# =============================================================================

def export_parameters_to_csv():
    """
    Export all parameters to CSV files for reference
    """
    # Basic parameters
    basic_params = []
    for params in [IL2_PARAMS, TNF_PARAMS, TGFB1_PARAMS, IFNG_PARAMS]:
        basic_params.append({
            'Cytokine': params.name,
            'MW (kDa)': params.molecular_weight,
            'D (μm²/s)': params.diffusion_coefficient,
            'η (1/s)': params.degradation_rate,
            't½ (min)': params.half_life_minutes,
            'Kd (nM)': params.Kd,
            'Hill n': params.hill_coefficient,
            'EC50 (nM)': params.EC50,
            'Max fold change': params.max_fold_change,
            'λ (μm)': params.signaling_range
        })
    
    df_basic = pd.DataFrame(basic_params)
    df_basic.to_csv('cytokine_parameters.csv', index=False)
    
    # Interaction matrices
    prod_matrix, cons_matrix = create_interaction_matrix()
    prod_matrix.to_csv('production_matrix.csv')
    cons_matrix.to_csv('consumption_matrix.csv')
    
    # Spatial properties
    spatial_props = calculate_signaling_properties()
    spatial_props.to_csv('spatial_properties.csv', index=False)
    
    print("Parameters exported to CSV files:")
    print("  - cytokine_parameters.csv")
    print("  - production_matrix.csv")
    print("  - consumption_matrix.csv")
    print("  - spatial_properties.csv")

# =============================================================================
# Display Functions
# =============================================================================

def display_all_parameters():
    """
    Display all cytokine parameters in formatted tables
    """
    print("="*80)
    print("COMPREHENSIVE CYTOKINE PARAMETERS")
    print("="*80)
    
    for params in [IL2_PARAMS, TNF_PARAMS, TGFB1_PARAMS, IFNG_PARAMS]:
        print(f"\n{params.name}")
        print("-"*40)
        print(f"Molecular weight: {params.molecular_weight} kDa")
        print(f"Diffusion coefficient: {params.diffusion_coefficient} μm²/s")
        print(f"Half-life: {params.half_life_minutes} minutes")
        print(f"Kd: {params.Kd*1000} pM")
        print(f"Signaling range: {params.signaling_range:.0f} μm")
        print(f"Primary function: {params.primary_function}")
        print(f"Reference: {params.diffusion_reference}")
    
    print("\n" + "="*80)
    print("SPATIAL SIGNALING PROPERTIES")
    print("="*80)
    
    spatial_props = calculate_signaling_properties()
    print(spatial_props.to_string(index=False))
    
    print("\n" + "="*80)
    print("CELL TYPE INTERACTION SUMMARY")
    print("="*80)
    
    prod_matrix, cons_matrix = create_interaction_matrix()
    
    print("\nProduction Matrix (relative rates):")
    print(prod_matrix)
    
    print("\nConsumption Matrix (relative rates):")
    print(cons_matrix)

# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("\nCytokine Parameter Reference Database")
    print("=====================================\n")
    
    # Display all parameters
    display_all_parameters()
    
    # Export to CSV
    print("\n" + "="*80)
    export_parameters_to_csv()
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print("\nThis reference provides:")
    print("1. Literature-based diffusion coefficients")
    print("2. Experimentally measured half-lives")
    print("3. Receptor binding affinities (Kd values)")
    print("4. Cell type-specific production/consumption rates")
    print("5. Calculated spatial signaling ranges")
    print("\nUse these parameters for accurate cytokine signaling simulations.")

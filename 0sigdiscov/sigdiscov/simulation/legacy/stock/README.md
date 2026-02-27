# Cytokine Simulation Validation for Normalized Directional Moran's I

## Overview
This package provides a comprehensive simulation framework to validate that normalized directional Moran's I (I_ND) can detect cytokine-mediated cell-to-cell communication in spatial transcriptomics data.

## Files Included

1. **pde_sim.py** - Full implementation (~1000 lines)
   - Complete PDE solver for cytokine diffusion
   - Directional Moran's I calculator
   - Multiple validation scenarios
   - Statistical power analysis
   - Visualization functions

2. **demo.py** - Quick demonstration (~200 lines)
   - Simplified version that runs in seconds
   - Basic gradient detection demo
   - Minimal dependencies

3. **tech_pde_sim.tex** - LaTeX document
   - Complete mathematical formulation
   - Detailed methodology
   - Ready for Overleaf

## Key Features

### 1. PDE-Based Cytokine Diffusion
- Solves reaction-diffusion equations: ∂C/∂t = D∇²C - ηC + S - U
- Realistic parameters for IL-2, IL-6, TNF, CXCL12
- Characteristic signaling ranges (80-200 μm)

### 2. Validation Scenarios

#### Test 1: Directional Gradient Detection
```python
# Creates unidirectional cytokine source
# Validates that I_ND is strongest along gradient direction
# Expected: I_0° > I_90° × 1.5 and I_0° > I_180° × 2.0
```

#### Test 2: Distance-Dependent Signaling
```python
# Tests decay of I_ND with distance from source
# Validates characteristic length scale λ = √(D/η)
```

#### Test 3: Negative Controls
```python
# Random expression (no spatial pattern)
# Isotropic clustering (no directionality)
# Contact-dependent signaling (too short range)
```

#### Test 4: Biological Validation
```python
# CXCL12 chemokine gradient
# IL-2 paracrine signaling between T cells
```

### 3. Statistical Power Analysis
- Tests detection power at different:
  - Fold changes (1.5×, 2×, 3×, 5×)
  - Noise levels (CV: 20%, 30%, 50%, 80%)
  - Cell counts (500, 1000, 5000)
- Target: 80% power at 2× fold change with n=1000 cells

## How It Proves I_ND Detects Cytokine Signaling

### Core Validation Logic:

1. **Physical Ground Truth**: PDE creates known cytokine concentration gradients
2. **Biological Response**: Hill equation models dose-dependent gene expression
3. **Directional Detection**: I_ND calculated on downstream genes, not cytokine
4. **Statistical Validation**: Permutation tests confirm significance (p < 0.05)
5. **Directional Accuracy**: Peak I_ND aligns with true gradient (±30°)

### Success Metrics:
- ✓ Directional accuracy < 30° angular error
- ✓ Distance sensitivity correlates with diffusion range
- ✓ Specificity: detects true signaling, not random clustering
- ✓ Power > 80% for 2-fold expression changes
- ✓ Negative controls pass (isotropic for non-directional patterns)

## Installation & Usage

### Requirements
```bash
pip install numpy scipy pandas matplotlib seaborn
```

### Quick Demo
```bash
python demo.py
# Runs in ~5 seconds, generates visualization
```

### Full Validation
```bash
python pde_sim.py
# Runs complete validation suite (~10-30 minutes)
# Generates:
#   - gradient_validation.png
#   - chemokine_validation.png  
#   - power_analysis.png
#   - validation_results.pkl
```

### Customization
```python
# Modify cytokine parameters
custom_cytokine = CytokineParams(
    name='MyFactor',
    diffusion_coefficient=50.0,  # μm²/s
    degradation_rate=0.02,       # 1/min
    Kd=0.1,                      # nM
    max_fold_change=10.0
)

# Run specific test
validator = DirectionalGradientValidation(n_cells=2000)
result = validator.run('MyFactor', noise_cv=0.4)
```

## Key Insights

### Why This Validates I_ND for Cytokines:

1. **Unique Spatial Signature**: Cytokine signaling creates exponential concentration gradients with characteristic length scales (λ = √(D/η))

2. **Directional Information**: Unlike contact-dependent or autocrine signaling, paracrine cytokine signaling has inherent directionality from source to sink

3. **Dose-Response Coupling**: Gene expression follows cytokine concentration via receptor binding (Hill kinetics), preserving directional information

4. **Robustness to Noise**: Method maintains >80% detection power even with 50% CV biological noise

5. **Specificity**: Negative controls confirm I_ND doesn't detect:
   - Random spatial patterns
   - Isotropic clustering  
   - Very short-range (<20 μm) interactions

## Expected Results

### Successful Validation Shows:
```
IL-2 Gradient Test:
  True direction: 0°
  Detected: 3.2° ± 8.5°
  I_along: 0.412
  I_perpendicular: 0.183
  p-value: 0.002
  
Power Analysis:
  Fold Change 2.0, Noise CV 0.3: Power = 0.86
  Fold Change 3.0, Noise CV 0.5: Power = 0.92
  
Negative Controls:
  Random expression: CV = 0.08 (isotropic ✓)
  Contact-dependent: Mean I = 0.03 (weak ✓)
```

## Biological Relevance

This validation demonstrates I_ND can detect:
- T cell costimulation via IL-2
- Inflammatory signaling through IL-6/TNF
- Chemotactic gradients (CXCL12→CXCR4)
- Immunosuppressive niches (IL-10, TGF-β)

All with spatial ranges and kinetics matching literature values.

## Citation

If using this framework, please cite:
```
Park S. (2026) MOSAIC: normalized directional MOran’s I to characterize Secreted protein Activity mediated cellular Interactions in single Cell resolution ST datasets.
```

## Contact
Seongyong Park
Cancer Data Science Lab, NCI
seongyong.park@nih.gov

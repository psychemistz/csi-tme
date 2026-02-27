# I_ND Simulation Framework

Simulation framework for validating the **Normalized Directional Moran's I (I_ND)** spatial statistic in controlled reaction-diffusion scenarios of cell-cell signaling.

## Overview

This simulation models intercellular communication in a 2D tissue:

1. **Sender cells** secrete a diffusible factor (e.g., cytokine, chemokine)
2. The factor diffuses through the tissue and is consumed by **receiver cells**
3. Receivers upregulate a **response gene** proportional to local factor concentration
4. **I_ND** measures the spatial correlation between sender factor expression and receiver response expression as a function of distance

The key question: does I_ND peak near the biophysical signaling range (lambda)?

## Biophysical Model

The concentration field follows a **Michaelis-Menten steady-state diffusion** model. For a point source secreting factor at rate Q, the 2D steady-state solution is:

```
C(r) = Q * exp(-r / lambda) / sqrt(r)
```

where the **characteristic decay length** (signaling range) is:

```
lambda = sqrt(D * Kd / (n_eff * k_max))
```

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| Diffusion coefficient | D | um^2/s |
| Dissociation constant | Kd | nM |
| Maximum uptake rate | k_max | molecules/cell/s |
| Effective receiver density | n_eff | cells/um^2 (= n_receivers * p_r) |

For multiple source positions, concentration fields are superimposed linearly.

## Simulation Pipeline

Each run sweeps over a range of **receiver fractions** (proportion of total cells that are receivers). For each fraction:

1. **Generate positions** -- Scatter `n_total` cells uniformly in a circular domain (radius `max_radius`)
2. **Assign cell types** -- Designate sender cells (active + silent) at configured positions; randomly sample receiver cells from the remaining population
3. **Factor expression** -- Generate sender gene expression using the configured factor model (deterministic or stochastic)
4. **Diffusion** -- Solve the steady-state concentration field at all cell positions given active sender locations and expression levels
5. **Response expression** -- Generate receiver gene expression using the configured response model, dependent on local concentration
6. **I_ND computation** -- Sweep distance parameter d across the domain, computing I_ND at each distance using ring or Gaussian annular kernels

## Expression Models

### Factor Models

| Model | Description | Key Parameters |
|-------|-------------|----------------|
| `deterministic` | Active senders get `HIGH_FACTOR`, others get `BASAL_FACTOR` | `HIGH_FACTOR`, `BASAL_FACTOR` |
| `stochastic` | Active senders switch on/off with probability `p_express`; on-cells drawn from gamma or lognormal | `p_express`, `expr_cv`, `use_gamma_dist` |
| `stochastic_ref` | On/off Bernoulli switch; on-cells: `F_high * LN(0, sigma_f^2)`; off-cells: `F_basal * LN(0, sigma_f_basal^2)` | `p_express`, `sigma_f`, `sigma_f_basal` |
| `bernoulli_mixture` | `F_i = S_i * F_high * LN(0, sigma_f^2) + (1-S_i) * F_basal * LN(0, sigma_f_b^2)` where `S_i ~ Bernoulli(p_s)` | `p_s`, `sigma_f`, `sigma_f_b` |

### Response Models

| Model | Description | Key Parameters |
|-------|-------------|----------------|
| `deterministic` | `R = R_basal * (1 + FC * Act)` where `Act = C/(Kd+C)` | `BASAL_RESPONSIVE`, `FOLD_CHANGE` |
| `stochastic_hill` | Bernoulli switch with Hill-function probability; responding cells get fold-change activation with lognormal noise | `p_respond_max`, `response_hill_coef`, `sigma_r` |
| `bernoulli_constant` | `B_j ~ Bernoulli(p_r)` constant probability; mixture of activated + basal with lognormal noise | `p_r`, `sigma_r`, `sigma_r_b` |
| `bernoulli_hill` | `B_j ~ Bernoulli(p_r_max * C^n / (K_p^n + C^n))` concentration-dependent probability | `p_r_max`, `K_p`, `hill_n`, `sigma_r`, `sigma_r_b` |

## Sender Configurations

| Config | Description |
|--------|-------------|
| `center` | All active senders placed at domain center |
| `center_silent_distributed` | Active senders at center; silent senders scattered throughout domain |
| `multi_fixed` | Active senders distributed across 5 fixed positions (C, W, E, N, S) at `position_offset` from center |
| `multi_random` | Active senders distributed across `n_active_positions` random positions with `min_separation` constraint |

## I_ND Methods

### Ring Kernel

A binary annular weight matrix. For each sender-receiver pair at distance d:

```
W(d) = 1   if  (distance - BW/2) < d <= (distance + BW/2)
W(d) = 0   otherwise
```

Row-normalized, then I_ND = cosine similarity between z-scored sender expression and weighted mean neighbor z-scored receiver expression.

### Gaussian Annular Kernel

Gaussian weights centered at sender (d=0) with `sigma = outer_radius / sigma_fraction`, masked to the annular region `[outer_radius - bandwidth, outer_radius]`. Otherwise computed identically to the ring kernel.

## Variance-Stabilizing Transforms (VST)

Applied before I_ND computation when `vst_method` is set. Optional zero-inflation (dropout) can be applied first.

| Method | Transform | Notes |
|--------|-----------|-------|
| `log1p` | `log(1+x)` centered by global mean | Standard scRNA-seq normalization |
| `pearson` | Pearson residuals from NB model | sctransform-style; `theta` dispersion parameter |
| `shifted_log` | `log(x+1) - log(1) - 0.5` | Ensures non-expressors map to negative values |

## `core.py` Reference

Public functions in `core.py`:

| Function | Purpose |
|----------|---------|
| `calculate_lambda` | Compute signaling range from biophysical parameters |
| `generate_circular_positions` | Uniform random positions in circular domain |
| `solve_concentration_field_MM` | Single-source concentration field (co-located senders) |
| `solve_concentration_field_MM_multipos` | Multi-source concentration via superposition |
| `compute_IND_ring` | I_ND with ring (annular) kernel |
| `compute_IND_gaussian_annular` | I_ND with Gaussian annular kernel |
| `get_5_positions` | Generate 5 cardinal positions (C, W, E, N, S) |
| `get_random_positions` | Generate N random separated positions |
| `distribute_active_senders` | Assign active senders across multiple positions |
| `generate_stochastic_expression` | Gamma/lognormal on/off stochastic factor model |
| `generate_stochastic_expression_ref` | Reference lognormal stochastic factor model |
| `generate_bernoulli_factor_expression` | Bernoulli mixture factor model |
| `hill_function` | Hill function for concentration-dependent probability |
| `generate_bernoulli_response_constant` | Constant-probability Bernoulli response model |
| `generate_bernoulli_response_hill` | Hill-function Bernoulli response model |
| `generate_stochastic_response` | Hybrid Bernoulli-Hill stochastic response model |
| `apply_vst_log1p` | Log1p VST with centering |
| `apply_vst_pearson` | Pearson residual VST |
| `apply_vst_shifted_log` | Shifted log VST |

## `unified_sim.py` Usage

```bash
# List all available presets
python unified_sim.py --list

# Run a preset with default settings
python unified_sim.py demo

# Run with a specific random seed and output directory
python unified_sim.py demo3b --seed 42 --output ./my_results
```

### Key Config Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_total` | 100000 | Total cells in the domain |
| `max_radius` | 5000 | Domain radius (um) |
| `n_active` | 20 | Number of active sender cells |
| `n_silent` | 0 | Number of silent sender cells |
| `receiver_fractions` | [0.02..0.80] | Receiver fractions to sweep |
| `BANDWIDTH` | 100 | Ring/annular kernel width |
| `test_distance_end` | 5000 | Max distance for I_ND sweep |
| `test_distance_step` | 25 | Distance step for I_ND sweep |

## Preset Reference

| Preset | Factor Model | Response Model | Sender Config | Key Differences |
|--------|-------------|----------------|---------------|-----------------|
| `demo` | deterministic | deterministic | center | Baseline: simple center source |
| `demo1` | deterministic | deterministic | center | Same as demo (alias) |
| `demo1b` | deterministic | deterministic | center | + Gaussian annular kernel comparison; fewer fractions |
| `demo1_stc` | bernoulli_mixture | bernoulli_constant | center | Stochastic factor + constant-prob stochastic response |
| `demo1_stc2` | bernoulli_mixture | bernoulli_hill | center | Stochastic factor + Hill-dependent stochastic response |
| `demo2a` | deterministic | deterministic | center_silent_distributed | 100 silent senders; HIGH_FACTOR=100, FC=50 |
| `demo2b` | deterministic | deterministic | center_silent_distributed | 20 silent senders; 80% silent receivers |
| `demo2c` | deterministic | deterministic | multi_random | 20 random positions; secretion_rate=100 |
| `demo3` | stochastic | deterministic | multi_random | 200 active, 100 silent; 1 random cluster; BW=20 |
| `demo3b` | stochastic_ref | stochastic_hill | multi_random | 200 active; stochastic both sides; Kd=1, FC=10 |
| `demo_det` | deterministic | deterministic | center | Biophysical params: D=1060, k_max=300, Kd=5; center=[2500,2500] |
| `demo_det_dec` | deterministic | deterministic | center_silent_distributed | Like demo_det + 200 silent senders distributed |
| `demo_vst` | deterministic | deterministic | center | VST normalization (log1p); zero-inflation 70%/50% |

## Legacy Files

The `legacy/` directory contains the original individual demo scripts that have been consolidated into `unified_sim.py`. Each preset in `unified_sim.py` reproduces the corresponding demo script's behavior exactly (pixel-identical output verified).

`legacy/stock/` contains even older prototype code (PDE solvers, 3-cell-type models, density tests) from earlier development iterations.

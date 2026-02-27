#!/usr/bin/env python3
"""
Unified Simulation Script for I_ND Validation Experiments.

Replaces 13 individual demo files with a single config-driven script.
Each demo is available as a named preset.

Usage:
    python unified_sim.py <preset> [--seed N] [--output DIR]
    python unified_sim.py --list          # show available presets

Examples:
    python unified_sim.py demo
    python unified_sim.py demo3b --seed 123 --output ./results
"""

import argparse
import os
from dataclasses import dataclass, field
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter1d

from core import (
    apply_vst_log1p,
    apply_vst_pearson,
    apply_vst_shifted_log,
    calculate_lambda,
    compute_IND_gaussian_annular,
    compute_IND_ring,
    distribute_active_senders,
    generate_bernoulli_factor_expression,
    generate_bernoulli_response_constant,
    generate_bernoulli_response_hill,
    generate_stochastic_expression,
    generate_stochastic_expression_ref,
    generate_stochastic_response,
    get_5_positions,
    get_random_positions,
    solve_concentration_field_MM,
    solve_concentration_field_MM_multipos,
)

white_red = LinearSegmentedColormap.from_list("white_red", ["white", "red"])


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class ExperimentConfig:
    name: str = "experiment"

    # Domain
    center: List[float] = field(default_factory=lambda: [0.0, 0.0])
    n_total: int = 100000
    max_radius: float = 5000.0

    # Senders
    n_active: int = 20
    n_silent: int = 0
    silent_expr_zero: bool = False  # True = silent senders get 0 expression

    # Receivers
    receiver_fractions: List[float] = field(
        default_factory=lambda: [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
    )
    receiver_silent_fraction: float = 0.0  # fraction of receivers that are silent
    ind_uses_active_receivers_only: bool = False

    # Biophysics
    k_max: float = 10.0
    Kd: float = 30.0
    D: float = 100.0

    # Gene expression parameters
    BASAL_FACTOR: float = 0.1
    HIGH_FACTOR: float = 10.0
    BASAL_RESPONSIVE: float = 0.1
    FOLD_CHANGE: float = 2.0
    sigma_f: float = 0.1
    sigma_r: float = 0.1

    # Factor expression model: deterministic | stochastic | stochastic_ref | bernoulli_mixture
    factor_model: str = "deterministic"

    # Bernoulli mixture params (factor_model='bernoulli_mixture')
    p_s: float = 1.0
    sigma_f_b: float = 0.1

    # Stochastic on/off params (factor_model='stochastic')
    p_express: float = 0.7
    expr_cv: float = 0.5
    use_gamma_dist: bool = True

    # Stochastic ref params (factor_model='stochastic_ref')
    sigma_f_basal: float = 0.1

    # Response model: deterministic | stochastic_hill | bernoulli_constant | bernoulli_hill
    response_model: str = "deterministic"

    # Bernoulli constant response params
    p_r: float = 1.0
    sigma_r_b: float = 0.1

    # Bernoulli Hill response params
    p_r_max: float = 1.0
    K_p: float = 10.0
    hill_n: float = 1.0

    # Stochastic Hill response params (response_model='stochastic_hill')
    p_respond_max: float = 0.9
    response_hill_coef: float = 1.0

    # Sender positions: center | center_silent_distributed | multi_fixed | multi_random
    sender_positions: str = "center"
    n_active_positions: int = 5
    position_offset: float = 3000.0
    min_separation: float = 500.0
    secretion_rate: float = 1.0
    active_threshold: float = 1.0

    # Whether to fix sender assignment across fractions
    fix_senders_across_fractions: bool = False

    # I_ND methods
    ind_methods: List[str] = field(default_factory=lambda: ["ring"])
    BANDWIDTH: float = 100.0
    sigma_fraction: float = 3.0  # for gaussian_annular
    test_distance_start: Optional[float] = None  # defaults to BANDWIDTH/2
    test_distance_end: float = 5000.0
    test_distance_step: float = 25.0

    # VST
    vst_method: Optional[str] = None  # None | 'log1p' | 'pearson' | 'shifted_log'
    zero_inflate_factor: float = 0.0
    zero_inflate_responsive: float = 0.0
    vst_active_threshold: float = 1.0

    # Concentration plot normalization
    vmin_conc: float = 1e-4
    vmax_conc: float = 1e4


# =============================================================================
# Preset Configurations
# =============================================================================

def _base_config(**overrides):
    cfg = ExperimentConfig(**overrides)
    return cfg


PRESETS = {}


def _register(name, **kwargs):
    kwargs['name'] = name
    PRESETS[name] = kwargs


# demo.py
_register('demo',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=10.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=0,
    sender_positions='center', factor_model='deterministic', response_model='deterministic',
    ind_methods=['ring'],
)

# demo1.py
_register('demo1',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=10.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=0,
    sender_positions='center', factor_model='deterministic', response_model='deterministic',
    ind_methods=['ring'],
)

# demo1b.py
_register('demo1b',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=10.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=0,
    sender_positions='center', factor_model='deterministic', response_model='deterministic',
    ind_methods=['ring', 'gaussian_annular'], sigma_fraction=3.0,
    receiver_fractions=[0.02, 0.05, 0.10, 0.20, 0.30, 0.50],
    test_distance_start=110.0, test_distance_end=3000.0,
)

# demo1_stc.py
_register('demo1_stc',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=10.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    BANDWIDTH=100.0, n_active=20, n_silent=0,
    sender_positions='center',
    factor_model='bernoulli_mixture',
    p_s=1.0, sigma_f=0.1, sigma_f_b=0.1,
    response_model='bernoulli_constant',
    p_r=1.0, sigma_r=0.01, sigma_r_b=0.1,
    ind_methods=['ring'],
)

# demo1_stc2.py
_register('demo1_stc2',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=10.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    BANDWIDTH=100.0, n_active=20, n_silent=0,
    sender_positions='center',
    factor_model='bernoulli_mixture',
    p_s=1.0, sigma_f=0.1, sigma_f_b=0.1,
    response_model='bernoulli_hill',
    p_r_max=1.0, K_p=1.0, hill_n=1.0, sigma_r=0.01, sigma_r_b=0.1,
    ind_methods=['ring'],
)

# demo2a.py
_register('demo2a',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=100.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=50.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=100,
    sender_positions='center_silent_distributed',
    factor_model='deterministic', response_model='deterministic',
    ind_methods=['ring'],
)

# demo2b.py
_register('demo2b',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=10.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=20,
    sender_positions='center_silent_distributed',
    factor_model='deterministic', response_model='deterministic',
    receiver_silent_fraction=0.8,
    ind_uses_active_receivers_only=True,
    ind_methods=['ring'],
)

# demo2c.py
_register('demo2c',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=1.5, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=40,
    sender_positions='multi_random',
    n_active_positions=20, min_separation=800.0, position_offset=3000.0,
    secretion_rate=100.0,
    factor_model='deterministic', response_model='deterministic',
    fix_senders_across_fractions=True,
    ind_methods=['ring'],
)

# demo3.py
_register('demo3',
    center=[0.0, 0.0], k_max=50.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=1.5, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=20.0,
    n_active=200, n_silent=100,
    sender_positions='multi_random',
    n_active_positions=1, min_separation=0.0, position_offset=3000.0,
    secretion_rate=0.5, active_threshold=0.2,
    factor_model='stochastic',
    p_express=0.5, expr_cv=0.5, use_gamma_dist=False,
    response_model='deterministic',
    fix_senders_across_fractions=True,
    ind_methods=['ring'],
)

# demo3b.py
_register('demo3b',
    center=[0.0, 0.0], k_max=10.0, Kd=1.0, D=100.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=10.0, BASAL_RESPONSIVE=0.1, FOLD_CHANGE=10.0,
    BANDWIDTH=20.0,
    sigma_f=0.1, sigma_r=0.1,
    sigma_f_basal=0.1, sigma_r_b=0.1,
    n_active=200, n_silent=100,
    sender_positions='multi_random',
    n_active_positions=1, min_separation=0.0, position_offset=3000.0,
    secretion_rate=1.0, active_threshold=0.2,
    factor_model='stochastic_ref',
    p_express=0.9,
    response_model='stochastic_hill',
    p_respond_max=1.0, response_hill_coef=1.0,
    fix_senders_across_fractions=True,
    ind_methods=['ring'],
)

# demo_det.py
_register('demo_det',
    center=[2500.0, 2500.0], k_max=300.0, Kd=5.0, D=1060.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=2.0, BASAL_RESPONSIVE=2.0, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=20, silent_expr_zero=True,
    sender_positions='center',
    factor_model='deterministic', response_model='deterministic',
    receiver_fractions=[0.10, 0.20, 0.30, 0.40, 0.50],
    ind_methods=['ring'],
)

# demo_det_dec.py
_register('demo_det_dec',
    center=[2500.0, 2500.0], k_max=300.0, Kd=5.0, D=1060.0,
    BASAL_FACTOR=0.1, HIGH_FACTOR=2.0, BASAL_RESPONSIVE=2.0, FOLD_CHANGE=2.0,
    sigma_f=0.1, sigma_r=0.1, BANDWIDTH=100.0,
    n_active=20, n_silent=200, silent_expr_zero=True,
    sender_positions='center_silent_distributed',
    factor_model='deterministic', response_model='deterministic',
    receiver_fractions=[0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80],
    ind_methods=['ring'],
)

# demo_vst.py
_register('demo_vst',
    center=[0.0, 0.0], k_max=10.0, Kd=30.0, D=100.0,
    BASAL_FACTOR=0.5, HIGH_FACTOR=50.0, BASAL_RESPONSIVE=0.5, FOLD_CHANGE=5.0,
    sigma_f=0.8, sigma_r=0.8, BANDWIDTH=100.0,
    n_active=20, n_silent=0,
    sender_positions='center',
    factor_model='deterministic', response_model='deterministic',
    vst_method='log1p',
    zero_inflate_factor=0.7, zero_inflate_responsive=0.5,
    vst_active_threshold=5.0,
    ind_methods=['ring'],
)


def make_config(preset_name):
    """Create an ExperimentConfig from a preset name."""
    if preset_name not in PRESETS:
        raise ValueError(f"Unknown preset: {preset_name}. Available: {list(PRESETS.keys())}")
    return ExperimentConfig(**PRESETS[preset_name])


# =============================================================================
# Simulation State
# =============================================================================

@dataclass
class SimState:
    """Bundles all state for a single receiver fraction iteration."""
    cfg: ExperimentConfig
    frac: float
    frac_str: str
    all_positions: np.ndarray
    sender_indices: np.ndarray
    active_indices: np.ndarray
    silent_indices: np.ndarray
    receiver_indices: np.ndarray
    factor_expr: np.ndarray
    responsive_expr: np.ndarray
    concentrations: np.ndarray
    lambda_val: float
    ind_curves: dict  # method_name -> list of dicts
    # Optional fields
    expressing_mask: Optional[np.ndarray] = None
    responding_mask: Optional[np.ndarray] = None
    response_probs: Optional[np.ndarray] = None
    mean_p_r: Optional[float] = None
    active_receiver_indices: Optional[np.ndarray] = None
    silent_receiver_indices: Optional[np.ndarray] = None
    factor_raw: Optional[np.ndarray] = None
    factor_vst: Optional[np.ndarray] = None
    responsive_raw: Optional[np.ndarray] = None
    responsive_vst: Optional[np.ndarray] = None
    position_dict: Optional[dict] = None
    active_assignments: Optional[dict] = None
    n_expressing: Optional[int] = None
    n_responding: Optional[int] = None


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_ind_summary(cfg, all_results, output_folder):
    """Plot I_ND summary curves for all receiver fractions."""
    methods = cfg.ind_methods

    if len(methods) == 1 and methods[0] == 'ring':
        _plot_ind_single(cfg, all_results, 'ring', output_folder)
    elif len(methods) == 2 and 'ring' in methods and 'gaussian_annular' in methods:
        _plot_ind_comparison(cfg, all_results, output_folder)
    else:
        for method in methods:
            _plot_ind_single(cfg, all_results, method, output_folder)


def _plot_ind_single(cfg, all_results, method, output_folder):
    """Plot I_ND curves for a single method."""
    plt.figure(figsize=(10, 8))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(cfg.receiver_fractions)))

    for f, col in zip(cfg.receiver_fractions, colors):
        state = all_results[f]
        curve = state.ind_curves[method]
        dists = np.array([x['distance'] for x in curve])
        vals = np.array([x['I_ND'] for x in curve])
        lambda_v = state.lambda_val

        vals_smooth = gaussian_filter1d(vals, sigma=1.5)

        label_parts = [f'{f*100:.0f}% Receivers (λ={lambda_v:.0f}μm)']
        if state.n_expressing is not None:
            label_parts[0] = f'{f*100:.0f}% Recv (λ={lambda_v:.0f}μm, n_on={state.n_expressing})'
        if state.n_responding is not None and cfg.response_model == 'stochastic_hill':
            label_parts[0] += f', resp={state.n_responding}'
        if state.mean_p_r is not None and cfg.response_model == 'bernoulli_hill':
            label_parts[0] = f'{f*100:.0f}% Recv (λ={lambda_v:.0f}μm, <p_r>={state.mean_p_r:.2f})'

        plt.plot(dists, vals_smooth, '-', color=col, linewidth=3, label=label_parts[0])
        plt.axvline(lambda_v, color=col, linestyle=':', linewidth=1.5, alpha=0.6)

    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.xlabel("Distance from Senders (μm)", fontsize=14, fontweight='bold')
    method_label = method.replace('_', ' ').title().replace(' ', ', ') if '_' in method else method.capitalize()
    plt.ylabel(f"I_ND ({method_label}, bw={cfg.BANDWIDTH:.0f}μm)", fontsize=14, fontweight='bold')
    plt.title(f"Diffusion Coefficient (D = {cfg.D:.0f} μm²/s)", fontsize=16, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.ylim(-1.1, 1.1)
    plt.xlim(0, cfg.test_distance_end)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"ind_summary_{method}.png"), dpi=300)
    plt.close()


def _plot_ind_comparison(cfg, all_results, output_folder):
    """Plot comparison between gaussian_annular and ring methods (demo1b style)."""
    fracs = cfg.receiver_fractions
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(fracs)))

    # Side-by-side comparison
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for method_idx, (method, ax) in enumerate(zip(['gaussian_annular', 'ring'], axes)):
        for f, col in zip(fracs, colors):
            state = all_results[f]
            curve = state.ind_curves[method]
            dists = np.array([x['distance'] for x in curve])
            vals = np.array([x['I_ND'] for x in curve])
            vals_smooth = gaussian_filter1d(vals, sigma=1.5)
            lambda_v = state.lambda_val
            ax.plot(dists, vals_smooth, '-', color=col, linewidth=2.5,
                    label=f'{f*100:.0f}% (λ={lambda_v:.0f}μm)')
            ax.axvline(lambda_v, color=col, linestyle=':', linewidth=1, alpha=0.5)
        ax.axhline(0, color='black', linestyle='-', alpha=0.3)
        ax.set_xlabel("Outer Radius (μm)", fontsize=13)
        ax.set_ylabel("I_ND", fontsize=13)
        title = "Gaussian ANNULAR" if method == 'gaussian_annular' else "Binary Ring"
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.2, 1.1)
        ax.set_xlim(0, cfg.test_distance_end)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, "ind_comparison.png"), dpi=300)
    plt.close()

    # Overlay
    plt.figure(figsize=(12, 8))
    for f, col in zip(fracs, colors):
        state = all_results[f]
        for method, ls in [('gaussian_annular', '-'), ('ring', '--')]:
            curve = state.ind_curves[method]
            dists = np.array([x['distance'] for x in curve])
            vals = np.array([x['I_ND'] for x in curve])
            vals_smooth = gaussian_filter1d(vals, sigma=1.5)
            label = f'{f*100:.0f}% Gaussian' if ls == '-' else None
            plt.plot(dists, vals_smooth, ls, color=col, linewidth=2.5,
                     label=label, alpha=0.7 if ls == '--' else 1.0)
        plt.axvline(state.lambda_val, color=col, linestyle=':', linewidth=1, alpha=0.4)
    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.xlabel("Outer Radius (μm)", fontsize=14, fontweight='bold')
    plt.ylabel("I_ND", fontsize=14, fontweight='bold')
    plt.title("Gaussian ANNULAR (solid) vs Binary Ring (dashed)", fontsize=16, fontweight='bold')
    plt.legend(fontsize=10, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.2, 1.1)
    plt.xlim(0, cfg.test_distance_end)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, "ind_overlay.png"), dpi=300)
    plt.close()

    # Difference
    plt.figure(figsize=(12, 6))
    for f, col in zip(fracs, colors):
        state = all_results[f]
        g_curve = state.ind_curves['gaussian_annular']
        r_curve = state.ind_curves['ring']
        dists = np.array([x['distance'] for x in g_curve])
        diff = np.array([g['I_ND'] - r['I_ND'] for g, r in zip(g_curve, r_curve)])
        diff_smooth = gaussian_filter1d(diff, sigma=1.5)
        plt.plot(dists, diff_smooth, '-', color=col, linewidth=2, label=f'{f*100:.0f}%')
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.xlabel("Outer Radius (μm)", fontsize=13)
    plt.ylabel("I_ND(Gaussian) - I_ND(Binary)", fontsize=13)
    plt.title("Difference: Gaussian ANNULAR - Binary Ring", fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, cfg.test_distance_end)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, "ind_difference.png"), dpi=300)
    plt.close()


def plot_cell_distribution(state, output_folder):
    """Plot spatial distribution of cells."""
    cfg = state.cfg
    fig, ax = plt.subplots(figsize=(10, 10))

    all_pos = state.all_positions
    si = state.sender_indices
    ri = state.receiver_indices

    other_indices = np.setdiff1d(np.arange(len(all_pos)), np.concatenate([si, ri]))
    ax.scatter(all_pos[other_indices, 0], all_pos[other_indices, 1],
               s=1, c='lightgray', alpha=0.3, label='Other cells')

    # Receiver plotting
    if state.active_receiver_indices is not None:
        ax.scatter(all_pos[state.silent_receiver_indices, 0],
                   all_pos[state.silent_receiver_indices, 1],
                   s=2, c='lightblue', alpha=0.5, label='Silent receivers')
        ax.scatter(all_pos[state.active_receiver_indices, 0],
                   all_pos[state.active_receiver_indices, 1],
                   s=2, c='blue', alpha=0.6, label='Active receivers')
    else:
        ax.scatter(all_pos[ri, 0], all_pos[ri, 1],
                   s=2, c='blue', alpha=0.5, label='Receivers')

    # Sender plotting
    if cfg.sender_positions in ('center_silent_distributed',) and cfg.n_silent > 0:
        ax.scatter(all_pos[state.silent_indices, 0], all_pos[state.silent_indices, 1],
                   s=20, c='red', marker='D', alpha=0.6, edgecolors='black',
                   linewidths=0.3, label='Silent senders', zorder=9)
        ax.scatter(all_pos[state.active_indices, 0], all_pos[state.active_indices, 1],
                   s=100, c='red', marker='*', edgecolors='black', linewidths=0.5,
                   label='Active senders', zorder=10)
    elif state.position_dict is not None:
        # Multi-position mode
        ax.scatter(all_pos[state.silent_indices, 0], all_pos[state.silent_indices, 1],
                   s=20, c='red', marker='D', alpha=0.6, edgecolors='black',
                   linewidths=0.3, label='Silent senders', zorder=9)
        if state.expressing_mask is not None:
            silent_active = state.active_indices[~state.expressing_mask]
            if len(silent_active) > 0:
                ax.scatter(all_pos[silent_active, 0], all_pos[silent_active, 1],
                           s=80, c='white', marker='o', edgecolors='red', linewidths=2,
                           label=f'Active-OFF ({len(silent_active)})', zorder=9)
        n_pos = len(state.position_dict)
        color_cycle = plt.cm.tab10(np.linspace(0, 1, max(n_pos, 10)))
        marker_cycle = ['*', '^', 'v', '>', '<', 's', 'p', 'h', 'D', 'o']
        for i, (pn, coord) in enumerate(state.position_dict.items()):
            n_at = state.active_assignments[pn]
            if n_at > 0:
                size = 100 + n_at * 30
                ax.scatter(coord[0], coord[1], s=size,
                           c=[color_cycle[i % len(color_cycle)]],
                           marker=marker_cycle[i % len(marker_cycle)],
                           edgecolors='black', linewidths=1,
                           label=f'{pn}: {n_at} senders', zorder=10)
    else:
        ax.scatter(all_pos[si, 0], all_pos[si, 1],
                   s=100, c='red', marker='*', edgecolors='black', linewidths=0.5,
                   label='Senders', zorder=10)

    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Cell Distribution ({state.frac_str} Receivers)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"cell_distribution_{state.frac_str}.png"), dpi=300)
    plt.close()


def plot_concentration_gradient(state, output_folder):
    """Plot concentration field as spatial gradient."""
    cfg = state.cfg
    center = np.array(cfg.center)
    fig, ax = plt.subplots(figsize=(10, 10))

    scatter = ax.scatter(state.all_positions[:, 0], state.all_positions[:, 1],
                         c=state.concentrations, s=1, cmap='hot', alpha=0.6,
                         norm=plt.matplotlib.colors.LogNorm(
                             vmin=cfg.vmin_conc, vmax=cfg.vmax_conc))

    # Mark sender positions
    if state.position_dict is not None:
        n_pos = len(state.position_dict)
        bright_colors = ['lime', 'yellow', 'cyan', 'magenta', 'white',
                         'springgreen', 'gold', 'deepskyblue', 'hotpink', 'lavender']
        marker_cycle = ['*', '^', 'v', '>', '<', 's', 'p', 'h', 'D', 'o']
        if cfg.n_silent > 0:
            ax.scatter(state.all_positions[state.silent_indices, 0],
                       state.all_positions[state.silent_indices, 1],
                       s=50, c='cyan', marker='D', alpha=0.8,
                       edgecolors='black', linewidths=0.5, label='Silent senders', zorder=9)
        for i, (pn, coord) in enumerate(state.position_dict.items()):
            n_at = state.active_assignments[pn]
            if n_at > 0:
                size = 150 + n_at * 40
                ax.scatter(coord[0], coord[1], s=size,
                           c=bright_colors[i % len(bright_colors)],
                           marker=marker_cycle[i % len(marker_cycle)],
                           edgecolors='black', linewidths=1.5,
                           label=f'{pn}: {n_at}', zorder=10)
                circle = plt.Circle(coord, state.lambda_val, fill=False,
                                    edgecolor=bright_colors[i % len(bright_colors)],
                                    linewidth=1.5, linestyle='--', alpha=0.6)
                ax.add_patch(circle)
    elif cfg.sender_positions == 'center_silent_distributed' and cfg.n_silent > 0:
        ax.scatter(state.all_positions[state.silent_indices, 0],
                   state.all_positions[state.silent_indices, 1],
                   s=50, c='cyan', marker='D', alpha=0.8,
                   edgecolors='black', linewidths=0.5, label='Silent senders', zorder=9)
        ax.scatter(state.all_positions[state.active_indices, 0],
                   state.all_positions[state.active_indices, 1],
                   s=100, c='cyan', marker='*', edgecolors='black', linewidths=1,
                   label='Active senders', zorder=10)
        circle = plt.Circle(center, state.lambda_val, fill=False,
                            edgecolor='cyan', linewidth=2, linestyle='--',
                            label=f'λ = {state.lambda_val:.0f} μm')
        ax.add_patch(circle)
    else:
        ax.scatter(center[0], center[1], s=200, c='cyan', marker='*',
                   edgecolors='black', linewidths=1, label='Senders', zorder=10)
        circle = plt.Circle(center, state.lambda_val, fill=False,
                            edgecolor='cyan', linewidth=2, linestyle='--',
                            label=f'λ = {state.lambda_val:.0f} μm')
        ax.add_patch(circle)

    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title(f'Concentration Gradient ({state.frac_str} Receivers)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"concentration_gradient_{state.frac_str}.png"), dpi=300)
    plt.close()


def plot_expression_distribution(state, output_folder):
    """Plot factor expression distribution (for stochastic sender models)."""
    if state.expressing_mask is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    ax1.hist(state.factor_expr, bins=50, color='gray', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Factor Expression', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('All Cells - Expression Distribution', fontsize=14, fontweight='bold')
    ax1.set_yscale('log')

    ax2 = axes[1]
    active_expr = state.factor_expr[state.active_indices]
    expressing_expr = active_expr[state.expressing_mask]
    silent_expr = active_expr[~state.expressing_mask]
    if len(silent_expr) > 0:
        ax2.hist(silent_expr, bins=20, color='lightcoral', alpha=0.7,
                 edgecolor='black', label=f'OFF ({len(silent_expr)})')
    if len(expressing_expr) > 0:
        ax2.hist(expressing_expr, bins=20, color='green', alpha=0.7,
                 edgecolor='black', label=f'ON ({len(expressing_expr)})')
        ax2.axvline(np.mean(expressing_expr), color='darkgreen', linestyle='--',
                    linewidth=2, label='Mean (ON)')
    ax2.set_xlabel('Factor Expression', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Active Senders - Stochastic Expression', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"expression_distribution_{state.frac_str}.png"), dpi=300)
    plt.close()


def plot_vst_distributions(state, output_folder):
    """Plot raw vs VST expression distributions (demo_vst style)."""
    if state.factor_vst is None:
        return
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    si = state.sender_indices
    ri = state.receiver_indices
    other = np.setdiff1d(np.arange(len(state.factor_raw)), np.concatenate([si, ri]))

    # Factor Raw
    ax = axes[0, 0]
    ax.hist(state.factor_raw[other], bins=50, alpha=0.5, label='Other', color='gray', density=True)
    ax.hist(state.factor_raw[si], bins=20, alpha=0.7, label='Senders', color='red', density=True)
    ax.set_xlabel('Raw Expression')
    ax.set_ylabel('Density')
    ax.set_title('Factor Gene - Raw Counts')
    ax.legend()
    ax.set_xlim(0, np.percentile(state.factor_raw, 99))

    # Factor VST
    ax = axes[0, 1]
    ax.hist(state.factor_vst[other], bins=50, alpha=0.5, label='Other', color='gray', density=True)
    ax.hist(state.factor_vst[si], bins=20, alpha=0.7, label='Senders', color='red', density=True)
    ax.axvline(0, color='black', linestyle='--', linewidth=1, label='Zero')
    ax.set_xlabel('VST Expression')
    ax.set_ylabel('Density')
    ax.set_title('Factor Gene - VST Normalized')
    ax.legend()

    # Response Raw
    ax = axes[1, 0]
    ax.hist(state.responsive_raw[other], bins=50, alpha=0.5, label='Other', color='gray', density=True)
    ax.hist(state.responsive_raw[ri], bins=50, alpha=0.7, label='Receivers', color='blue', density=True)
    ax.set_xlabel('Raw Expression')
    ax.set_ylabel('Density')
    ax.set_title('Responsive Gene - Raw Counts')
    ax.legend()
    ax.set_xlim(0, np.percentile(state.responsive_raw, 99))

    # Response VST
    ax = axes[1, 1]
    ax.hist(state.responsive_vst[other], bins=50, alpha=0.5, label='Other', color='gray', density=True)
    ax.hist(state.responsive_vst[ri], bins=50, alpha=0.7, label='Receivers', color='blue', density=True)
    ax.axvline(0, color='black', linestyle='--', linewidth=1, label='Zero')
    ax.set_xlabel('VST Expression')
    ax.set_ylabel('Density')
    ax.set_title('Responsive Gene - VST Normalized')
    ax.legend()

    plt.suptitle(f'Expression Distributions ({state.frac_str} Receivers)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"expression_dist_{state.frac_str}.png"), dpi=300)
    plt.close()


def plot_response_distribution(state, output_folder):
    """Plot response gene diagnostic (demo3b style)."""
    if state.responding_mask is None or state.response_probs is None:
        return
    cfg = state.cfg
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    C = state.concentrations[state.receiver_indices]
    receiver_expr = state.responsive_expr[state.receiver_indices]

    # Top-left: P(respond) vs concentration
    ax1 = axes[0, 0]
    n_bins = 50
    bin_edges = np.logspace(np.log10(max(C.min(), 1e-6)), np.log10(max(C.max(), 1e-5)), n_bins + 1)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    observed_rates = []
    expected_probs = []
    for i in range(n_bins):
        mask = (C >= bin_edges[i]) & (C < bin_edges[i + 1])
        if np.sum(mask) > 5:
            observed_rates.append(np.mean(state.responding_mask[mask]))
            expected_probs.append(np.mean(state.response_probs[mask]))
        else:
            observed_rates.append(np.nan)
            expected_probs.append(np.nan)
    ax1.scatter(bin_centers, observed_rates, s=50, c='blue', alpha=0.7, label='Observed', zorder=5)
    ax1.plot(bin_centers, expected_probs, 'r-', linewidth=2, label='Expected (Hill)', zorder=4)
    ax1.axvline(cfg.Kd, color='green', linestyle='--', linewidth=2, label=f'Kd={cfg.Kd}')
    ax1.set_xscale('log')
    ax1.set_xlabel('Concentration', fontsize=12)
    ax1.set_ylabel('P(respond)', fontsize=12)
    ax1.set_title('Response Probability vs Concentration', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.05, 1.05)

    # Top-right: expression histogram
    ax2 = axes[0, 1]
    responding_expr = receiver_expr[state.responding_mask]
    non_responding_expr = receiver_expr[~state.responding_mask]
    if len(non_responding_expr) > 0:
        ax2.hist(non_responding_expr, bins=30, color='lightcoral', alpha=0.7,
                 edgecolor='black', label=f'Non-responding ({len(non_responding_expr)})')
    if len(responding_expr) > 0:
        ax2.hist(responding_expr, bins=30, color='green', alpha=0.7,
                 edgecolor='black', label=f'Responding ({len(responding_expr)})')
    ax2.set_xlabel('Response Gene Expression', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Response Expression Distribution', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.set_yscale('log')

    # Bottom-left: scatter
    ax3 = axes[1, 0]
    ax3.scatter(C[~state.responding_mask], receiver_expr[~state.responding_mask],
                s=5, c='lightcoral', alpha=0.3, label='Non-responding')
    ax3.scatter(C[state.responding_mask], receiver_expr[state.responding_mask],
                s=5, c='green', alpha=0.5, label='Responding')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('Concentration', fontsize=12)
    ax3.set_ylabel('Response Expression', fontsize=12)
    ax3.set_title('Expression vs Concentration', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10, loc='lower right')
    ax3.grid(True, alpha=0.3)

    # Bottom-right: summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    n_resp = np.sum(state.responding_mask)
    n_tot = len(state.responding_mask)
    summary = (f"  Total Receivers: {n_tot:,}\n"
               f"  Responding: {n_resp:,} ({100 * n_resp / n_tot:.1f}%)\n"
               f"  Mean P(respond): {np.mean(state.response_probs):.3f}\n"
               f"  Conc range: {C.min():.2e} - {C.max():.2e}")
    ax4.text(0.1, 0.9, summary, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"response_distribution_{state.frac_str}.png"), dpi=300)
    plt.close()


# =============================================================================
# VST Dispatch
# =============================================================================

def _apply_vst(raw_expr, method):
    """Apply VST transform by method name."""
    if method == 'log1p':
        return apply_vst_log1p(raw_expr)
    elif method == 'pearson':
        return apply_vst_pearson(raw_expr)
    elif method == 'shifted_log':
        return apply_vst_shifted_log(raw_expr)
    else:
        raise ValueError(f"Unknown VST method: {method}")


# =============================================================================
# Main Simulation Loop
# =============================================================================

def run_experiment(cfg, output_folder='./output', seed=42):
    """Run the full experiment defined by cfg."""
    np.random.seed(seed)
    os.makedirs(output_folder, exist_ok=True)

    center = np.array(cfg.center)
    n_total = cfg.n_total
    max_radius = cfg.max_radius
    domain_area = np.pi * max_radius ** 2
    n_senders = cfg.n_active + cfg.n_silent

    print("=" * 70)
    print(f"EXPERIMENT: {cfg.name}")
    print("=" * 70)
    print(f"  D={cfg.D}, k_max={cfg.k_max}, Kd={cfg.Kd}")
    print(f"  Factor model: {cfg.factor_model}, Response model: {cfg.response_model}")
    print(f"  Sender positions: {cfg.sender_positions}")
    print(f"  I_ND methods: {cfg.ind_methods}")
    if cfg.vst_method:
        print(f"  VST: {cfg.vst_method}")
    print("=" * 70)

    # Compute test distances
    d_start = cfg.test_distance_start if cfg.test_distance_start is not None else cfg.BANDWIDTH / 2
    test_distances = np.arange(d_start, cfg.test_distance_end, cfg.test_distance_step)

    # Generate positions
    all_positions_orig = generate_positions(n_total, center, max_radius)

    # Setup multi-position senders if needed
    position_dict = None
    active_assignments = None
    sender_pos_list = None

    if cfg.sender_positions in ('multi_fixed', 'multi_random'):
        if cfg.sender_positions == 'multi_random':
            position_dict = get_random_positions(cfg.n_active_positions, center, max_radius, cfg.min_separation)
        else:
            position_dict = get_5_positions(center, cfg.position_offset)
        print(f"\n  {len(position_dict)} sender positions:")
        for pn, coord in position_dict.items():
            print(f"    {pn}: ({coord[0]:.0f}, {coord[1]:.0f})")

    # Fixed sender setup (if fix_senders_across_fractions)
    fixed_sender_indices = None
    fixed_active_indices = None
    fixed_silent_indices = None

    if cfg.fix_senders_across_fractions:
        all_positions = all_positions_orig.copy()
        all_indices = np.arange(n_total)
        fixed_sender_indices = np.random.choice(all_indices, n_senders, replace=False)
        fixed_active_indices = fixed_sender_indices[:cfg.n_active]
        fixed_silent_indices = fixed_sender_indices[cfg.n_active:]

        if position_dict is not None:
            active_assignments, sender_pos_list = distribute_active_senders(cfg.n_active, position_dict)
            for i, (pn, coord) in enumerate(sender_pos_list):
                all_positions[fixed_active_indices[i]] = coord
            print(f"\n  Active sender distribution (fixed):")
            for pn, count in active_assignments.items():
                if count > 0:
                    print(f"    {pn}: {count} senders")

    all_results = {}

    for frac in cfg.receiver_fractions:
        print(f"\nProcessing {frac * 100:.0f}% Receivers...")

        # Position setup
        if cfg.fix_senders_across_fractions:
            all_pos = all_positions.copy() if cfg.sender_positions not in ('multi_fixed', 'multi_random') else all_positions
            sender_indices = fixed_sender_indices
            active_indices = fixed_active_indices
            silent_indices = fixed_silent_indices
        else:
            all_pos = all_positions_orig.copy()
            all_indices = np.arange(n_total)
            sender_indices = np.random.choice(all_indices, n_senders, replace=False)
            active_indices = sender_indices[:cfg.n_active]
            silent_indices = sender_indices[cfg.n_active:]

        # Position senders
        if cfg.sender_positions == 'center':
            all_pos[active_indices] = center
            all_pos[silent_indices] = center
        elif cfg.sender_positions == 'center_silent_distributed':
            all_pos[active_indices] = center
            # silent remain at random positions
        elif cfg.sender_positions in ('multi_fixed', 'multi_random'):
            if not cfg.fix_senders_across_fractions:
                if position_dict is not None:
                    active_assignments, sender_pos_list = distribute_active_senders(cfg.n_active, position_dict)
                    for i, (pn, coord) in enumerate(sender_pos_list):
                        all_pos[active_indices[i]] = coord

        # Select receivers
        all_indices = np.arange(n_total)
        non_sender_indices = np.setdiff1d(all_indices, sender_indices)
        n_receivers = int(n_total * frac)
        receiver_indices = np.random.choice(non_sender_indices, n_receivers, replace=False)

        # Split receivers into active/silent if needed
        active_receiver_indices = None
        silent_receiver_indices = None
        if cfg.receiver_silent_fraction > 0:
            n_active_recv = int(n_receivers * (1 - cfg.receiver_silent_fraction))
            active_receiver_indices = receiver_indices[:n_active_recv]
            silent_receiver_indices = receiver_indices[n_active_recv:]

        # ===== Factor Expression =====
        expressing_mask = None

        if cfg.vst_method:
            # VST mode: generate raw counts then transform
            factor_raw = cfg.BASAL_FACTOR * np.random.lognormal(0, cfg.sigma_f, n_total)
            # Zero-inflate
            if cfg.zero_inflate_factor > 0:
                zero_mask = np.random.rand(n_total) < cfg.zero_inflate_factor
                zero_mask[sender_indices] = False
                factor_raw[zero_mask] = 0.0
            factor_raw[active_indices] = cfg.HIGH_FACTOR * np.random.lognormal(0, cfg.sigma_f, cfg.n_active)
            if cfg.n_silent > 0:
                factor_raw[silent_indices] = cfg.BASAL_FACTOR * np.random.lognormal(0, cfg.sigma_f, cfg.n_silent)
            factor_vst = _apply_vst(factor_raw, cfg.vst_method)
            factor_expr_for_diffusion = factor_raw
            factor_expr_for_ind = factor_vst
        elif cfg.factor_model == 'deterministic':
            factor_expr = cfg.BASAL_FACTOR * np.random.lognormal(0, cfg.sigma_f, n_total)
            factor_expr[active_indices] = cfg.HIGH_FACTOR * np.random.lognormal(0, cfg.sigma_f, cfg.n_active)
            if cfg.n_silent > 0:
                if cfg.silent_expr_zero:
                    factor_expr[silent_indices] = 0.0
                else:
                    factor_expr[silent_indices] = cfg.BASAL_FACTOR * np.random.lognormal(0, cfg.sigma_f, cfg.n_silent)
            factor_expr_for_diffusion = factor_expr
            factor_expr_for_ind = factor_expr
            factor_raw = None
            factor_vst = None
        elif cfg.factor_model == 'bernoulli_mixture':
            factor_expr, active_sender_idx = generate_bernoulli_factor_expression(
                n_total, sender_indices, cfg.p_s,
                cfg.HIGH_FACTOR, cfg.BASAL_FACTOR, cfg.sigma_f, cfg.sigma_f_b
            )
            factor_expr_for_diffusion = factor_expr
            factor_expr_for_ind = factor_expr
            factor_raw = None
            factor_vst = None
        elif cfg.factor_model == 'stochastic':
            factor_expr, expressing_mask = generate_stochastic_expression(
                n_senders=n_total, n_active=cfg.n_active, active_indices=active_indices,
                basal_level=cfg.BASAL_FACTOR, high_level=cfg.HIGH_FACTOR,
                sigma_noise=cfg.sigma_f, p_express=cfg.p_express,
                expr_cv=cfg.expr_cv, use_gamma=cfg.use_gamma_dist
            )
            factor_expr_for_diffusion = factor_expr
            factor_expr_for_ind = factor_expr
            factor_raw = None
            factor_vst = None
        elif cfg.factor_model == 'stochastic_ref':
            factor_expr, expressing_mask = generate_stochastic_expression_ref(
                n_senders=n_total, n_active=cfg.n_active, active_indices=active_indices,
                F_basal=cfg.BASAL_FACTOR, F_high=cfg.HIGH_FACTOR,
                sigma_f=cfg.sigma_f, sigma_f_basal=cfg.sigma_f_basal,
                p_express=cfg.p_express
            )
            factor_expr_for_diffusion = factor_expr
            factor_expr_for_ind = factor_expr
            factor_raw = None
            factor_vst = None
        else:
            raise ValueError(f"Unknown factor_model: {cfg.factor_model}")

        if not cfg.vst_method:
            factor_expr = factor_expr_for_diffusion  # alias

        # ===== Diffusion =====
        n_density = n_receivers / domain_area

        if cfg.sender_positions in ('multi_fixed', 'multi_random') and position_dict is not None:
            p_r_eff = cfg.p_r if cfg.response_model == 'bernoulli_constant' else 1.0
            if cfg.response_model == 'stochastic_hill':
                p_r_eff = cfg.p_respond_max
            concentrations, lambda_val = solve_concentration_field_MM_multipos(
                all_pos[sender_indices], factor_expr_for_diffusion[sender_indices],
                all_pos, n_density, cfg.D, cfg.k_max, cfg.Kd,
                position_dict=position_dict, active_assignments=active_assignments,
                secretion_rate=cfg.secretion_rate, p_r=p_r_eff,
                active_threshold=cfg.active_threshold
            )
        else:
            if cfg.response_model == 'bernoulli_constant':
                n_density_eff = n_density * cfg.p_r
            elif cfg.response_model == 'bernoulli_hill':
                n_density_eff = n_density * cfg.p_r_max
            else:
                n_density_eff = n_density
            at = cfg.vst_active_threshold if cfg.vst_method else cfg.active_threshold
            concentrations, lambda_val = solve_concentration_field_MM(
                all_pos[sender_indices], factor_expr_for_diffusion[sender_indices],
                all_pos, n_density_eff, cfg.D, cfg.k_max, cfg.Kd,
                active_threshold=at
            )

        # ===== Response Expression =====
        responding_mask = None
        response_probs = None
        mean_p_r = None

        if cfg.vst_method:
            responsive_raw = cfg.BASAL_RESPONSIVE * np.random.lognormal(0, cfg.sigma_r, n_total)
            if cfg.zero_inflate_responsive > 0:
                zero_mask_r = np.random.rand(n_total) < cfg.zero_inflate_responsive
                zero_mask_r[receiver_indices] = False
                responsive_raw[zero_mask_r] = 0.0
            activation = concentrations[receiver_indices] / (cfg.Kd + concentrations[receiver_indices])
            responsive_raw[receiver_indices] = cfg.BASAL_RESPONSIVE * (1 + cfg.FOLD_CHANGE * activation)
            responsive_vst = _apply_vst(responsive_raw, cfg.vst_method)
            responsive_expr_for_ind = responsive_vst
            responsive_expr = responsive_raw
        elif cfg.response_model == 'deterministic':
            responsive_expr = cfg.BASAL_RESPONSIVE * np.random.lognormal(0, cfg.sigma_r, n_total)
            ri_for_activation = active_receiver_indices if active_receiver_indices is not None else receiver_indices
            activation = concentrations[ri_for_activation] / (cfg.Kd + concentrations[ri_for_activation])
            responsive_expr[ri_for_activation] = cfg.BASAL_RESPONSIVE * (1 + cfg.FOLD_CHANGE * activation)
            responsive_expr_for_ind = responsive_expr
            responsive_raw = None
            responsive_vst = None
        elif cfg.response_model == 'bernoulli_constant':
            responsive_expr, responding_recv_idx = generate_bernoulli_response_constant(
                n_total, receiver_indices, cfg.p_r, concentrations,
                cfg.BASAL_RESPONSIVE, cfg.FOLD_CHANGE, cfg.Kd, cfg.sigma_r, cfg.sigma_r_b
            )
            responsive_expr_for_ind = responsive_expr
            responsive_raw = None
            responsive_vst = None
        elif cfg.response_model == 'bernoulli_hill':
            responsive_expr, responding_recv_idx, p_r_values, mean_p_r = generate_bernoulli_response_hill(
                n_total, receiver_indices, concentrations,
                cfg.BASAL_RESPONSIVE, cfg.FOLD_CHANGE, cfg.Kd, cfg.sigma_r, cfg.sigma_r_b,
                p_r_max=cfg.p_r_max, K_p=cfg.K_p, hill_n=cfg.hill_n
            )
            responsive_expr_for_ind = responsive_expr
            responsive_raw = None
            responsive_vst = None
        elif cfg.response_model == 'stochastic_hill':
            responsive_expr, responding_mask, response_probs = generate_stochastic_response(
                n_total=n_total, receiver_indices=receiver_indices,
                concentrations=concentrations,
                basal_responsive=cfg.BASAL_RESPONSIVE, fold_change=cfg.FOLD_CHANGE,
                sigma_r=cfg.sigma_r, sigma_r_basal=cfg.sigma_r_b,
                Kd=cfg.Kd, p_respond_max=cfg.p_respond_max,
                hill_coef=cfg.response_hill_coef
            )
            responsive_expr_for_ind = responsive_expr
            responsive_raw = None
            responsive_vst = None
        else:
            raise ValueError(f"Unknown response_model: {cfg.response_model}")

        if not cfg.vst_method and cfg.response_model not in ('bernoulli_constant', 'bernoulli_hill', 'stochastic_hill'):
            responsive_raw = None
            responsive_vst = None

        # ===== Compute I_ND =====
        ri_for_ind = receiver_indices
        if cfg.ind_uses_active_receivers_only and active_receiver_indices is not None:
            ri_for_ind = active_receiver_indices

        ind_curves = {}
        for method in cfg.ind_methods:
            curve = []
            for d in test_distances:
                if method == 'ring':
                    val, n_conn = compute_IND_ring(
                        sender_indices, ri_for_ind, all_pos,
                        factor_expr_for_ind, responsive_expr_for_ind,
                        d, bandwidth=cfg.BANDWIDTH
                    )
                elif method == 'gaussian_annular':
                    val, n_conn = compute_IND_gaussian_annular(
                        sender_indices, ri_for_ind, all_pos,
                        factor_expr_for_ind, responsive_expr_for_ind,
                        outer_radius=d, bandwidth=cfg.BANDWIDTH,
                        sigma_fraction=cfg.sigma_fraction
                    )
                else:
                    raise ValueError(f"Unknown I_ND method: {method}")
                curve.append({'distance': d, 'I_ND': val, 'n_conn': n_conn})
            ind_curves[method] = curve

        n_expressing = int(np.sum(expressing_mask)) if expressing_mask is not None else None
        n_responding = int(np.sum(responding_mask)) if responding_mask is not None else None

        state = SimState(
            cfg=cfg, frac=frac, frac_str=f"{int(frac * 100)}pct",
            all_positions=all_pos,
            sender_indices=sender_indices, active_indices=active_indices,
            silent_indices=silent_indices, receiver_indices=receiver_indices,
            factor_expr=factor_expr_for_ind if not cfg.vst_method else factor_expr_for_diffusion,
            responsive_expr=responsive_expr_for_ind if not cfg.vst_method else responsive_expr,
            concentrations=concentrations, lambda_val=lambda_val,
            ind_curves=ind_curves,
            expressing_mask=expressing_mask,
            responding_mask=responding_mask,
            response_probs=response_probs,
            mean_p_r=mean_p_r,
            active_receiver_indices=active_receiver_indices,
            silent_receiver_indices=silent_receiver_indices,
            factor_raw=factor_raw if cfg.vst_method else None,
            factor_vst=factor_vst if cfg.vst_method else None,
            responsive_raw=responsive_raw if cfg.vst_method else None,
            responsive_vst=responsive_vst if cfg.vst_method else None,
            position_dict=position_dict,
            active_assignments=active_assignments,
            n_expressing=n_expressing,
            n_responding=n_responding,
        )

        all_results[frac] = state

        print(f"  λ = {lambda_val:.0f} μm")
        if n_expressing is not None:
            print(f"  Expressing: {n_expressing}/{cfg.n_active}")
        if n_responding is not None:
            print(f"  Responding: {n_responding}/{n_receivers}")

        # Per-fraction plots
        plot_cell_distribution(state, output_folder)
        plot_concentration_gradient(state, output_folder)
        if expressing_mask is not None:
            plot_expression_distribution(state, output_folder)
        if cfg.vst_method:
            plot_vst_distributions(state, output_folder)
        if responding_mask is not None and response_probs is not None:
            plot_response_distribution(state, output_folder)

    # Summary plots
    plot_ind_summary(cfg, all_results, output_folder)

    print(f"\nDone. Results in {output_folder}")
    return all_results


def generate_positions(n_total, center, max_radius):
    """Generate uniformly distributed random positions in a circular domain."""
    angles = np.random.rand(n_total) * 2 * np.pi
    radii = np.sqrt(np.random.rand(n_total)) * max_radius
    return np.column_stack([
        center[0] + radii * np.cos(angles),
        center[1] + radii * np.sin(angles)
    ])


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Unified I_ND simulation experiments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Available presets:\n  " + "\n  ".join(sorted(PRESETS.keys()))
    )
    parser.add_argument('preset', nargs='?', help='Preset name to run')
    parser.add_argument('--list', action='store_true', help='List available presets')
    parser.add_argument('--seed', type=int, default=42, help='Random seed (default: 42)')
    parser.add_argument('--output', type=str, default='./output', help='Output directory')

    args = parser.parse_args()

    if args.list or args.preset is None:
        print("Available presets:")
        for name in sorted(PRESETS.keys()):
            print(f"  {name}")
        return

    cfg = make_config(args.preset)
    run_experiment(cfg, output_folder=args.output, seed=args.seed)


if __name__ == "__main__":
    main()

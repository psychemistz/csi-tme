#!/usr/bin/env python3
"""
Core simulation functions for I_ND validation experiments.

Shared computational functions used across all experiment variants (v5 series).

Handles:
- Lambda calculation for reaction-diffusion systems
- Concentration field solving (single-source and multi-source)
- Normalized Directional Moran's I (I_ND) computation
- Position generation and sender distribution
- Stochastic expression models
"""

import numpy as np
from scipy import spatial


# =============================================================================
# Lambda & Position Generation
# =============================================================================

def calculate_lambda(D, n_receivers, k_max, Kd, p_r=1.0):
    """Calculate characteristic decay length lambda.

    Parameters
    ----------
    D : float
        Diffusion coefficient (um^2/s)
    n_receivers : float
        Receiver density (cells/um^2)
    k_max : float
        Maximum uptake rate
    Kd : float
        Dissociation constant
    p_r : float, optional
        Probability of responding (for effective density). Default 1.0.

    Returns
    -------
    float
        Lambda (um)
    """
    n_eff = n_receivers * p_r
    if n_eff <= 0:
        return np.inf
    return np.sqrt(D * Kd / (n_eff * k_max))


def generate_circular_positions(n_total, center, max_radius):
    """Generate uniformly distributed random positions in a circular domain.

    Parameters
    ----------
    n_total : int
        Number of positions to generate
    center : array-like
        Center of the circular domain [x, y]
    max_radius : float
        Radius of the domain

    Returns
    -------
    ndarray
        Shape (n_total, 2) array of positions
    """
    angles = np.random.rand(n_total) * 2 * np.pi
    radii = np.sqrt(np.random.rand(n_total)) * max_radius
    return np.column_stack([
        center[0] + radii * np.cos(angles),
        center[1] + radii * np.sin(angles)
    ])


# =============================================================================
# Concentration Field Solvers
# =============================================================================

def solve_concentration_field_MM(sender_positions, sender_factor_expr,
                                  cell_positions, n_receivers,
                                  D, k_max=300.0, Kd=5.0,
                                  active_threshold=1.0):
    """Calculate factor concentration field from a single source (co-located senders).

    Uses 2D steady-state analytical solution: C(r) = Q * exp(-r/lambda) / sqrt(r)

    Parameters
    ----------
    sender_positions : ndarray
        Positions of sender cells
    sender_factor_expr : ndarray
        Expression levels of sender cells
    cell_positions : ndarray
        Positions of all cells to compute concentration at
    n_receivers : float
        Receiver density
    D : float
        Diffusion coefficient
    k_max : float
        Maximum uptake rate
    Kd : float
        Dissociation constant
    active_threshold : float
        Expression threshold to consider a sender "active". Default 1.0.

    Returns
    -------
    tuple
        (concentrations, lambda_val)
    """
    lambda_val = calculate_lambda(D, n_receivers, k_max, Kd)
    concentrations = np.zeros(len(cell_positions))

    # Identify active senders
    active_mask = sender_factor_expr > active_threshold
    active_pos = sender_positions[active_mask]
    active_expr = sender_factor_expr[active_mask]

    if len(active_pos) > 0:
        # All active senders are at center
        center_pos = active_pos[0]
        total_factor = np.sum(active_expr)

        for i, pos in enumerate(cell_positions):
            r = np.linalg.norm(pos - center_pos)
            if r < 1e-3:
                concentrations[i] = total_factor * 100
            else:
                concentrations[i] = total_factor * np.exp(-r / lambda_val) / np.sqrt(r)

    return concentrations, lambda_val


def solve_concentration_field_MM_multipos(sender_positions, sender_factor_expr,
                                           cell_positions, n_receivers,
                                           D, k_max=300.0, Kd=5.0,
                                           position_dict=None, active_assignments=None,
                                           secretion_rate=1.0, p_r=1.0,
                                           active_threshold=1.0):
    """Calculate concentration field from MULTIPLE source positions using superposition.

    Parameters
    ----------
    sender_positions : ndarray
        Positions of all senders (including silent)
    sender_factor_expr : ndarray
        Expression levels of all senders
    cell_positions : ndarray
        Positions of all cells
    n_receivers : float
        Receiver density
    D : float
        Diffusion coefficient
    k_max : float
        Maximum uptake rate
    Kd : float
        Dissociation constant
    position_dict : dict
        Dictionary of positions {name: coordinate}
    active_assignments : dict
        Dictionary of {position_name: count}
    secretion_rate : float
        Scaling factor for expression -> concentration
    p_r : float
        Probability of responding (for effective density)
    active_threshold : float
        Expression threshold to consider a sender "active"

    Returns
    -------
    tuple
        (concentrations, lambda_val)
    """
    lambda_val = calculate_lambda(D, n_receivers, k_max, Kd, p_r=p_r)
    concentrations = np.zeros(len(cell_positions))

    active_mask = sender_factor_expr > active_threshold
    active_expr = sender_factor_expr[active_mask]
    active_pos = sender_positions[active_mask]

    if len(active_pos) > 0 and position_dict is not None:
        unique_positions = {}

        for i, (pos, expr) in enumerate(zip(active_pos, active_expr)):
            pos_tuple = tuple(pos)
            if pos_tuple not in unique_positions:
                unique_positions[pos_tuple] = {'pos': pos, 'total_expr': 0}
            unique_positions[pos_tuple]['total_expr'] += expr

        for pos_tuple, data in unique_positions.items():
            source_pos = data['pos']
            total_factor = data['total_expr'] * secretion_rate

            for i, cell_pos in enumerate(cell_positions):
                r = np.linalg.norm(cell_pos - source_pos)
                if r < 1e-3:
                    concentrations[i] += total_factor * 100
                else:
                    concentrations[i] += total_factor * np.exp(-r / lambda_val) / np.sqrt(r)

    return concentrations, lambda_val


# =============================================================================
# I_ND Computation
# =============================================================================

def compute_IND_ring(sender_indices, receiver_indices,
                     all_positions, factor_expr, responsive_expr,
                     distance, bandwidth=100):
    """Compute I_ND using a ring (annular) neighborhood centered at `distance`.

    The ring spans [distance - bandwidth/2, distance + bandwidth/2].

    Parameters
    ----------
    sender_indices : ndarray
        Indices of sender cells
    receiver_indices : ndarray
        Indices of receiver cells
    all_positions : ndarray
        Positions of all cells
    factor_expr : ndarray
        Factor (sender gene) expression for all cells
    responsive_expr : ndarray
        Response (receiver gene) expression for all cells
    distance : float
        Center distance of the ring
    bandwidth : float
        Width of the ring

    Returns
    -------
    tuple
        (I_ND, total_connections)
    """
    # 1. Global Normalization
    mu_f, sigma_f = np.mean(factor_expr), np.std(factor_expr)
    mu_r, sigma_r = np.mean(responsive_expr), np.std(responsive_expr)

    z_s = (factor_expr[sender_indices] - mu_f) / (sigma_f + 1e-10)
    z_r = (responsive_expr[receiver_indices] - mu_r) / (sigma_r + 1e-10)

    sender_pos = all_positions[sender_indices]
    receiver_pos = all_positions[receiver_indices]

    # 2. Build Weight Matrix (Ring)
    dists = spatial.distance.cdist(sender_pos, receiver_pos)

    half_bw = bandwidth / 2.0
    lower = distance - half_bw
    upper = distance + half_bw

    W = ((dists > lower) & (dists <= upper)).astype(float)

    # Row-normalize
    row_sums = W.sum(axis=1, keepdims=True)

    # Check connections
    total_connections = np.sum(W)
    if total_connections == 0:
        return 0.0, 0

    row_sums[row_sums == 0] = 1.0
    W_tilde = W / row_sums

    # 3. Compute Metrics
    mean_neighbor_z = np.dot(W_tilde, z_r)
    numerator = np.dot(z_s, mean_neighbor_z)
    norm_z_s = np.linalg.norm(z_s)
    norm_W_z_r = np.linalg.norm(mean_neighbor_z)

    if norm_z_s > 1e-10 and norm_W_z_r > 1e-10:
        I_ND = numerator / (norm_z_s * norm_W_z_r)
    else:
        I_ND = 0.0

    return I_ND, int(total_connections)


def compute_IND_gaussian_annular(sender_indices, receiver_indices,
                                  all_positions, factor_expr, responsive_expr,
                                  outer_radius, bandwidth=100, sigma_fraction=3.0):
    """Compute I_ND using Gaussian annular kernel.

    Gaussian is centered at sender (d=0), sigma = outer_radius / sigma_fraction.
    Weights are applied within [outer_radius - bandwidth, outer_radius].

    Parameters
    ----------
    sender_indices, receiver_indices, all_positions, factor_expr, responsive_expr
        Same as compute_IND_ring
    outer_radius : float
        Outer edge of the annular region
    bandwidth : float
        Width of the annular region
    sigma_fraction : float
        sigma = outer_radius / sigma_fraction

    Returns
    -------
    tuple
        (I_ND, total_connections)
    """
    mu_f, sigma_f = np.mean(factor_expr), np.std(factor_expr)
    mu_r, sigma_r = np.mean(responsive_expr), np.std(responsive_expr)

    z_s = (factor_expr[sender_indices] - mu_f) / (sigma_f + 1e-10)
    z_r = (responsive_expr[receiver_indices] - mu_r) / (sigma_r + 1e-10)

    sender_pos = all_positions[sender_indices]
    receiver_pos = all_positions[receiver_indices]

    dists = spatial.distance.cdist(sender_pos, receiver_pos)

    # Annular region
    inner_radius = outer_radius - bandwidth

    # Gaussian centered at d=0 (sender position)
    sigma = outer_radius / sigma_fraction

    W = np.exp(-dists**2 / (2 * sigma**2))

    # Apply annular mask
    W[(dists <= inner_radius) | (dists > outer_radius)] = 0

    # Row-normalize
    row_sums = W.sum(axis=1, keepdims=True)
    total_connections = np.sum(W > 1e-6)

    if total_connections == 0:
        return 0.0, 0

    row_sums[row_sums == 0] = 1.0
    W_tilde = W / row_sums

    # Compute I_ND
    mean_neighbor_z = np.dot(W_tilde, z_r)
    numerator = np.dot(z_s, mean_neighbor_z)
    norm_z_s = np.linalg.norm(z_s)
    norm_W_z_r = np.linalg.norm(mean_neighbor_z)

    if norm_z_s > 1e-10 and norm_W_z_r > 1e-10:
        I_ND = numerator / (norm_z_s * norm_W_z_r)
    else:
        I_ND = 0.0

    return I_ND, int(total_connections)


# =============================================================================
# Multi-Position Sender Helpers
# =============================================================================

def get_5_positions(center, offset_distance):
    """Generate 5 fixed positions: Center, West, East, North, South.

    Parameters
    ----------
    center : array-like
        Center position [x, y]
    offset_distance : float
        Distance from center for W, E, N, S positions

    Returns
    -------
    dict
        {name: coordinate} dictionary
    """
    return {
        'C': np.array(center, dtype=float),
        'W': np.array([center[0] - offset_distance, center[1]]),
        'E': np.array([center[0] + offset_distance, center[1]]),
        'N': np.array([center[0], center[1] + offset_distance]),
        'S': np.array([center[0], center[1] - offset_distance]),
    }


def get_random_positions(n_positions, center, max_radius, min_separation=500.0):
    """Generate n random positions within the domain with minimum separation.

    Positions are placed within 80% of max_radius to avoid edge effects.

    Parameters
    ----------
    n_positions : int
        Number of positions to generate
    center : array-like
        Center of the domain [x, y]
    max_radius : float
        Maximum radius of the circular domain
    min_separation : float
        Minimum distance between any two positions

    Returns
    -------
    dict
        {label: coordinate} dictionary
    """
    positions = {}
    coords_list = []

    effective_radius = max_radius * 0.8

    attempts = 0
    max_attempts = 1000

    while len(positions) < n_positions and attempts < max_attempts:
        angle = np.random.rand() * 2 * np.pi
        r = np.sqrt(np.random.rand()) * effective_radius

        new_pos = np.array([
            center[0] + r * np.cos(angle),
            center[1] + r * np.sin(angle)
        ])

        if len(coords_list) == 0:
            valid = True
        else:
            distances = [np.linalg.norm(new_pos - existing) for existing in coords_list]
            valid = min(distances) >= min_separation

        if valid:
            pos_label = f'P{len(positions) + 1}'
            positions[pos_label] = new_pos
            coords_list.append(new_pos)

        attempts += 1

    if len(positions) < n_positions:
        print(f"  Warning: Only placed {len(positions)}/{n_positions} positions with min_separation={min_separation}")

    return positions


def distribute_active_senders(n_active, position_dict):
    """Distribute n_active senders across positions.

    Each position receives at least one sender.
    Remaining senders are assigned randomly.

    Parameters
    ----------
    n_active : int
        Number of active senders to distribute
    position_dict : dict
        Dictionary of positions {name: coordinate}

    Returns
    -------
    tuple
        (assignments dict, sender_positions list of (name, coord) tuples)
    """
    position_names = list(position_dict.keys())
    P = len(position_names)

    if n_active < P:
        raise ValueError("n_active must be >= number of positions to ensure at least one per position.")

    # Initialize: each position gets one sender
    assignments = {name: 1 for name in position_names}
    sender_positions = [(name, position_dict[name]) for name in position_names]

    # Remaining senders to distribute
    remaining = n_active - P

    for _ in range(remaining):
        chosen = np.random.choice(position_names)
        assignments[chosen] += 1
        sender_positions.append((chosen, position_dict[chosen]))

    return assignments, sender_positions


# =============================================================================
# Stochastic Expression Models
# =============================================================================

def generate_stochastic_expression(n_senders, n_active, active_indices,
                                   basal_level, high_level, sigma_noise,
                                   p_express=0.7, expr_cv=0.5,
                                   use_gamma=True):
    """Generate stochastic factor expression with on/off switching.

    Each active sender has probability p_express of being "on".
    Expression levels drawn from gamma or lognormal distribution.

    Parameters
    ----------
    n_senders : int
        Total number of cells
    n_active : int
        Number of active senders
    active_indices : ndarray
        Indices of active senders
    basal_level : float
        Basal expression level
    high_level : float
        Mean expression level when "on"
    sigma_noise : float
        Noise level for basal expression (lognormal sigma)
    p_express : float
        Probability of being "on"
    expr_cv : float
        Coefficient of variation when "on"
    use_gamma : bool
        If True, use gamma distribution; else lognormal

    Returns
    -------
    tuple
        (factor_expr, expressing_mask)
    """
    factor_expr = basal_level * np.random.lognormal(0, sigma_noise, n_senders)

    expressing_mask = np.random.rand(n_active) < p_express
    n_expressing = np.sum(expressing_mask)

    if n_expressing > 0:
        expressing_indices = active_indices[expressing_mask]

        if use_gamma:
            shape = 1.0 / (expr_cv ** 2)
            scale = high_level * (expr_cv ** 2)
            expr_values = np.random.gamma(shape, scale, n_expressing)
        else:
            sigma_ln = np.sqrt(np.log(1 + expr_cv**2))
            mu_ln = np.log(high_level) - sigma_ln**2 / 2
            expr_values = np.random.lognormal(mu_ln, sigma_ln, n_expressing)

        # Ensure minimum expression above basal
        expr_values = np.maximum(expr_values, basal_level * 2)
        factor_expr[expressing_indices] = expr_values

    return factor_expr, expressing_mask


def generate_stochastic_expression_ref(n_senders, n_active, active_indices,
                                        F_basal, F_high, sigma_f, sigma_f_basal,
                                        p_express=0.7):
    """Generate stochastic factor expression using reference model.

    F_i = S_i * F_high * LogNormal(0, sigma_f^2) + (1-S_i) * F_basal * LogNormal(0, sigma_f_basal^2)
    S_i ~ Bernoulli(p_express)

    Parameters
    ----------
    n_senders : int
        Total number of cells
    n_active : int
        Number of active senders
    active_indices : ndarray
        Indices of active senders
    F_basal : float
        Basal expression level
    F_high : float
        High expression level when "on"
    sigma_f : float
        Lognormal sigma for expressing cells
    sigma_f_basal : float
        Lognormal sigma for non-expressing/basal cells
    p_express : float
        Probability of being "on"

    Returns
    -------
    tuple
        (factor_expr, expressing_mask)
    """
    factor_expr = F_basal * np.random.lognormal(0, sigma_f_basal, n_senders)

    expressing_mask = np.random.rand(n_active) < p_express
    n_expressing = np.sum(expressing_mask)

    if n_expressing > 0:
        expressing_indices = active_indices[expressing_mask]
        factor_expr[expressing_indices] = (
            F_high * np.random.lognormal(0, sigma_f, n_expressing)
        )

    return factor_expr, expressing_mask


def hill_function(C, p_r_max, K_p, n):
    """Hill function for concentration-dependent probability.

    p_r(C) = p_r_max * C^n / (K_p^n + C^n)

    Parameters
    ----------
    C : array-like
        Concentration values
    p_r_max : float
        Maximum probability (saturation level)
    K_p : float
        Concentration at half-maximal probability
    n : float
        Hill coefficient (n > 1: switch-like, n = 1: hyperbolic, n < 1: gradual)

    Returns
    -------
    ndarray
        Probability values in [0, p_r_max]
    """
    C = np.asarray(C)
    C_n = np.power(np.maximum(C, 0), n)
    K_n = np.power(K_p, n)
    return p_r_max * C_n / (K_n + C_n)


def generate_bernoulli_factor_expression(n_total, sender_indices, p_s,
                                          F_high, F_basal, sigma_f, sigma_f_b):
    """Generate factor expression using Bernoulli mixture model.

    F_i = S_i * F_high * LN(0, sigma_f^2) + (1-S_i) * F_basal * LN(0, sigma_f_b^2)

    where S_i ~ Bernoulli(p_s) for sender cells, S_i = 0 for non-senders.

    Parameters
    ----------
    n_total : int
        Total number of cells
    sender_indices : ndarray
        Indices of sender cells
    p_s : float
        Probability of sender being active
    F_high : float
        High expression level when "on"
    F_basal : float
        Basal expression level
    sigma_f : float
        Lognormal sigma for active expression
    sigma_f_b : float
        Lognormal sigma for basal expression

    Returns
    -------
    tuple
        (factor_expr, active_sender_indices)
    """
    factor_expr = np.zeros(n_total)

    # Non-sender cells: always basal expression
    non_sender_mask = np.ones(n_total, dtype=bool)
    non_sender_mask[sender_indices] = False
    n_non_senders = np.sum(non_sender_mask)
    factor_expr[non_sender_mask] = F_basal * np.random.lognormal(0, sigma_f_b, n_non_senders)

    # Sender cells: stochastic S_i ~ Bernoulli(p_s)
    n_senders = len(sender_indices)
    S = np.random.binomial(1, p_s, n_senders)

    # Mixture: S_i * high + (1-S_i) * basal
    high_expr = F_high * np.random.lognormal(0, sigma_f, n_senders)
    basal_expr = F_basal * np.random.lognormal(0, sigma_f_b, n_senders)
    factor_expr[sender_indices] = S * high_expr + (1 - S) * basal_expr

    active_sender_indices = sender_indices[S == 1]
    return factor_expr, active_sender_indices


def generate_bernoulli_response_constant(n_total, receiver_indices, p_r,
                                          concentrations, R_basal, FC, Kd,
                                          sigma_r, sigma_r_b):
    """Generate response expression with constant Bernoulli probability.

    R_j = B_j * R_basal * (1 + FC * Act_j) * LN(0, sigma_r^2)
        + (1-B_j) * R_basal * LN(0, sigma_r_b^2)

    where B_j ~ Bernoulli(p_r) for receiver cells.

    Parameters
    ----------
    n_total : int
        Total number of cells
    receiver_indices : ndarray
        Indices of receiver cells
    p_r : float
        Probability of receiver responding
    concentrations : ndarray
        Concentration at each cell position
    R_basal : float
        Basal response expression
    FC : float
        Fold change
    Kd : float
        Dissociation constant
    sigma_r : float
        Lognormal sigma for responding cells
    sigma_r_b : float
        Lognormal sigma for non-responding cells

    Returns
    -------
    tuple
        (responsive_expr, responding_receiver_indices)
    """
    responsive_expr = np.zeros(n_total)

    # Non-receiver cells: always basal expression
    non_receiver_mask = np.ones(n_total, dtype=bool)
    non_receiver_mask[receiver_indices] = False
    n_non_receivers = np.sum(non_receiver_mask)
    responsive_expr[non_receiver_mask] = R_basal * np.random.lognormal(0, sigma_r_b, n_non_receivers)

    # Receiver cells: stochastic B_j ~ Bernoulli(p_r)
    n_receivers = len(receiver_indices)
    B = np.random.binomial(1, p_r, n_receivers)

    C_receivers = concentrations[receiver_indices]
    activation = C_receivers / (Kd + C_receivers)

    activated_expr = R_basal * (1 + FC * activation) * np.random.lognormal(0, sigma_r, n_receivers)
    basal_expr = R_basal * np.random.lognormal(0, sigma_r_b, n_receivers)
    responsive_expr[receiver_indices] = B * activated_expr + (1 - B) * basal_expr

    responding_receiver_indices = receiver_indices[B == 1]
    return responsive_expr, responding_receiver_indices


def generate_bernoulli_response_hill(n_total, receiver_indices, concentrations,
                                      R_basal, FC, Kd, sigma_r, sigma_r_b,
                                      p_r_max=1.0, K_p=10.0, hill_n=1.0):
    """Generate response expression with concentration-dependent Bernoulli probability.

    p_r(C_j) = p_r_max * C_j^n / (K_p^n + C_j^n)  (Hill function)
    B_j ~ Bernoulli(p_r(C_j))

    R_j = B_j * R_basal * (1 + FC * Act_j) * LN(0, sigma_r^2)
        + (1-B_j) * R_basal * LN(0, sigma_r_b^2)

    Parameters
    ----------
    n_total : int
        Total number of cells
    receiver_indices : ndarray
        Indices of receiver cells
    concentrations : ndarray
        Concentration at each cell position
    R_basal : float
        Basal response expression
    FC : float
        Fold change
    Kd : float
        Dissociation constant
    sigma_r : float
        Lognormal sigma for responding cells
    sigma_r_b : float
        Lognormal sigma for non-responding cells
    p_r_max : float
        Maximum receiving probability
    K_p : float
        Concentration at half-maximal probability
    hill_n : float
        Hill coefficient

    Returns
    -------
    tuple
        (responsive_expr, responding_receiver_indices, p_r_values, mean_p_r)
    """
    responsive_expr = np.zeros(n_total)

    # Non-receiver cells: always basal expression
    non_receiver_mask = np.ones(n_total, dtype=bool)
    non_receiver_mask[receiver_indices] = False
    n_non_receivers = np.sum(non_receiver_mask)
    responsive_expr[non_receiver_mask] = R_basal * np.random.lognormal(0, sigma_r_b, n_non_receivers)

    # Receiver cells: concentration-dependent p_r via Hill function
    n_receivers = len(receiver_indices)
    C_receivers = concentrations[receiver_indices]

    p_r_values = hill_function(C_receivers, p_r_max, K_p, hill_n)
    B = np.random.binomial(1, p_r_values)

    activation = C_receivers / (Kd + C_receivers)

    activated_expr = R_basal * (1 + FC * activation) * np.random.lognormal(0, sigma_r, n_receivers)
    basal_expr = R_basal * np.random.lognormal(0, sigma_r_b, n_receivers)
    responsive_expr[receiver_indices] = B * activated_expr + (1 - B) * basal_expr

    responding_receiver_indices = receiver_indices[B == 1]
    mean_p_r = np.mean(p_r_values)
    return responsive_expr, responding_receiver_indices, p_r_values, mean_p_r


# =============================================================================
# VST Normalization Functions
# =============================================================================

def apply_vst_log1p(raw_expr):
    """Log1p VST normalization with centering.

    Applies log1p transform then centers by global mean.
    Non-expressing cells get small negative values.

    Parameters
    ----------
    raw_expr : ndarray
        Raw expression values

    Returns
    -------
    ndarray
        VST-normalized expression
    """
    log_expr = np.log1p(raw_expr)
    return log_expr - np.mean(log_expr)


def apply_vst_pearson(raw_expr, theta=100, clip_value=10):
    """Pearson residual VST normalization (sctransform-style).

    Parameters
    ----------
    raw_expr : ndarray
        Raw expression values
    theta : float
        Dispersion parameter for negative binomial model
    clip_value : float
        Clip extreme residuals

    Returns
    -------
    ndarray
        Pearson residuals
    """
    mu = np.mean(raw_expr) + 1e-10
    variance = mu + (mu ** 2) / theta
    residuals = (raw_expr - mu) / np.sqrt(variance)
    return np.clip(residuals, -clip_value, clip_value)


def apply_vst_shifted_log(raw_expr, pseudocount=1.0):
    """Shifted log normalization.

    Ensures non-expressors get negative values.
    raw=0 maps to -0.5.

    Parameters
    ----------
    raw_expr : ndarray
        Raw expression values
    pseudocount : float
        Pseudocount added before log

    Returns
    -------
    ndarray
        Shifted log-normalized expression
    """
    log_expr = np.log(raw_expr + pseudocount)
    baseline = np.log(pseudocount)
    shift = baseline + 0.5
    return log_expr - shift


def generate_stochastic_response(n_total, receiver_indices, concentrations,
                                 basal_responsive, fold_change,
                                 sigma_r, sigma_r_basal, Kd,
                                 p_respond_max=0.9, hill_coef=1.0):
    """Generate stochastic response gene expression dependent on concentration.

    Hybrid model:
    1. B_j ~ Bernoulli(p_max * C^n / (Kd^n + C^n))
    2. If B_j=1: R_j = R_basal * (1 + FC * Act) * LN(0, sigma_r^2)
    3. If B_j=0: R_j = R_basal * LN(0, sigma_r_basal^2)

    Parameters
    ----------
    n_total : int
        Total number of cells
    receiver_indices : ndarray
        Indices of receiver cells
    concentrations : ndarray
        Concentration at each cell position
    basal_responsive : float
        Basal expression level
    fold_change : float
        Fold change constant
    sigma_r : float
        Lognormal sigma for responding cells
    sigma_r_basal : float
        Lognormal sigma for non-responding cells
    Kd : float
        Dissociation constant
    p_respond_max : float
        Maximum probability of responding
    hill_coef : float
        Hill coefficient

    Returns
    -------
    tuple
        (responsive_expr, responding_mask, response_probs)
    """
    n_receivers = len(receiver_indices)

    responsive_expr = basal_responsive * np.random.lognormal(0, sigma_r_basal, n_total)

    C = concentrations[receiver_indices]

    C_n = np.power(C, hill_coef)
    Kd_n = np.power(Kd, hill_coef)
    response_probs = p_respond_max * C_n / (Kd_n + C_n)

    responding_mask = np.random.rand(n_receivers) < response_probs
    n_responding = np.sum(responding_mask)

    non_responding_mask = ~responding_mask
    n_non_responding = np.sum(non_responding_mask)
    if n_non_responding > 0:
        non_responding_indices = receiver_indices[non_responding_mask]
        responsive_expr[non_responding_indices] = (
            basal_responsive * np.random.lognormal(0, sigma_r_basal, n_non_responding)
        )

    if n_responding > 0:
        responding_indices = receiver_indices[responding_mask]
        C_responding = C[responding_mask]

        activation = C_responding / (Kd + C_responding)
        mean_expr = basal_responsive * (1 + fold_change * activation)

        responsive_expr[responding_indices] = (
            mean_expr * np.random.lognormal(0, sigma_r, n_responding)
        )

    return responsive_expr, responding_mask, response_probs

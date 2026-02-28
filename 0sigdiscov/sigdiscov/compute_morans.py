#!/usr/bin/env python
"""
compute_morans.py - Spatial Moran's I Analysis for Cell-Cell Communication

Computes cross-correlation (normalized directional Moran's I) between factor genes
in sender cells and all target genes in receiver cells, across multiple radii and
cell type pairs.

Two output modes:
  matrix (default): HDF5 matrix file (pairs x radii x factors x targets) for production
  csv:              Per-pair CSV with permutation-based p-values for exploration

Requires GPU (CuPy) for computation.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import time
from datetime import datetime, timedelta
import anndata as ad
from scipy import sparse
import h5py
from tqdm import tqdm
import gc
import warnings
import logging
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional, Any
import threading
import gzip

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

@dataclass
class Config:
    """Analysis configuration with auto batch sizing."""

    RANDOM_SEED: int = 42
    RADII: List[float] = field(default_factory=lambda: [10, 20, 30, 50, 100, 200, 300, 500])
    MIN_EXPR_QUANTILE: float = 0.25
    MIN_CELLS_THRESHOLD: int = 10
    MIN_CONNECTIONS_THRESHOLD: int = 10

    # Precision
    USE_FLOAT16: bool = False
    FACTOR_BATCH_SIZE: Optional[int] = None  # Auto-set

    # Annular parameters
    USE_ANNULAR_EDGES: bool = False
    ANNULAR_WIDTH: float = 20.0

    # HDF5 output
    COMPRESSION: str = 'gzip'
    COMPRESSION_LEVEL: int = 1

    # GPU memory
    GPU_MEMORY_FRACTION: float = 0.8
    MIN_BATCH_SIZE: int = 50
    MAX_BATCH_SIZE: int = 2000
    GPU_AVAILABLE_MEM: float = 30e9
    CLEAR_GPU_EVERY_N: int = 10

    # CSV mode
    N_PERMUTATIONS: int = 10000
    FDR_ALPHA: float = 0.05
    GENE_BATCH_SIZE: int = 500

    def get_sigma(self, radius: float) -> float:
        return radius / 3.0

    def get_inner_radius(self, outer_radius: float) -> Optional[float]:
        if not self.USE_ANNULAR_EDGES:
            return None
        return max(0, outer_radius - self.ANNULAR_WIDTH)

    @staticmethod
    def generate_regular_radii(max_distance: float = 5000, spacing: float = 20) -> List[float]:
        return list(range(int(spacing), int(max_distance) + 1, int(spacing)))

    def calculate_optimal_batch_size(self, n_factors: int, n_cells: int,
                                     n_targets: int, available_gpu_gb: float,
                                     logger: logging.Logger) -> int:
        bytes_per_element = 2 if self.USE_FLOAT16 else 4
        memory_per_factor = (
            n_cells * bytes_per_element +
            n_cells +
            n_cells * bytes_per_element * 10
        ) / 1e9

        usable_memory = available_gpu_gb * self.GPU_MEMORY_FRACTION
        max_factors_in_memory = int(usable_memory / memory_per_factor)

        if n_factors <= self.MAX_BATCH_SIZE and n_factors <= max_factors_in_memory:
            optimal_batch = n_factors
            logger.info(f"All {n_factors} factors fit in memory - single batch")
        else:
            optimal_batch = min(max_factors_in_memory, self.MAX_BATCH_SIZE)
            optimal_batch = max(optimal_batch, self.MIN_BATCH_SIZE)
            n_batches = (n_factors + optimal_batch - 1) // optimal_batch
            logger.info(f"Will process {n_factors} factors in {n_batches} batches of {optimal_batch}")

        logger.info(f"Memory estimate: {memory_per_factor:.3f}GB per factor, "
                     f"{usable_memory:.1f}GB available")
        return optimal_batch


# ============================================================================
# LOGGING
# ============================================================================

def setup_logging(output_dir: Path, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger('compute_morans')
    logger.setLevel(logging.DEBUG)
    if logger.handlers:
        logger.handlers.clear()

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    fh = logging.FileHandler(output_dir / f'analysis_log_{timestamp}.log')
    fh.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)

    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ============================================================================
# GPU SETUP
# ============================================================================

def setup_gpu(gpu_id: str, gpu_type: str, config: Config, logger: logging.Logger):
    """Initialize GPU and return cupy module. Sets config.GPU_AVAILABLE_MEM."""
    os.environ['CUDA_VISIBLE_DEVICES'] = gpu_id
    try:
        import cupy as cp
        mempool = cp.get_default_memory_pool()

        if gpu_type == 'auto':
            try:
                mem_info = cp.cuda.runtime.memGetInfo()
                total_bytes = mem_info[1]
                gpu_gb = int(total_bytes / (1024**3) * 0.9)
                logger.info(f"Auto-detected GPU memory: {total_bytes/(1024**3):.1f}GB")
            except Exception:
                gpu_gb = 35
                logger.warning("Could not auto-detect GPU memory, using 35GB default")
        elif gpu_type.lower() == 'a100':
            gpu_gb = 75
            logger.info("GPU type: A100 (80GB total)")
        elif gpu_type.lower() == 'v100':
            gpu_gb = 35
            logger.info("GPU type: V100 (40GB total)")
        else:
            try:
                gpu_gb = int(float(gpu_type) * 0.9)
                logger.info(f"GPU memory specified: {gpu_type}GB total")
            except Exception:
                gpu_gb = 35
                logger.warning(f"Could not parse GPU type '{gpu_type}', using 35GB default")

        mempool.set_limit(size=gpu_gb * 1024**3)
        config.GPU_AVAILABLE_MEM = gpu_gb * 0.9 * 1e9
        logger.info(f"GPU initialized with {gpu_gb}GB memory limit")
        return cp, mempool

    except Exception as e:
        logger.error(f"GPU setup failed: {e}")
        sys.exit(1)


# ============================================================================
# FACTOR LIST READER
# ============================================================================

def read_factor_list(factor_file_path: str, logger: logging.Logger) -> List[str]:
    """Read factor gene names from TSV file (first row header, columns 1+)."""
    logger.info(f"Reading factor list from {factor_file_path}")
    try:
        if factor_file_path.endswith('.gz'):
            with gzip.open(factor_file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', nrows=1)
        else:
            df = pd.read_csv(factor_file_path, sep='\t', nrows=1)
        factor_names = df.columns[1:].tolist()
        logger.info(f"Found {len(factor_names)} factors in file")
        logger.info(f"First 5 factors: {factor_names[:5]}")
        return factor_names
    except Exception as e:
        logger.error(f"Error reading factor file: {e}")
        raise


# ============================================================================
# CELL PAIR PARSING
# ============================================================================

def parse_cell_pairs(cell_pairs_str: str) -> List[Tuple[str, str]]:
    """Parse cell pairs from 'sender1->receiver1,sender2->receiver2' string."""
    pairs = []
    for pair in cell_pairs_str.split(','):
        if '->' in pair:
            sender, receiver = pair.split('->', 1)
            pairs.append((sender.strip(), receiver.strip()))
    return pairs


def resolve_pairs(adata: ad.AnnData, cell_type_col: str, config: Config,
                  cell_pairs_str: Optional[str],
                  sender_types: Optional[List[str]],
                  receiver_types: Optional[List[str]],
                  all_pairs: bool,
                  logger: logging.Logger) -> Tuple[List[Tuple[str, str]], List[str]]:
    """Resolve cell type pairs from CLI args. Returns (pairs, cell_types_used)."""

    # Get valid cell types (enough cells)
    all_cts = []
    for ct in adata.obs[cell_type_col].unique():
        if (adata.obs[cell_type_col] == ct).sum() >= config.MIN_CELLS_THRESHOLD:
            all_cts.append(ct)
    logger.info(f"Valid cell types (>={config.MIN_CELLS_THRESHOLD} cells): {len(all_cts)}")

    if cell_pairs_str:
        pairs = parse_cell_pairs(cell_pairs_str)
        # Validate
        valid_pairs = []
        for s, r in pairs:
            if s in all_cts and r in all_cts:
                valid_pairs.append((s, r))
            else:
                logger.warning(f"Skipping pair {s}->{r}: cell type not found or too few cells")
        cell_types_used = list(set(s for s, r in valid_pairs) | set(r for s, r in valid_pairs))
        return valid_pairs, cell_types_used

    if sender_types and receiver_types:
        pairs = []
        for s in sender_types:
            for r in receiver_types:
                if s in all_cts and r in all_cts and s != r:
                    pairs.append((s, r))
        cell_types_used = list(set(sender_types) | set(receiver_types))
        return pairs, cell_types_used

    # Default: all non-self pairs
    pairs = [(s, r) for s in all_cts for r in all_cts if s != r]
    return pairs, all_cts


# ============================================================================
# WEIGHT MATRIX BUILDER
# ============================================================================

class WeightBuilder:
    """GPU-accelerated Gaussian weight matrix construction."""

    def __init__(self, cp, config: Config):
        self.cp = cp
        self.config = config

    def build(self, sender_coords_gpu, receiver_coords_gpu,
              radius: float, inner_radius: Optional[float] = None):
        """Build row-normalized sparse weight matrix. Returns (W_normalized, S0) or (None, 0)."""
        cp = self.cp
        n_senders = len(sender_coords_gpu)
        n_receivers = len(receiver_coords_gpu)

        dtype = cp.float16 if self.config.USE_FLOAT16 else cp.float32

        mem_per_chunk = n_receivers * 4 * (2 if self.config.USE_FLOAT16 else 4)
        available_mem = self.config.GPU_AVAILABLE_MEM
        chunk_size = min(n_senders, int(available_mem / mem_per_chunk / 4))

        sigma = self.config.get_sigma(radius)
        gaussian_factor = -1.0 / (2 * sigma * sigma)
        radius_sq = radius * radius
        inner_radius_sq = (inner_radius * inner_radius) if inner_radius else 0

        rows_list, cols_list, weights_list = [], [], []

        for chunk_start in range(0, n_senders, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_senders)
            sender_chunk = sender_coords_gpu[chunk_start:chunk_end].astype(dtype)

            diff_x = sender_chunk[:, 0:1] - receiver_coords_gpu[:, 0].astype(dtype).T
            diff_y = sender_chunk[:, 1:2] - receiver_coords_gpu[:, 1].astype(dtype).T
            dist_sq = diff_x * diff_x + diff_y * diff_y

            mask = dist_sq <= radius_sq
            if inner_radius:
                mask = mask & (dist_sq > inner_radius_sq)

            weights_chunk = cp.exp(dist_sq * gaussian_factor) * mask
            nonzero_mask = weights_chunk > 1e-6

            if cp.any(nonzero_mask):
                rows, cols = cp.where(nonzero_mask)
                rows_list.append(rows + chunk_start)
                cols_list.append(cols)
                weights_list.append(weights_chunk[nonzero_mask].astype(cp.float32))

            del diff_x, diff_y, dist_sq, mask, weights_chunk

        if rows_list:
            all_rows = cp.concatenate(rows_list)
            all_cols = cp.concatenate(cols_list)
            all_weights = cp.concatenate(weights_list)

            W_sparse = cp.sparse.coo_matrix(
                (all_weights, (all_rows, all_cols)),
                shape=(n_senders, n_receivers), dtype=cp.float32
            ).tocsr()

            row_sums = cp.array(W_sparse.sum(axis=1)).ravel()
            row_sums_inv = cp.where(row_sums > 0, 1.0 / row_sums, 0)
            D_inv = cp.sparse.diags(row_sums_inv, dtype=cp.float32, format='csr')
            W_normalized = D_inv @ W_sparse

            return W_normalized, float(W_sparse.sum())
        else:
            return None, 0.0


# ============================================================================
# MORAN'S I COMPUTATION
# ============================================================================

class MoransIComputer:
    """Batched Moran's I matrix computation on GPU."""

    def __init__(self, cp, config: Config):
        self.cp = cp
        self.config = config

    def compute_matrix(self, all_factor_expr_gpu, all_target_expr_gpu,
                       W_normalized, quantile_masks_gpu):
        """Compute Moran's I for all factor-target pairs at one radius/pair.
        Returns (n_factors, n_targets) matrix."""
        cp = self.cp
        n_factors = all_factor_expr_gpu.shape[1]
        n_targets = all_target_expr_gpu.shape[1]

        spatial_lags_all = W_normalized @ all_target_expr_gpu

        dtype = cp.float16 if self.config.USE_FLOAT16 else cp.float32
        morans_matrix = cp.zeros((n_factors, n_targets), dtype=dtype)

        batch_size = self.config.FACTOR_BATCH_SIZE or n_factors

        for f_start in range(0, n_factors, batch_size):
            f_end = min(f_start + batch_size, n_factors)
            for i in range(f_start, f_end):
                high_expr_mask = quantile_masks_gpu[i]
                if cp.sum(high_expr_mask) < self.config.MIN_CELLS_THRESHOLD:
                    continue

                high_expr_indices = cp.where(high_expr_mask)[0]
                factor_high = all_factor_expr_gpu[high_expr_indices, i]
                spatial_lags_high = spatial_lags_all[high_expr_indices, :]

                factor_norm = cp.linalg.norm(factor_high)
                if factor_norm > 0:
                    factor_normalized = factor_high / factor_norm
                    correlations = cp.dot(factor_normalized, spatial_lags_high)
                    spatial_norms = cp.linalg.norm(spatial_lags_high, axis=0)
                    valid_mask = spatial_norms > 0
                    correlations[valid_mask] /= spatial_norms[valid_mask]
                    morans_matrix[i, :] = correlations

        return morans_matrix


# ============================================================================
# HDF5 OUTPUT HANDLER
# ============================================================================

class HDF5OutputHandler:
    """Thread-safe HDF5 output for matrix mode."""

    def __init__(self, output_path: Path, config: Config):
        self.output_path = output_path
        self.config = config
        self.file = None
        self.lock = threading.Lock()

    def __enter__(self):
        self.file = h5py.File(self.output_path, 'w')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()

    def initialize(self, n_pairs: int, n_radii: int, n_factors: int, n_targets: int,
                   pairs: List[Tuple[str, str]], radii: List[float],
                   factor_names: List[str], target_names: List[str]):
        dtype = np.float16 if self.config.USE_FLOAT16 else np.float32

        # Chunk layout: (1, n_radii, n_factors, target_block) — optimized for
        # downstream pair-level reads ds[pidx, :, :, :] in extract_signatures.py.
        # Eliminates the need for rechunking (previously 15+ hours per quantile).
        self.file.create_dataset(
            'morans_i_matrices',
            shape=(n_pairs, n_radii, n_factors, n_targets),
            dtype=dtype,
            chunks=(1, n_radii, n_factors, min(100, n_targets)),
            compression=self.config.COMPRESSION,
            compression_opts=self.config.COMPRESSION_LEVEL
        )

        self.file.attrs['n_pairs'] = n_pairs
        self.file.attrs['n_radii'] = n_radii
        self.file.attrs['n_factors'] = n_factors
        self.file.attrs['n_targets'] = n_targets

        cell_types_set = list(dict.fromkeys(
            ct for pair in pairs for ct in pair
        ))
        self.file.create_dataset('cell_types', data=np.array(cell_types_set, dtype='S'))
        self.file.create_dataset('radii', data=radii)
        self.file.create_dataset('factor_genes', data=np.array(factor_names, dtype='S'))
        self.file.create_dataset('target_genes', data=np.array(target_names, dtype='S'))

        pair_strs = [f"{s}->{r}" for s, r in pairs]
        self.file.create_dataset('pairs', data=np.array(pair_strs, dtype='S'))

    def write_pair(self, pair_idx: int, pair_data: np.ndarray):
        """Write all radii for one pair at once. pair_data shape: (n_radii, n_factors, n_targets)."""
        with self.lock:
            self.file['morans_i_matrices'][pair_idx, :, :, :] = pair_data
            self.file.flush()


# ============================================================================
# PERMUTATION TESTING (CSV MODE)
# ============================================================================

def prepare_permutation_indices(n_cells: int, n_permutations: int,
                                random_seed: int, cp) -> Any:
    """Pre-generate permutation indices on GPU."""
    cp.random.seed(random_seed)
    n_cells = int(n_cells)
    dtype = cp.uint16 if n_cells < 65536 else cp.int32
    perm_indices = cp.empty((n_cells, n_permutations), dtype=dtype)
    batch_size = min(1000, n_permutations)
    for i in range(0, n_permutations, batch_size):
        end_i = min(i + batch_size, n_permutations)
        for j in range(i, end_i):
            perm_indices[:, j] = cp.random.permutation(n_cells).astype(dtype)
    return perm_indices


def compute_permutation_test(factor_expr_all: np.ndarray,
                             target_expr_receiver: np.ndarray,
                             W_normalized, W_sparse,
                             factor_gene: str,
                             target_gene_names: List[str],
                             min_expr_quantile: float,
                             n_permutations: int,
                             random_seed: int,
                             global_mean: float,
                             global_std: float,
                             target_means, target_stds,
                             radius: float,
                             inner_radius: Optional[float],
                             cp, logger) -> List[Dict]:
    """Run permutation test for one factor gene against all targets at one radius.

    Returns list of result dicts for CSV output.
    """
    n_senders = len(factor_expr_all)
    n_targets = len(target_gene_names)

    # Quantile filtering for observed statistic
    if min_expr_quantile > 0:
        threshold = np.quantile(factor_expr_all, min_expr_quantile)
        high_expr_indices = np.where(factor_expr_all > threshold)[0]
        n_high_expr = len(high_expr_indices)
        high_expr_mask = np.zeros(n_senders, dtype=bool)
        high_expr_mask[high_expr_indices] = True
        factor_expr_observed = factor_expr_all * high_expr_mask
    else:
        high_expr_mask = np.ones(n_senders, dtype=bool)
        factor_expr_observed = factor_expr_all
        n_high_expr = n_senders

    if factor_expr_observed.std() == 0:
        return []

    # Normalize factor
    factor_z_observed = (factor_expr_observed - global_mean) / (global_std + 1e-10)
    factor_z_all = (factor_expr_all - global_mean) / (global_std + 1e-10)

    factor_z_observed_gpu = cp.asarray(factor_z_observed, dtype=cp.float32)
    factor_z_all_gpu = cp.asarray(factor_z_all, dtype=cp.float32)

    # Normalize targets
    target_expr_gpu = cp.asarray(target_expr_receiver, dtype=cp.float32)
    target_means_gpu = cp.asarray(target_means, dtype=cp.float32)
    target_stds_gpu = cp.asarray(target_stds, dtype=cp.float32)
    target_stds_gpu = cp.maximum(target_stds_gpu, 1e-10)
    target_z_gpu = (target_expr_gpu - target_means_gpu) / target_stds_gpu

    # Spatial lags
    spatial_lags_norm = W_normalized @ target_z_gpu

    # Observed Moran's I
    covariances = cp.dot(factor_z_observed_gpu, spatial_lags_norm)
    norm_factor = cp.linalg.norm(factor_z_observed_gpu)
    norms_spatial_lag = cp.linalg.norm(spatial_lags_norm, axis=0)
    observed_I_all = covariances / (norm_factor * norms_spatial_lag + 1e-10)

    # Permutation test
    perm_indices = prepare_permutation_indices(n_senders, n_permutations, random_seed, cp)

    perm_I_all = cp.zeros((n_permutations, n_targets), dtype=cp.float32)

    # Random cell selections for each permutation
    cp.random.seed(random_seed + 1)
    random_selections = cp.zeros((n_senders, n_permutations), dtype=cp.bool_)
    for p in range(n_permutations):
        random_indices = cp.random.choice(n_senders, n_high_expr, replace=False)
        random_selections[random_indices, p] = True

    perm_batch_size = min(1000, n_permutations)
    for perm_start in range(0, n_permutations, perm_batch_size):
        perm_end = min(perm_start + perm_batch_size, n_permutations)
        batch_perms = perm_indices[:, perm_start:perm_end]

        for i in range(perm_end - perm_start):
            perm_idx = batch_perms[:, i]
            random_mask = random_selections[:, perm_start + i]
            perm_factor_z = factor_z_all_gpu[perm_idx]
            perm_factor_selected = perm_factor_z * random_mask

            perm_cov = cp.dot(perm_factor_selected, spatial_lags_norm)
            perm_norm = cp.linalg.norm(perm_factor_selected)

            if perm_norm > 0:
                perm_I_all[perm_start + i, :] = perm_cov / (perm_norm * norms_spatial_lag + 1e-10)

    # Statistics
    perm_means = cp.mean(perm_I_all, axis=0)
    perm_stds = cp.std(perm_I_all, axis=0) + 1e-10
    z_scores = (observed_I_all - perm_means) / perm_stds

    # Annulus info
    if inner_radius is not None:
        annulus_str = f"{inner_radius:.0f}-{radius:.0f}"
        edge_type = 'annular'
    else:
        annulus_str = f"0-{radius:.0f}"
        edge_type = 'circular'

    results = []
    for i in range(n_targets):
        if norms_spatial_lag[i] == 0:
            continue
        observed_I = float(observed_I_all[i])
        perm_I_gene = perm_I_all[:, i]
        n_extreme = int(cp.sum(cp.abs(perm_I_gene) >= cp.abs(observed_I)))
        p_value = (n_extreme + 1) / (n_permutations + 1)

        results.append({
            'factor_gene': factor_gene,
            'target_gene': target_gene_names[i],
            'interaction': f"{factor_gene}->{target_gene_names[i]}",
            'morans_i_normalized': observed_I,
            'zscore': float(z_scores[i]),
            'pvalue': p_value,
            'perm_mean': float(perm_means[i]),
            'perm_std': float(perm_stds[i]),
            'effect_size': float(z_scores[i]),
            'n_connections': int(W_normalized.nnz),
            'S0': float(W_sparse.sum()),
            'n_high_expr': n_high_expr,
            'n_total_senders': n_senders,
            'radius': radius,
            'inner_radius': inner_radius if inner_radius else 0,
            'annulus_range': annulus_str,
            'edge_type': edge_type,
        })

    return results


def apply_fdr_correction(results_df: pd.DataFrame, fdr_alpha: float = 0.05,
                         logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    """Apply BH FDR correction to p-values."""
    from statsmodels.stats.multitest import multipletests

    if logger:
        logger.info(f"Applying FDR correction (alpha={fdr_alpha})...")

    _, pvals_corrected, _, _ = multipletests(
        results_df['pvalue'].fillna(1.0), alpha=fdr_alpha, method='fdr_bh'
    )
    results_df['pvalue_fdr'] = pvals_corrected
    results_df['significant_fdr'] = pvals_corrected < fdr_alpha

    n_sig = results_df['significant_fdr'].sum()
    if logger:
        logger.info(f"FDR correction: {n_sig}/{len(results_df)} significant results")
    return results_df


# ============================================================================
# DATA PREPARATION (SHARED)
# ============================================================================

def prepare_data(adata: ad.AnnData, factor_names: Optional[List[str]],
                 config: Config, logger: logging.Logger):
    """Filter genes, normalize, extract factors. Returns dict of precomputed data."""
    n_cells = adata.shape[0]

    # Filter low-expression genes
    if sparse.issparse(adata.X):
        X_csc = adata.X.tocsc()
        gene_counts = np.diff(X_csc.indptr)
    else:
        gene_counts = (adata.X > 0).sum(axis=0)
    min_cells = max(10, n_cells * 0.001)
    keep_genes = gene_counts >= min_cells
    adata_filtered = adata[:, keep_genes].copy()
    n_genes = adata_filtered.shape[1]
    logger.info(f"Genes: {n_genes}/{adata.shape[1]} kept (min {min_cells:.0f} cells)")

    gene_names = adata_filtered.var_names.tolist()

    # Determine factors
    if factor_names is not None:
        factor_indices = []
        factor_names_found = []
        for factor in factor_names:
            if factor in gene_names:
                factor_indices.append(gene_names.index(factor))
                factor_names_found.append(factor)
        n_factors = len(factor_indices)
        logger.info(f"Found {n_factors}/{len(factor_names)} factors in gene list")
        if n_factors == 0:
            logger.error("No factors found in gene list!")
            sys.exit(1)
    else:
        factor_indices = list(range(n_genes))
        factor_names_found = gene_names
        n_factors = n_genes

    # Dense + normalize
    if sparse.issparse(adata_filtered.X):
        logger.info("Converting sparse to dense...")
        X_dense = adata_filtered.X.toarray().astype(np.float32)
    else:
        X_dense = adata_filtered.X.astype(np.float32)

    X_mean = X_dense.mean(axis=0, keepdims=True)
    X_std = X_dense.std(axis=0, keepdims=True)
    X_normalized = (X_dense - X_mean) / (X_std + 1e-10)

    if factor_names is not None:
        X_factors = X_dense[:, factor_indices]
        X_normalized_factors = X_normalized[:, factor_indices]
    else:
        X_factors = X_dense
        X_normalized_factors = X_normalized

    quantile_thresholds = np.quantile(X_factors, config.MIN_EXPR_QUANTILE, axis=0)

    return {
        'adata_filtered': adata_filtered,
        'gene_names': gene_names,
        'n_genes': n_genes,
        'factor_indices': factor_indices,
        'factor_names_found': factor_names_found,
        'n_factors': n_factors,
        'X_dense': X_dense,
        'X_normalized': X_normalized,
        'X_factors': X_factors,
        'X_normalized_factors': X_normalized_factors,
        'X_mean': X_mean.ravel(),
        'X_std': X_std.ravel(),
        'quantile_thresholds': quantile_thresholds,
    }


# ============================================================================
# MATRIX MODE PIPELINE
# ============================================================================

def run_matrix_mode(adata: ad.AnnData, pairs: List[Tuple[str, str]],
                    cell_type_col: str, config: Config,
                    cp, mempool, data: dict,
                    output_dir: Path, logger: logging.Logger) -> Dict[str, Any]:
    """Run matrix (HDF5) output mode."""

    n_pairs = len(pairs)
    n_factors = data['n_factors']
    n_genes = data['n_genes']
    radii = config.RADII
    n_radii = len(radii)

    # Auto batch size
    if config.FACTOR_BATCH_SIZE is None:
        config.FACTOR_BATCH_SIZE = config.calculate_optimal_batch_size(
            n_factors, adata.shape[0], n_genes,
            config.GPU_AVAILABLE_MEM / 1e9, logger
        )

    # Cache cell data
    cell_types_used = list(set(s for s, r in pairs) | set(r for s, r in pairs))
    cell_coords_cache = {}
    cell_expr_cache = {}
    cell_factor_cache = {}
    for ct in cell_types_used:
        mask = (adata.obs[cell_type_col] == ct).values
        coords = data['adata_filtered'].obsm['spatial'][mask].astype(np.float32)
        cell_coords_cache[ct] = coords
        cell_expr_cache[ct] = data['X_normalized'][mask, :]
        cell_factor_cache[ct] = data['X_normalized_factors'][mask, :]

    # Quantile thresholds on GPU
    quantile_thresholds_gpu = cp.asarray(data['quantile_thresholds'], dtype=cp.float32)

    # Initialize output
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = output_dir / f'morans_matrix_{timestamp}.h5'

    weight_builder = WeightBuilder(cp, config)
    morans_computer = MoransIComputer(cp, config)

    logger.info("=" * 80)
    logger.info(f"MATRIX MODE: {n_pairs} pairs x {n_radii} radii x {n_factors} factors x {n_genes} targets")
    logger.info(f"Total values: {n_pairs * n_radii * n_factors * n_genes:,}")
    logger.info("=" * 80)

    start_time = time.time()

    with HDF5OutputHandler(output_file, config) as h5:
        h5.initialize(n_pairs, n_radii, n_factors, n_genes,
                      pairs, radii,
                      data['factor_names_found'],
                      data['adata_filtered'].var_names.tolist())

        dtype = np.float16 if config.USE_FLOAT16 else np.float32
        total_ops = n_pairs * n_radii
        with tqdm(total=total_ops, desc="Computing matrices") as pbar:
            for pair_idx, (sender, receiver) in enumerate(pairs):
                # Buffer all radii for this pair, then write once (chunk-aligned)
                pair_buf = np.zeros((n_radii, n_factors, n_genes), dtype=dtype)

                sender_coords_gpu = cp.asarray(cell_coords_cache[sender])
                receiver_coords_gpu = cp.asarray(cell_coords_cache[receiver])
                sender_factors_gpu = cp.asarray(cell_factor_cache[sender], dtype=cp.float32)
                receiver_expr_gpu = cp.asarray(cell_expr_cache[receiver], dtype=cp.float32)

                masks = sender_factors_gpu > quantile_thresholds_gpu[cp.newaxis, :]
                masks = masks.T

                for radius_idx, radius in enumerate(radii):
                    inner_radius = config.get_inner_radius(radius) if config.USE_ANNULAR_EDGES else None

                    W_normalized, S0 = weight_builder.build(
                        sender_coords_gpu, receiver_coords_gpu, radius, inner_radius
                    )

                    if W_normalized is not None and W_normalized.nnz >= config.MIN_CONNECTIONS_THRESHOLD:
                        morans_matrix = morans_computer.compute_matrix(
                            sender_factors_gpu, receiver_expr_gpu, W_normalized, masks
                        )
                        pair_buf[radius_idx] = cp.asnumpy(morans_matrix)
                        del W_normalized

                    pbar.update(1)
                    elapsed = time.time() - start_time
                    done = radius_idx + 1 + pair_idx * n_radii
                    rate = done / elapsed
                    eta = (total_ops - done) / rate
                    pbar.set_postfix_str(f"P={pair_idx+1}/{n_pairs} R={radius:.0f}um | ETA={timedelta(seconds=int(eta))}")

                del sender_coords_gpu, receiver_coords_gpu, sender_factors_gpu, receiver_expr_gpu, masks
                h5.write_pair(pair_idx, pair_buf)
                del pair_buf

                if (pair_idx + 1) % config.CLEAR_GPU_EVERY_N == 0:
                    mempool.free_all_blocks()
                    cp._default_memory_pool.free_all_blocks()
                    gc.collect()

    runtime = (time.time() - start_time) / 60
    logger.info("=" * 80)
    logger.info(f"COMPLETE - Runtime: {runtime:.1f} min ({runtime/60:.1f} hrs)")
    logger.info(f"Output: {output_file}")
    logger.info(f"Shape: {n_pairs} x {n_radii} x {n_factors} x {n_genes}")
    logger.info("=" * 80)

    return {'success': True, 'output_file': str(output_file), 'runtime_minutes': runtime}


# ============================================================================
# CSV MODE PIPELINE
# ============================================================================

def run_csv_mode(adata: ad.AnnData, pairs: List[Tuple[str, str]],
                 cell_type_col: str, config: Config,
                 cp, mempool, data: dict,
                 output_dir: Path, logger: logging.Logger) -> Dict[str, Any]:
    """Run CSV output mode with permutation testing."""

    radii = config.RADII
    n_radii = len(radii)
    n_pairs = len(pairs)
    factor_names = data['factor_names_found']
    n_factors = data['n_factors']
    n_genes = data['n_genes']
    gene_names = data['gene_names']

    weight_builder = WeightBuilder(cp, config)

    logger.info("=" * 80)
    logger.info(f"CSV MODE: {n_pairs} pairs x {n_radii} radii x {n_factors} factors")
    logger.info(f"Permutations: {config.N_PERMUTATIONS}")
    logger.info("=" * 80)

    all_results = []
    start_time = time.time()

    total_ops = n_pairs * n_radii * n_factors
    with tqdm(total=total_ops, desc="Permutation testing") as pbar:
        for radius in radii:
            inner_radius = config.get_inner_radius(radius) if config.USE_ANNULAR_EDGES else None

            for sender, receiver in pairs:
                sender_mask = (adata.obs[cell_type_col] == sender).values
                receiver_mask = (adata.obs[cell_type_col] == receiver).values

                sender_coords = data['adata_filtered'].obsm['spatial'][sender_mask].astype(np.float32)
                receiver_coords = data['adata_filtered'].obsm['spatial'][receiver_mask].astype(np.float32)

                sender_coords_gpu = cp.asarray(sender_coords)
                receiver_coords_gpu = cp.asarray(receiver_coords)

                W_normalized, S0 = weight_builder.build(
                    sender_coords_gpu, receiver_coords_gpu, radius, inner_radius
                )

                if W_normalized is None or W_normalized.nnz < config.MIN_CONNECTIONS_THRESHOLD:
                    pbar.update(n_factors)
                    continue

                # Build sparse W for S0 (reuse from builder)
                W_sparse = W_normalized  # already have S0

                # Extract expression data
                if sparse.issparse(data['adata_filtered'].X):
                    X_csc = data['adata_filtered'].X.tocsc()
                    receiver_expr = X_csc[receiver_mask, :].toarray()
                else:
                    receiver_expr = data['adata_filtered'].X[receiver_mask, :]

                # Global stats for targets
                target_means = data['X_mean']
                target_stds = data['X_std']

                for factor_gene in factor_names:
                    factor_idx = gene_names.index(factor_gene)

                    if sparse.issparse(data['adata_filtered'].X):
                        factor_expr = data['adata_filtered'].X[sender_mask, factor_idx].toarray().ravel()
                    else:
                        factor_expr = data['adata_filtered'].X[sender_mask, factor_idx]

                    factor_mean = data['X_mean'][factor_idx]
                    factor_std = data['X_std'][factor_idx]

                    results = compute_permutation_test(
                        factor_expr_all=factor_expr,
                        target_expr_receiver=receiver_expr,
                        W_normalized=W_normalized,
                        W_sparse=W_normalized,  # for nnz/sum
                        factor_gene=factor_gene,
                        target_gene_names=gene_names,
                        min_expr_quantile=config.MIN_EXPR_QUANTILE,
                        n_permutations=config.N_PERMUTATIONS,
                        random_seed=config.RANDOM_SEED,
                        global_mean=factor_mean,
                        global_std=factor_std,
                        target_means=target_means,
                        target_stds=target_stds,
                        radius=radius,
                        inner_radius=inner_radius,
                        cp=cp,
                        logger=logger,
                    )

                    for r in results:
                        r['sender_type'] = sender
                        r['receiver_type'] = receiver
                    all_results.extend(results)

                    pbar.update(1)

                del W_normalized
                mempool.free_all_blocks()

    if not all_results:
        logger.warning("No results produced!")
        return {'success': False}

    results_df = pd.DataFrame(all_results)
    results_df = apply_fdr_correction(results_df, config.FDR_ALPHA, logger)

    # Save
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = output_dir / f'morans_results_{timestamp}.csv'
    results_df.to_csv(output_file, index=False)

    runtime = (time.time() - start_time) / 60
    logger.info("=" * 80)
    logger.info(f"COMPLETE - Runtime: {runtime:.1f} min")
    logger.info(f"Output: {output_file}")
    logger.info(f"Total results: {len(results_df)}")
    n_sig = results_df['significant_fdr'].sum()
    logger.info(f"Significant (FDR<{config.FDR_ALPHA}): {n_sig}")
    logger.info("=" * 80)

    return {'success': True, 'output_file': str(output_file), 'runtime_minutes': runtime}


# ============================================================================
# CLI
# ============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Spatial Moran\'s I Analysis for Cell-Cell Communication',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Production run (H5 matrix output)
  python compute_morans.py --data dataset/coad_spatial.h5ad --output results/out \\
    --factors dataset/db/factors.tsv.gz --max-distance 5000 --spacing 20

  # Single factor exploration
  python compute_morans.py --data dataset/coad_spatial.h5ad --output results/test \\
    --factor-genes TGFB1 --cell-pairs "macrophage->CAF" --radii 100 200 500

  # CSV mode with permutation testing
  python compute_morans.py --data dataset/coad_spatial.h5ad --output results/explore \\
    --output-format csv --factor-genes TGFB1 IFNG \\
    --cell-pairs "macrophage->CAF" --n-permutations 10000

  # Quick test
  python compute_morans.py --data dataset/coad_spatial.h5ad --output /tmp/test \\
    --factor-genes TGFB1 --cell-pairs "macrophage->CAF" --fast

  # List cell types
  python compute_morans.py --data dataset/coad_spatial.h5ad --list-cell-types
        """
    )

    # Required
    parser.add_argument('--data', type=str, required=True,
                        help='Path to h5ad file')
    parser.add_argument('--output', type=str, required=True,
                        help='Output directory')

    # Output control
    parser.add_argument('--output-format', type=str, choices=['matrix', 'csv'],
                        default='matrix',
                        help='Output format: matrix (H5) or csv with p-values (default: matrix)')

    # Cell type selection
    parser.add_argument('--cell-type-col', type=str, default='cell_type',
                        help='Column name for cell types (default: cell_type)')
    parser.add_argument('--cell-pairs', type=str, default=None,
                        help='Cell type pairs: "sender1->receiver1,sender2->receiver2"')
    parser.add_argument('--sender-types', type=str, nargs='+', default=None,
                        help='List of sender cell types')
    parser.add_argument('--receiver-types', type=str, nargs='+', default=None,
                        help='List of receiver cell types')
    parser.add_argument('--all-pairs', action='store_true',
                        help='All non-self pairs (default when no pairs specified)')

    # Factor selection
    parser.add_argument('--factors', type=str, default=None,
                        help='Path to factor TSV/TSV.gz file')
    parser.add_argument('--factor-genes', type=str, nargs='+', default=None,
                        help='Explicit factor gene names')

    # Spatial parameters
    parser.add_argument('--max-distance', type=float, default=5000,
                        help='Maximum distance for regular radii (default: 5000)')
    parser.add_argument('--spacing', type=float, default=20,
                        help='Spacing for regular radii (default: 20)')
    parser.add_argument('--radii', type=float, nargs='+', default=None,
                        help='Custom radii list (overrides max-distance/spacing)')
    parser.add_argument('--use-annular', action='store_true',
                        help='Enable annular (ring-shaped) edges')
    parser.add_argument('--annular-width', type=float, default=20,
                        help='Annular ring width (default: 20)')

    # Expression filtering
    parser.add_argument('--min-expr-quantile', type=float, default=0.25,
                        help='Minimum expression quantile for factor filtering (default: 0.25)')

    # Statistical (csv mode)
    parser.add_argument('--n-permutations', type=int, default=10000,
                        help='Number of permutations for csv mode (default: 10000)')
    parser.add_argument('--fdr-alpha', type=float, default=0.05,
                        help='FDR threshold (default: 0.05)')

    # Performance
    parser.add_argument('--gpu', type=str, default='0',
                        help='GPU device ID (default: 0)')
    parser.add_argument('--gpu-type', type=str, default='auto',
                        help='GPU type: auto, a100, v100, or memory in GB (default: auto)')
    parser.add_argument('--factor-batch-size', type=int, default=None,
                        help='Override auto batch size')
    parser.add_argument('--use-float16', action='store_true',
                        help='Use half precision (reduces memory, slight accuracy loss)')

    # Utility
    parser.add_argument('--list-cell-types', action='store_true',
                        help='List available cell types and exit')
    parser.add_argument('--fast', action='store_true',
                        help='Fast mode: reduced radii + small permutations for testing')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable debug logging')

    return parser


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = build_parser()
    args = parser.parse_args()

    data_path = Path(args.data)
    if not data_path.exists():
        print(f"Error: Data file not found: {data_path}")
        return 1

    # List cell types mode
    if args.list_cell_types:
        adata = ad.read_h5ad(data_path)
        col = args.cell_type_col
        if col not in adata.obs.columns:
            print(f"Error: Column '{col}' not found. Available: {list(adata.obs.columns)}")
            return 1
        cts = adata.obs[col].value_counts()
        print(f"\nCell types in '{col}' ({len(cts)} types, {adata.shape[0]} total cells):\n")
        for ct, count in cts.items():
            print(f"  {ct:40s}  {count:6d} cells")
        return 0

    # Setup output
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(output_dir, args.verbose)

    logger.info("=" * 80)
    logger.info("SPATIAL MORAN'S I ANALYSIS")
    logger.info("=" * 80)

    # Load data
    logger.info(f"Loading data from {data_path}...")
    adata = ad.read_h5ad(data_path)
    logger.info(f"Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")

    # Determine factors
    factor_names = None
    if args.factors:
        factor_names = read_factor_list(args.factors, logger)
    elif args.factor_genes:
        factor_names = args.factor_genes
        logger.info(f"Using {len(factor_names)} explicit factor genes: {factor_names}")

    # Configure
    config = Config()
    config.MIN_EXPR_QUANTILE = args.min_expr_quantile
    config.USE_ANNULAR_EDGES = args.use_annular
    config.ANNULAR_WIDTH = args.annular_width
    config.N_PERMUTATIONS = args.n_permutations
    config.FDR_ALPHA = args.fdr_alpha

    if args.use_float16:
        config.USE_FLOAT16 = True

    if args.factor_batch_size is not None:
        config.FACTOR_BATCH_SIZE = args.factor_batch_size

    # Radii
    if args.radii:
        config.RADII = args.radii
    else:
        config.RADII = Config.generate_regular_radii(args.max_distance, args.spacing)
        logger.info(f"Radii: {len(config.RADII)} distances from {args.spacing} to {args.max_distance}")

    # Fast mode overrides
    if args.fast:
        config.RADII = [100, 200, 500, 1000]
        config.N_PERMUTATIONS = 100
        logger.info("FAST MODE: radii=[100,200,500,1000], permutations=100")

    # Setup GPU
    cp, mempool = setup_gpu(args.gpu, args.gpu_type, config, logger)

    # Resolve cell pairs
    pairs, cell_types_used = resolve_pairs(
        adata, args.cell_type_col, config,
        args.cell_pairs, args.sender_types, args.receiver_types,
        args.all_pairs, logger
    )
    if not pairs:
        logger.error("No valid cell type pairs found!")
        return 1
    logger.info(f"Cell type pairs: {len(pairs)}")
    for s, r in pairs[:10]:
        logger.info(f"  {s} -> {r}")
    if len(pairs) > 10:
        logger.info(f"  ... and {len(pairs) - 10} more")

    # Prepare data
    data = prepare_data(adata, factor_names, config, logger)

    # Run
    if args.output_format == 'matrix':
        result = run_matrix_mode(adata, pairs, args.cell_type_col, config,
                                 cp, mempool, data, output_dir, logger)
    else:
        result = run_csv_mode(adata, pairs, args.cell_type_col, config,
                              cp, mempool, data, output_dir, logger)

    if result.get('success'):
        logger.info(f"Analysis completed successfully!")
        logger.info(f"Output: {result['output_file']}")
        return 0
    else:
        logger.error("Analysis failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

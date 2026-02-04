"""Batch processing utilities for multi-sample SPLISOSM analysis.

Provides parallel and serial processing of multiple spatial transcriptomics
samples, with checkpointing and resumption support.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any, Callable, Union
import os
import json

import numpy as np
from numpy.typing import NDArray

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

from .analysis.single_sample import SplisomAnalyzer, SampleResult
from .io.visium import load_visium_tsv, find_sample_files


@dataclass
class BatchConfig:
    """Configuration for batch processing.

    Attributes
    ----------
    k_neighbors : int
        ICAR kernel neighbors
    rho : float
        ICAR autocorrelation
    min_spots_expressed : int
        Minimum spots for gene testing
    use_gpu : bool
        Enable GPU acceleration
    n_workers : int
        Number of parallel workers (0 = serial)
    checkpoint_dir : str
        Directory for checkpoints
    save_intermediate : bool
        Save results after each sample
    """
    k_neighbors: int = 20
    rho: float = 0.99
    min_spots_expressed: int = 50
    use_gpu: bool = False
    n_workers: int = 0
    checkpoint_dir: Optional[str] = None
    save_intermediate: bool = True


def discover_samples(
    base_dir: str,
    pattern: str = '*/*_counts.tsv',
    recursive: bool = True
) -> List[Path]:
    """Find all Visium samples in a directory structure.

    Parameters
    ----------
    base_dir : str
        Root directory to search
    pattern : str
        Glob pattern for sample files
    recursive : bool
        Search recursively

    Returns
    -------
    samples : list
        List of Path objects to sample files
    """
    base_path = Path(base_dir)

    if recursive:
        samples = list(base_path.glob(f'**/{pattern}'))
    else:
        samples = list(base_path.glob(pattern))

    # Sort by path for reproducibility
    samples = sorted(samples)

    return samples


def sample_id_from_path(path: Path) -> str:
    """Extract sample ID from file path.

    Uses the pattern: STUDY_SAMPLEID_counts.tsv
    or falls back to filename stem.

    Parameters
    ----------
    path : Path
        Path to sample file

    Returns
    -------
    sample_id : str
        Extracted sample identifier
    """
    stem = path.stem  # filename without extension

    # Try to parse STUDY_SAMPLEID_counts pattern
    if stem.endswith('_counts'):
        stem = stem[:-7]  # Remove '_counts'

    # Include parent folder for uniqueness
    parent = path.parent.name
    if parent:
        return f"{parent}/{stem}"

    return stem


def _analyze_single_sample(
    sample_path: Path,
    config: BatchConfig,
    isoform_mapping: Optional[dict] = None
) -> Optional[SampleResult]:
    """Analyze a single sample (worker function).

    Parameters
    ----------
    sample_path : Path
        Path to sample file
    config : BatchConfig
        Analysis configuration
    isoform_mapping : dict, optional
        Pre-computed isoform mapping

    Returns
    -------
    result : SampleResult or None
        Analysis result, or None if failed
    """
    from .analysis.single_sample import build_isoform_mapping

    try:
        sample_id = sample_id_from_path(sample_path)

        # Load sample
        counts, coords, genes, spot_ids = load_visium_tsv(str(sample_path))

        # Build isoform mapping if not provided
        if isoform_mapping is None:
            isoform_mapping = build_isoform_mapping(genes)

        if len(isoform_mapping) == 0:
            return None

        # Create analyzer and run
        analyzer = SplisomAnalyzer(
            k_neighbors=config.k_neighbors,
            rho=config.rho,
            min_spots_expressed=config.min_spots_expressed,
            use_gpu=config.use_gpu
        )

        result = analyzer.analyze_sample(
            counts=counts,
            coords=coords,
            genes=genes,
            isoform_mapping=isoform_mapping,
            sample_id=sample_id,
            show_progress=False
        )

        return result

    except Exception as e:
        import warnings
        warnings.warn(f"Failed to analyze {sample_path}: {e}")
        return None


def analyze_samples_serial(
    sample_paths: List[Union[str, Path]],
    config: Optional[BatchConfig] = None,
    checkpoint_dir: Optional[str] = None,
    callback: Optional[Callable[[SampleResult], None]] = None,
    show_progress: bool = True
) -> List[SampleResult]:
    """Run SPLISOSM on multiple samples serially.

    Supports checkpointing for resumption after interruption.

    Parameters
    ----------
    sample_paths : list
        Paths to sample files
    config : BatchConfig, optional
        Analysis configuration
    checkpoint_dir : str, optional
        Directory for checkpoint files (enables resumption)
    callback : callable, optional
        Function called after each sample (receives SampleResult)
    show_progress : bool
        Show progress bar

    Returns
    -------
    results : list
        List of SampleResult objects
    """
    if config is None:
        config = BatchConfig()

    sample_paths = [Path(p) for p in sample_paths]
    results = []

    # Load checkpoint if exists
    completed_samples = set()
    if checkpoint_dir:
        checkpoint_path = Path(checkpoint_dir)
        checkpoint_path.mkdir(parents=True, exist_ok=True)
        checkpoint_file = checkpoint_path / 'completed.json'

        if checkpoint_file.exists():
            with open(checkpoint_file, 'r') as f:
                completed_samples = set(json.load(f))

    # Progress tracking
    if show_progress:
        try:
            from tqdm import tqdm
            sample_iter = tqdm(sample_paths, desc="Analyzing samples")
        except ImportError:
            sample_iter = sample_paths
    else:
        sample_iter = sample_paths

    for path in sample_iter:
        sample_id = sample_id_from_path(path)

        # Skip if already completed
        if sample_id in completed_samples:
            # Try to load cached result
            if checkpoint_dir:
                result_file = Path(checkpoint_dir) / f"{sample_id.replace('/', '_')}.json"
                if result_file.exists():
                    # Skip loading for now, just note completion
                    pass
            continue

        # Analyze sample
        result = _analyze_single_sample(path, config)

        if result is not None:
            results.append(result)

            # Save checkpoint
            if checkpoint_dir:
                completed_samples.add(sample_id)
                with open(Path(checkpoint_dir) / 'completed.json', 'w') as f:
                    json.dump(list(completed_samples), f)

                # Save individual result
                if config.save_intermediate:
                    result_file = Path(checkpoint_dir) / f"{sample_id.replace('/', '_')}.csv"
                    result.to_dataframe().to_csv(result_file, index=False)

            # Callback
            if callback:
                callback(result)

    return results


def analyze_samples_parallel(
    sample_paths: List[Union[str, Path]],
    config: Optional[BatchConfig] = None,
    n_workers: Optional[int] = None,
    show_progress: bool = True
) -> List[SampleResult]:
    """Run SPLISOSM on multiple samples in parallel.

    Uses multiprocessing for CPU parallelism. Note: GPU acceleration
    should use serial mode with use_gpu=True instead.

    Parameters
    ----------
    sample_paths : list
        Paths to sample files
    config : BatchConfig, optional
        Analysis configuration
    n_workers : int, optional
        Number of worker processes (default: CPU count)
    show_progress : bool
        Show progress bar

    Returns
    -------
    results : list
        List of SampleResult objects
    """
    from multiprocessing import Pool, cpu_count

    if config is None:
        config = BatchConfig()

    if n_workers is None:
        n_workers = cpu_count()

    # GPU not compatible with multiprocessing
    if config.use_gpu:
        import warnings
        warnings.warn(
            "GPU acceleration is not compatible with parallel processing. "
            "Disabling GPU. Use analyze_samples_serial() with use_gpu=True instead."
        )
        config.use_gpu = False

    sample_paths = [Path(p) for p in sample_paths]

    # Prepare arguments
    args = [(path, config) for path in sample_paths]

    # Run parallel
    results = []
    with Pool(n_workers) as pool:
        if show_progress:
            try:
                from tqdm import tqdm
                imap_iter = tqdm(
                    pool.starmap(_analyze_single_sample, args),
                    total=len(args),
                    desc=f"Analyzing samples ({n_workers} workers)"
                )
            except ImportError:
                imap_iter = pool.starmap(_analyze_single_sample, args)
        else:
            imap_iter = pool.starmap(_analyze_single_sample, args)

        for result in imap_iter:
            if result is not None:
                results.append(result)

    return results


def combine_results(
    results: List[SampleResult],
    include_sample_id: bool = True
) -> "pd.DataFrame":
    """Combine results from multiple samples into a single DataFrame.

    Parameters
    ----------
    results : list
        List of SampleResult objects
    include_sample_id : bool
        Include sample_id column

    Returns
    -------
    combined : DataFrame
        Combined results with all genes from all samples
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas required. Install with: pip install pandas")

    dfs = []
    for result in results:
        df = result.to_dataframe()
        if include_sample_id:
            df['sample_id'] = result.sample_id
        dfs.append(df)

    if len(dfs) == 0:
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)

    # Sort by gene for consistency
    combined = combined.sort_values(['gene', 'sample_id'] if include_sample_id else ['gene'])

    return combined.reset_index(drop=True)


def summarize_batch_results(
    results: List[SampleResult],
    fdr_threshold: float = 0.05
) -> Dict[str, Any]:
    """Generate summary statistics for batch analysis.

    Parameters
    ----------
    results : list
        List of SampleResult objects
    fdr_threshold : float
        FDR threshold for counting significant genes

    Returns
    -------
    summary : dict
        Summary statistics
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas required. Install with: pip install pandas")

    combined = combine_results(results)

    if len(combined) == 0:
        return {
            'n_samples': 0,
            'n_genes_total': 0,
            'n_genes_significant': 0,
            'genes_significant_multiple': []
        }

    # Apply BH correction across all samples
    from scipy import stats as scipy_stats

    n_tests = len(combined)
    sorted_pvals = np.sort(combined['pvalue_ir'].values)
    thresholds = fdr_threshold * np.arange(1, n_tests + 1) / n_tests

    if np.any(sorted_pvals <= thresholds):
        k = np.max(np.where(sorted_pvals <= thresholds)[0]) + 1
        sig_threshold = sorted_pvals[k - 1]
    else:
        sig_threshold = 0

    sig_mask = combined['pvalue_ir'] <= sig_threshold

    # Count significant genes per sample
    sig_genes = combined[sig_mask]['gene'].value_counts()
    genes_multiple = sig_genes[sig_genes > 1].index.tolist()

    return {
        'n_samples': len(results),
        'n_genes_total': combined['gene'].nunique(),
        'n_tests': n_tests,
        'n_significant': sig_mask.sum(),
        'significant_genes': combined[sig_mask]['gene'].unique().tolist(),
        'genes_significant_multiple': genes_multiple,
        'per_sample_stats': {
            r.sample_id: {
                'n_spots': r.n_spots,
                'n_genes_tested': r.n_genes_tested
            }
            for r in results
        }
    }

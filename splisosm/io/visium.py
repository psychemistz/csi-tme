"""Load Visium spatial transcriptomics data from TSV format.

Handles the specific format from:
/data/parks34/projects/0sigdiscov/moran_i/datasets/visium_large/preprocessData/

Format: TSV with spot coordinates as columns (`row x col` e.g., "1x1", "3x5"),
genes as rows.
"""

import re
from pathlib import Path
from typing import List, Tuple, Optional, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray


def parse_coordinate(coord_str: str) -> Tuple[int, int]:
    """Parse spot coordinate string to (row, col) tuple.

    Parameters
    ----------
    coord_str : str
        Coordinate string like "1x1", "3x5"

    Returns
    -------
    row, col : tuple of int
    """
    match = re.match(r'(\d+)x(\d+)', coord_str)
    if match:
        return int(match.group(1)), int(match.group(2))
    raise ValueError(f"Cannot parse coordinate: {coord_str}")


def load_visium_tsv(
    filepath: Union[str, Path],
    gene_column: Optional[str] = None
) -> Tuple[NDArray, NDArray, List, List]:
    """Load Visium data from TSV file.

    Parameters
    ----------
    filepath : str or Path
        Path to TSV file
    gene_column : str, optional
        Name of column containing gene names.
        If None, assumes first column or row index.

    Returns
    -------
    counts : ndarray
        Count matrix of shape (n_spots, n_genes)
    coords : ndarray
        Spatial coordinates of shape (n_spots, 2) as (row, col)
    genes : list
        Gene names
    spot_ids : list
        Spot coordinate strings
    """
    filepath = Path(filepath)

    # Read TSV
    df = pd.read_csv(filepath, sep='\t', index_col=0)

    # Gene names from index
    genes = df.index.tolist()

    # Spot coordinates from columns
    spot_ids = df.columns.tolist()
    coords = np.array([parse_coordinate(s) for s in spot_ids])

    # Count matrix: transpose so spots are rows
    counts = df.values.T.astype(np.float64)

    return counts, coords, genes, spot_ids


def load_multiple_samples(
    filepaths: List[Union[str, Path]],
    return_sample_ids: bool = True
) -> Tuple[List[NDArray], List[NDArray], List, List[List]]:
    """Load multiple Visium samples.

    Parameters
    ----------
    filepaths : list
        Paths to TSV files
    return_sample_ids : bool
        Whether to return sample IDs

    Returns
    -------
    counts_list : list of ndarray
        Count matrices for each sample
    coords_list : list of ndarray
        Coordinates for each sample
    genes : list
        Gene names (assumes same across samples)
    spot_ids_list : list of list
        Spot IDs for each sample
    """
    counts_list = []
    coords_list = []
    spot_ids_list = []
    genes = None

    for filepath in filepaths:
        counts, coords, sample_genes, spot_ids = load_visium_tsv(filepath)

        if genes is None:
            genes = sample_genes
        # Note: genes should match across samples in same dataset

        counts_list.append(counts)
        coords_list.append(coords)
        spot_ids_list.append(spot_ids)

    return counts_list, coords_list, genes, spot_ids_list


def extract_isoform_counts(
    counts: NDArray,
    genes: list,
    gene_name: str,
    isoform_pattern: Optional[str] = None
) -> Tuple[NDArray, list]:
    """Extract isoform counts for a specific gene.

    Parameters
    ----------
    counts : ndarray
        Count matrix of shape (n_spots, n_genes)
    genes : list
        Gene names
    gene_name : str
        Base gene name to extract
    isoform_pattern : str, optional
        Regex pattern for isoform naming.
        Default: matches gene_name followed by isoform suffix

    Returns
    -------
    isoform_counts : ndarray
        Counts of shape (n_spots, n_isoforms)
    isoform_names : list
        Names of matched isoforms
    """
    if isoform_pattern is None:
        # Match gene name with optional suffix (e.g., "VEGF", "VEGF-A", "VEGF_001")
        isoform_pattern = f'^{re.escape(gene_name)}(?:[-_.].*)?$'

    pattern = re.compile(isoform_pattern)

    # Find matching genes
    isoform_indices = []
    isoform_names = []
    for i, g in enumerate(genes):
        if pattern.match(g):
            isoform_indices.append(i)
            isoform_names.append(g)

    if not isoform_indices:
        raise ValueError(f"No isoforms found matching pattern for {gene_name}")

    isoform_counts = counts[:, isoform_indices]

    return isoform_counts, isoform_names


def find_sample_files(
    base_dir: Union[str, Path],
    study_pattern: str = '*',
    file_suffix: str = '_counts.tsv'
) -> List[Path]:
    """Find all sample files in a directory structure.

    Expected structure:
    base_dir/
        STUDY_NAME/
            SAMPLE_counts.tsv

    Parameters
    ----------
    base_dir : str or Path
        Base directory containing study folders
    study_pattern : str
        Glob pattern for study folders
    file_suffix : str
        Suffix for count files

    Returns
    -------
    filepaths : list of Path
    """
    base_dir = Path(base_dir)
    filepaths = sorted(base_dir.glob(f'{study_pattern}/*{file_suffix}'))
    return filepaths


def get_study_info(filepath: Union[str, Path]) -> dict:
    """Extract study and sample info from filepath.

    Parameters
    ----------
    filepath : str or Path
        Path like "base/BRCA_2021_Wu/1142243F_counts.tsv"

    Returns
    -------
    info : dict
        Contains 'study', 'sample_id', 'cancer_type'
    """
    filepath = Path(filepath)
    study = filepath.parent.name
    sample_id = filepath.stem.replace('_counts', '')

    # Extract cancer type from study name (e.g., "BRCA_2021_Wu" -> "BRCA")
    cancer_type = study.split('_')[0] if '_' in study else study

    return {
        'study': study,
        'sample_id': sample_id,
        'cancer_type': cancer_type,
        'filepath': str(filepath)
    }

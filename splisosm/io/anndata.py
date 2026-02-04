"""AnnData integration for SPLISOSM.

Provides conversion between SPLISOSM data structures and AnnData/scanpy
for interoperability with the broader single-cell/spatial ecosystem.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

try:
    import anndata as ad
    import scanpy as sc
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False


def check_anndata():
    """Check if AnnData is available."""
    if not ANNDATA_AVAILABLE:
        raise ImportError(
            "anndata and scanpy are required for AnnData integration. "
            "Install with: pip install anndata scanpy"
        )


def from_visium_tsv(
    counts: NDArray,
    coords: NDArray,
    genes: list,
    spot_ids: list
) -> "ad.AnnData":
    """Convert Visium TSV data to AnnData.

    Parameters
    ----------
    counts : ndarray
        Count matrix (n_spots, n_genes)
    coords : ndarray
        Spatial coordinates (n_spots, 2)
    genes : list
        Gene names
    spot_ids : list
        Spot identifiers

    Returns
    -------
    adata : AnnData
    """
    check_anndata()

    adata = ad.AnnData(
        X=counts,
        obs={'spot_id': spot_ids},
        var={'gene_name': genes}
    )

    # Add spatial coordinates
    adata.obsm['spatial'] = coords
    adata.obs['array_row'] = coords[:, 0]
    adata.obs['array_col'] = coords[:, 1]

    return adata


def to_splisosm_format(
    adata: "ad.AnnData",
    spatial_key: str = 'spatial'
) -> Tuple[NDArray, NDArray, list]:
    """Extract SPLISOSM-compatible arrays from AnnData.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    spatial_key : str
        Key in adata.obsm containing spatial coordinates

    Returns
    -------
    counts : ndarray
        Count matrix
    coords : ndarray
        Spatial coordinates
    genes : list
        Gene names
    """
    check_anndata()

    # Get count matrix
    if hasattr(adata.X, 'toarray'):
        counts = adata.X.toarray()
    else:
        counts = np.array(adata.X)

    # Get coordinates
    coords = adata.obsm[spatial_key]

    # Get gene names
    if 'gene_name' in adata.var.columns:
        genes = adata.var['gene_name'].tolist()
    else:
        genes = adata.var_names.tolist()

    return counts, coords, genes


def add_splisosm_results(
    adata: "ad.AnnData",
    results: dict,
    key_prefix: str = 'splisosm'
) -> None:
    """Add SPLISOSM results to AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data object (modified in place)
    results : dict
        Results dictionary with gene names as keys
    key_prefix : str
        Prefix for result keys in adata.var
    """
    check_anndata()

    # Initialize result columns
    n_genes = adata.n_vars
    pvalues = np.ones(n_genes)
    statistics = np.zeros(n_genes)

    gene_list = adata.var_names.tolist()

    for gene, result in results.items():
        if gene in gene_list:
            idx = gene_list.index(gene)
            pvalues[idx] = result.get('pvalue', 1.0)
            statistics[idx] = result.get('statistic', 0.0)

    adata.var[f'{key_prefix}_pvalue'] = pvalues
    adata.var[f'{key_prefix}_statistic'] = statistics


def load_10x_visium(
    path: str,
    count_file: str = 'filtered_feature_bc_matrix.h5'
) -> "ad.AnnData":
    """Load 10x Visium data using scanpy.

    Parameters
    ----------
    path : str
        Path to spaceranger output directory
    count_file : str
        Name of count file

    Returns
    -------
    adata : AnnData
    """
    check_anndata()

    adata = sc.read_visium(path, count_file=count_file)
    return adata

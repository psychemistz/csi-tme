"""Gene lists for CSI-TME cytokine/secreted protein analysis.

This module provides curated gene lists from the RESEARCH.md documentation
for analyzing cytokine isoform patterns in the tumor microenvironment.
"""

from typing import List, Dict, Set, Optional

# Cytokine and secreted protein targets for isoform analysis
CYTOKINE_TARGETS: Dict[str, List[str]] = {
    'interleukins': [
        'IL33',   # Full-length vs delta-exon3,4 (secretion competence)
        'IL15',   # LSP vs SSP (signal peptide variants)
        'IL4',    # Canonical vs IL-4delta2 (antagonist)
        'IL18',   # Pro-form vs mature (requires caspase-1)
        'IL37',   # Multiple isoforms with different activities
        'IL16',   # Pro-IL16 vs secreted
    ],
    'growth_factors': [
        'VEGFA',  # 3' UTR variants (detectable), internal splicing (limited)
        'TGFB1',  # TGF-beta family, gene-level ratios
        'TGFB2',
        'TGFB3',
        'EGF',    # Membrane-bound vs soluble precursor
        'FGF2',   # Multiple isoforms from alternative initiation
    ],
    'chemokines': [
        'CXCL12',  # alpha, beta, gamma (C-terminal, ECM affinity)
        'CCL2',    # MCP-1, APA variants
        'CCL5',    # RANTES, potential APA
        'CXCL8',   # IL-8, 3' UTR variants
    ],
    'death_ligands': [
        'TNFSF10',  # TRAIL: alpha vs short (dominant negative)
        'FASLG',    # FasL: membrane vs soluble
        'TNF',      # TNF-alpha: membrane vs soluble (ADAM17-dependent)
    ],
    'cytokine_inhibitors': [
        'IL18BP',   # IL-18 binding protein
        'IL1RN',    # IL-1 receptor antagonist (isoforms 1-4)
    ],
    'proteases': [
        'MMP2',      # Gelatinase A
        'MMP9',      # Gelatinase B
        'ADAMTS17',  # Secreted metalloprotease
        'CTSB',      # Cathepsin B
    ],
    'convertases': [
        'PCSK5',  # Proprotein convertase 5
        'PCSK6',  # Proprotein convertase 6 (PACE4)
        'FURIN',  # Ubiquitous convertase
    ],
}

# RNA-binding proteins that regulate cytokine splicing
RBP_REGULATORS: List[str] = [
    # FOX family (tissue-specific splicing)
    'RBFOX1',
    'RBFOX2',
    'RBFOX3',
    # CELF family (CUG-BP, ETR-3 like)
    'CELF1',
    'CELF2',
    'CELF4',
    'CELF5',
    # QKI (quaking)
    'QKI',
    # NOVA family (neuronal splicing)
    'NOVA1',
    'NOVA2',
    # PTB family (polypyrimidine tract binding)
    'PTBP1',
    'PTBP2',
    # ELAV/Hu family (AU-rich element binding)
    'ELAVL1',  # HuR
    'ELAVL2',  # HuB
    'ELAVL3',  # HuC
    'ELAVL4',  # HuD
    # SR proteins
    'SRSF1',
    'SRSF2',
    'SRSF3',
]

# Hypoxia markers for VEGF correlation analysis
HYPOXIA_MARKERS: List[str] = [
    'HIF1A',   # Hypoxia-inducible factor 1 alpha
    'LDHA',    # Lactate dehydrogenase A
    'CA9',     # Carbonic anhydrase IX
    'PGK1',    # Phosphoglycerate kinase 1
    'GLUT1',   # Glucose transporter 1 (SLC2A1)
    'SLC2A1',  # Alias for GLUT1
    'BNIP3',   # BCL2 interacting protein 3
    'LOX',     # Lysyl oxidase
]

# Stromal markers for tumor-stroma interface analysis
STROMAL_MARKERS: List[str] = [
    'ACTA2',   # Alpha-smooth muscle actin (CAFs)
    'FAP',     # Fibroblast activation protein
    'PDGFRB',  # PDGF receptor beta
    'COL1A1',  # Collagen type I
    'COL1A2',
    'VIM',     # Vimentin
    'S100A4',  # FSP1
    'THY1',    # CD90
]

# Immune cell markers for spatial immune context
IMMUNE_MARKERS: Dict[str, List[str]] = {
    't_cells': ['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'],
    'b_cells': ['CD19', 'CD79A', 'MS4A1'],  # MS4A1 = CD20
    'macrophages': ['CD68', 'CD163', 'MSR1', 'MRC1'],
    'dendritic': ['CD1C', 'CLEC9A', 'XCR1'],
    'nk_cells': ['NCAM1', 'NKG7', 'GNLY'],
}


def get_all_cytokines() -> List[str]:
    """Get all cytokine target genes as a flat list.

    Returns
    -------
    genes : list
        All cytokine genes from all categories
    """
    genes = []
    for category_genes in CYTOKINE_TARGETS.values():
        genes.extend(category_genes)
    return list(set(genes))  # Remove duplicates


def get_category_genes(category: str) -> List[str]:
    """Get genes for a specific category.

    Parameters
    ----------
    category : str
        Category name (e.g., 'interleukins', 'growth_factors')

    Returns
    -------
    genes : list
        Genes in the category

    Raises
    ------
    ValueError
        If category is not found
    """
    if category not in CYTOKINE_TARGETS:
        raise ValueError(
            f"Unknown category: {category}. "
            f"Available: {list(CYTOKINE_TARGETS.keys())}"
        )
    return CYTOKINE_TARGETS[category].copy()


def filter_genes_by_category(
    genes: List[str],
    category: Optional[str] = None,
    include_all_cytokines: bool = False
) -> List[str]:
    """Filter a gene list to specific categories.

    Parameters
    ----------
    genes : list
        Gene names to filter
    category : str, optional
        Specific category to include. If None, include all cytokines.
    include_all_cytokines : bool
        If True and category is None, include all cytokine targets

    Returns
    -------
    filtered : list
        Genes that match the category filter
    """
    if category is not None:
        target_genes = set(get_category_genes(category))
    elif include_all_cytokines:
        target_genes = set(get_all_cytokines())
    else:
        return genes  # No filtering

    return [g for g in genes if g in target_genes]


def get_rbp_regulators() -> List[str]:
    """Get all RBP regulator genes.

    Returns
    -------
    genes : list
        RNA-binding protein genes
    """
    return RBP_REGULATORS.copy()


def get_hypoxia_markers() -> List[str]:
    """Get hypoxia marker genes.

    Returns
    -------
    genes : list
        Hypoxia marker genes
    """
    return HYPOXIA_MARKERS.copy()


def get_stromal_markers() -> List[str]:
    """Get stromal/CAF marker genes.

    Returns
    -------
    genes : list
        Stromal marker genes
    """
    return STROMAL_MARKERS.copy()


def get_immune_markers(cell_type: Optional[str] = None) -> List[str]:
    """Get immune cell marker genes.

    Parameters
    ----------
    cell_type : str, optional
        Specific cell type (e.g., 't_cells', 'macrophages').
        If None, return all immune markers.

    Returns
    -------
    genes : list
        Immune marker genes
    """
    if cell_type is not None:
        if cell_type not in IMMUNE_MARKERS:
            raise ValueError(
                f"Unknown cell type: {cell_type}. "
                f"Available: {list(IMMUNE_MARKERS.keys())}"
            )
        return IMMUNE_MARKERS[cell_type].copy()

    genes = []
    for markers in IMMUNE_MARKERS.values():
        genes.extend(markers)
    return list(set(genes))


def genes_in_data(
    gene_list: List[str],
    available_genes: List[str]
) -> List[str]:
    """Filter gene list to those present in the data.

    Parameters
    ----------
    gene_list : list
        Desired genes
    available_genes : list
        Genes present in the dataset

    Returns
    -------
    present : list
        Genes from gene_list that are in available_genes
    """
    available_set = set(available_genes)
    return [g for g in gene_list if g in available_set]


def summarize_gene_lists() -> Dict[str, int]:
    """Get summary counts of all gene lists.

    Returns
    -------
    summary : dict
        Category names mapped to gene counts
    """
    summary = {}
    for category, genes in CYTOKINE_TARGETS.items():
        summary[f'cytokine_{category}'] = len(genes)
    summary['rbp_regulators'] = len(RBP_REGULATORS)
    summary['hypoxia_markers'] = len(HYPOXIA_MARKERS)
    summary['stromal_markers'] = len(STROMAL_MARKERS)
    summary['immune_markers'] = len(get_immune_markers())
    return summary

# SPLISOSM Datasets

This document describes the datasets used for CSI-TME analysis, sourced from the SPLISOSM paper and associated Zenodo repositories.

## Data Sources Overview

| Repository | DOI | Contents | Size |
|------------|-----|----------|------|
| Zenodo 15556390 | [10.5281/zenodo.15556390](https://zenodo.org/records/15556390) | Raw quantification data | 16.8 GB |
| Zenodo 16905935 | [10.5281/zenodo.16905935](https://zenodo.org/records/16905935) | Processed isoform-level data | 10.1 GB |

---

## Zenodo 15556390: Raw and Processed Quantification Data

Raw gene- and isoform-level quantification outputs (SpaceRanger, Sierra, etc.) for reproducing paper analyses.

### Files

| File | Size | Description | Data Source |
|------|------|-------------|-------------|
| `sit_nar_23.tar.gz` | 7.6 GB | SiT (Visium-ONT) mouse olfactory bulb & CBS | [GSE153859](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153859) |
| `gbm_visium_cell_24.tar.gz` | 6.3 GB | Human glioblastoma Visium-SR | [PRJNA994130](https://www.ncbi.nlm.nih.gov/sra/PRJNA994130) |
| `gbm_ont_nc_23.tar.gz` | 1.8 GB | Human glioma Visium-ONT | Ren et al. 2023 |
| `human_dlpfc.tar.gz` | 626.6 MB | Human DLPFC Visium-SR | [Globus](https://research.libd.org/globus/jhpce_HumanPilot10x/index.html) |
| `simulation_data.tar.gz` | 307.4 MB | Synthetic benchmark data | Generated |
| `visium_mouse_cbs.tar.gz` | 106.1 MB | Mouse CBS Visium-SR | [10X Genomics](https://www.10xgenomics.com/datasets/adult-mouse-brain-coronal-section-fresh-frozen-1-standard) |
| `slideseqv2_mouse_hippo.tar.gz` | 80.2 MB | Mouse hippocampus Slide-seqV2 | [Broad Institute](https://singlecell.broadinstitute.org/single_cell/study/SCP815) |

---

## Zenodo 16905935: Processed Isoform-Level Data

Processed spatial transcriptomics datasets with isoform-level quantification from multiple platforms (10X Visium-SR, Visium-ONT, Xenium Prime 5K). Isoform quantification uses Sierra for short-read data and direct outputs for long-read data.

### Sample Summary

**Mouse CBS (4 samples):**
- 2 Visium-ONT samples: CBS1, CBS2 (healthy adult coronal brain section)
- 1 Visium-SR sample: CBS (healthy adult)
- 1 Xenium Prime 5K sample: CBS (healthy adult)

**Human Brain (36 samples):**
- 12 Visium-SR DLPFC samples: 151507-151676 (postmortem healthy)
- 13 Visium-SR glioma samples: MGH258, ZH881inf, ZH881T1, ZH916bulk, ZH916inf, ZH916T1, ZH1007inf, ZH1007nec, ZH1019inf, ZH1019T1, ZH8811Abulk, ZH8811Bbulk, ZH8812bulk
- 11 Visium-ONT glioma samples: DMG_1-5, GBM_1-6

### Interactive Notebooks
- `sp_tx_diversity_human.ipynb` - Human dataset visualization
- `sp_tx_diversity_mouse_cbs.ipynb` - Mouse CBS visualization

---

## Local Data Structure

**Location:** `/data/parks34/projects/3csitme/data/`

```
data/
├── zenodo_15556390/                    # Raw quantification data (16.8 GB)
│   ├── README.md
│   ├── sit_nar_23.tar.gz
│   ├── gbm_visium_cell_24.tar.gz
│   ├── gbm_ont_nc_23.tar.gz
│   ├── human_dlpfc.tar.gz
│   ├── simulation_data.tar.gz
│   ├── visium_mouse_cbs.tar.gz
│   └── slideseqv2_mouse_hippo.tar.gz
│
└── zenodo_16905935/                    # Processed isoform-level data (10.1 GB)
    ├── README.md
    ├── sp_tx_diversity_human.ipynb
    ├── sp_tx_diversity_mouse_cbs.ipynb
    │
    ├── human_dlpfc/                    # 12 DLPFC samples
    │   ├── 151507/
    │   │   ├── iso.quant.h5ad
    │   │   ├── gene.quant.h5ad
    │   │   ├── sv_results.csv
    │   │   ├── du_results.csv
    │   │   └── svs.exon.bed
    │   ├── 151508/
    │   └── ... (151509-151676)
    │
    ├── human_glioma_ont/               # 11 ONT glioma samples
    │   ├── DMG_1/
    │   │   ├── iso.quant.h5ad
    │   │   ├── gene.quant.h5ad
    │   │   ├── sv_results.csv
    │   │   └── svs_iso.gtf
    │   ├── DMG_2/ ... DMG_5/
    │   └── GBM_1/ ... GBM_6/
    │
    ├── human_glioma_sr/                # 13 SR glioma samples
    │   ├── MGH258/
    │   ├── ZH881inf/
    │   └── ... (13 samples)
    │
    ├── mouse_cbs/                      # Mouse cerebellum
    │   ├── subset_ensembl_102.gtf.gz   # Reference annotation
    │   ├── mouse_sit/                  # Visium-ONT (SiT)
    │   │   ├── cbs1/
    │   │   ├── cbs2/
    │   │   └── cbs_svs.txid.gtf
    │   ├── mouse_visium/               # Visium-SR
    │   │   └── cbs/
    │   └── mouse_xenium_5k/            # Xenium Prime 5K
    │       └── cbs/
    │
    └── splisosm_test_results/          # Pre-computed test results
        ├── sv_res/                     # 41 spatial variability files
        └── du_res/                     # 15 differential usage files
```

---

## Per-Sample File Formats

### AnnData Files (`.h5ad`)

**`iso.quant.h5ad`** - Isoform-level quantification
| Layer | Description |
|-------|-------------|
| `counts` | Raw isoform counts |
| `log1p` | Log-normalized counts |
| `ratios_obs` | Observed isoform ratios |
| `ratios_smoothed` | Spatially smoothed ratios |

**`gene.quant.h5ad`** - Gene-level quantification
| Layer | Description |
|-------|-------------|
| `counts` | Raw gene counts |
| `log1p` | Log-normalized counts |
| `adata.obs` | Spot annotations (cell types, regions) |

### Results Files

**`sv_results.csv`** - Spatial variability test results

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `n_iso` | Number of isoforms |
| `pct_spot_on` | Fraction of spots with expression |
| `count_avg` | Average count per spot |
| `count_std` | Standard deviation of counts |
| `perplexity` | Isoform diversity measure |
| `major_ratio_avg` | Average ratio of dominant isoform |
| `pvalue_hsic-ir` | **HSIC-IR p-value** (isoform ratio variation) |
| `padj_hsic-ir` | BH-adjusted HSIC-IR p-value |
| `pvalue_hsic-ic` | HSIC-IC p-value (combined) |
| `padj_hsic-ic` | BH-adjusted HSIC-IC p-value |
| `pvalue_hsic-gc` | HSIC-GC p-value (gene expression) |
| `padj_hsic-gc` | BH-adjusted HSIC-GC p-value |
| `pvalue_spark-x` | SPARK-X p-value (comparison) |
| `padj_spark-x` | BH-adjusted SPARK-X p-value |

**`du_results.csv`** - Differential usage test results (subset of samples)

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `covariate` | Covariate gene (e.g., RBP) |
| `pvalue_hsic` | HSIC p-value for association |
| `pvalue_hsic-gp` | HSIC with GP regression |
| `pvalue_glm` | GLM-based p-value |
| `pvalue_glmm` | GLMM-based p-value |

### Isoform Structure Files

**`svs_iso.gtf`** (ONT samples) - Full-length isoform structures
- Spatially variable isoforms filtered by `padj_hsic-ir`
- Mouse: Mus_musculus.GRCm38.102 reference
- Human: Gencode v32 (hg38) + novel isoforms

**`svs.exon.bed`** (SR samples) - TREND region structures
- BED12 format for 3' transcript diversity events
- Genome: mm10 (mouse), hg38 (human)

---

## Dataset Details by Tissue

### Human DLPFC (Visium-SR)
- **Source:** spatialLIBD / Maynard et al. 2021
- **Samples:** 12 (151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676)
- **Tissue:** Postmortem healthy dorsolateral prefrontal cortex
- **Isoform method:** Sierra (TREND quantification)

### Human Glioma - ONT (Visium-ONT)
- **Source:** Ren et al. 2023 (Nature Communications)
- **Samples:** 11 (DMG_1-5, GBM_1-6)
- **Tissue:** Surgically removed glioma (5 diffuse midline glioma, 6 glioblastoma)
- **Isoform method:** Direct from long-read

### Human Glioma - SR (Visium-SR)
- **Source:** PRJNA994130 / Cell 2024
- **Samples:** 13 (MGH258, ZH881inf, ZH881T1, ZH916bulk, ZH916inf, ZH916T1, ZH1007inf, ZH1007nec, ZH1019inf, ZH1019T1, ZH8811Abulk, ZH8811Bbulk, ZH8812bulk)
- **Tissue:** Surgically removed human glioma
- **Isoform method:** Sierra (TREND quantification)

### Mouse Cerebellum
- **Source:** Fu et al. 2023 / 10X Genomics
- **Platforms:** Visium-ONT (SiT), Visium-SR, Xenium Prime 5K
- **Samples:** CBS1, CBS2 (ONT), CBS (SR), CBS (Xenium)
- **Tissue:** Healthy adult mouse coronal brain section
- **Isoform method:** Direct (ONT) / Sierra (SR)

### Mouse Hippocampus (Slide-seqV2)
- **Source:** Broad Institute SCP815
- **Samples:** Puck_200115_08
- **Tissue:** Adult mouse hippocampus
- **Platform:** Slide-seqV2 (short-read)

---

## Preprocessing Pipeline

### Short-Read Data (Visium-SR, Slide-seqV2)
1. SpaceRanger alignment and gene-level quantification
2. Sierra (v0.99.27) for transcript 3' end diversity (TREND)
3. SPLISOSM spatial variability testing

### Long-Read Data (Visium-ONT / SiT)
1. Original study pipelines for full-length isoform quantification
2. Direct isoform-level count matrices
3. SPLISOSM spatial variability testing

### Xenium Data
1. Custom scripts for per-codeword transcript density
2. Spatial binning to match Visium resolution
3. SPLISOSM analysis

---

## References

**SPLISOSM:**
- Su J, et al. A computational framework for mapping isoform landscape and regulatory mechanisms from spatial transcriptomics data. *bioRxiv* 2025. doi: [10.1101/2025.05.02.651907](https://doi.org/10.1101/2025.05.02.651907)
- GitHub: [JiayuSuPKU/SPLISOSM](https://github.com/JiayuSuPKU/SPLISOSM)
- Paper repo: [JiayuSuPKU/SPLISOSM_paper](https://github.com/JiayuSuPKU/SPLISOSM_paper)

**Original Data Sources:**
- Maynard et al. (2021) - spatialLIBD DLPFC. *Nat Neurosci* 24:425-436
- Ren et al. (2023) - Glioma Visium-ONT. *Nat Commun* 14:1-14
- Fu et al. (2023) - Mouse SiT long-read spatial
- 10X Genomics - Mouse CBS Visium reference dataset

**Tools:**
- Sierra: [VCCRI/Sierra](https://github.com/VCCRI/Sierra) - TREND quantification from 3' scRNA-seq

---

## License

All datasets are available under **Creative Commons Attribution 4.0 International (CC BY 4.0)**.

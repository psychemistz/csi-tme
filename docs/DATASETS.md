# SPLISOSM Datasets

This document describes the datasets used for CSI-TME analysis, sourced from the SPLISOSM paper and associated Zenodo repositories.

## Data Sources

### Zenodo Repository 1: Raw and Processed Data
**DOI:** [10.5281/zenodo.15556390](https://zenodo.org/records/15556390)

Raw gene- and isoform-level quantification data (outputs of 10X SpaceRanger, Sierra, etc.) for reproducing paper analyses.

| File | Size | Description |
|------|------|-------------|
| `sit_nar_23.tar.gz` | 7.6 GB | SIT NAR 2023 dataset |
| `gbm_visium_cell_24.tar.gz` | 6.3 GB | GBM Visium Cell 2024 dataset |
| `gbm_ont_nc_23.tar.gz` | 1.8 GB | GBM ONT Nature Communications 2023 |
| `human_dlpfc.tar.gz` | 626.6 MB | Human dorsolateral prefrontal cortex |
| `simulation_data.tar.gz` | 307.4 MB | Simulation data for benchmarking |
| `visium_mouse_cbs.tar.gz` | 106.1 MB | Mouse cerebellum Visium |
| `slideseqv2_mouse_hippo.tar.gz` | 80.2 MB | Mouse hippocampus Slide-seqV2 |

**Total:** 16.8 GB

---

### Zenodo Repository 2: Processed Isoform-Level Data
**DOI:** [10.5281/zenodo.16905935](https://zenodo.org/records/16905935)

Processed spatial transcriptomics datasets with isoform-level analysis from multiple platforms (10X Visium-SR, Visium-ONT, Xenium Prime 5K).

| File | Size | Description |
|------|------|-------------|
| `human_dlpfc.zip` | 5.5 GB | Human DLPFC processed data |
| `human_glioma_sr.zip` | 2.3 GB | Human glioma short-read (Visium) |
| `human_glioma_ont.zip` | 1.3 GB | Human glioma long-read (ONT) |
| `mouse_cbs.zip` | 858.2 MB | Mouse cerebellum |
| `splisosm_test_results.zip` | 60.0 MB | Per-sample test results |
| `sp_tx_diversity_human.ipynb` | 32.5 MB | Human analysis notebook |
| `sp_tx_diversity_mouse_cbs.ipynb` | 16.0 MB | Mouse CBS analysis notebook |

**Total:** 10.1 GB

---

## Local Data Structure

### Test Results (Downloaded)
**Location:** `/data/parks34/projects/3csitme/data/splisosm_test_results/`

Pre-computed SPLISOSM analysis results for all samples:

```
splisosm_test_results/
├── sv_res/          # Spatial variability results (41 files)
│   ├── ONT-*.sv_results.csv    # Long-read samples
│   └── SR-*.sv_results.csv     # Short-read samples
└── du_res/          # Differential usage results (15 files)
    ├── ONT-*.du_results.csv
    └── SR-*.du_results.csv
```

### Planned Data Locations
```
/data/parks34/projects/3csitme/data/
├── splisosm_test_results/     # ✓ Downloaded (60 MB)
├── zenodo_15556390/           # Pending - Raw data (16.8 GB)
└── zenodo_16905935/           # Pending - Processed data (10.1 GB)
```

---

## File Formats

### Spatial Variability Results (`sv_res/*.csv`)

Per-sample results testing isoform spatial patterns.

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
| `pvalue_spark-x` | SPARK-X p-value (comparison method) |
| `padj_spark-x` | BH-adjusted SPARK-X p-value |

### Differential Usage Results (`du_res/*.csv`)

Per-sample results testing isoform-covariate associations.

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `covariate` | Covariate gene (e.g., RBP) |
| `pvalue_hsic` | HSIC p-value for association |
| `pvalue_hsic-gp` | HSIC with GP regression |
| `pvalue_glm` | GLM-based p-value |
| `pvalue_glmm` | GLMM-based p-value |

---

## Sample Naming Convention

### Prefix Codes
- **ONT-**: Long-read sequencing (Oxford Nanopore Technology)
- **SR-**: Short-read sequencing (10X Visium)

### Sample Types
| Pattern | Species | Tissue |
|---------|---------|--------|
| `ONT-CBS*` | Mouse | Cerebellum |
| `ONT-DMG*` | Human | Diffuse midline glioma |
| `ONT-GBM*` | Human | Glioblastoma |
| `SR-15XXXX` | Human | DLPFC (Visium) |
| `SR-CBS` | Mouse | Cerebellum |

---

## Datasets by Tissue/Platform

### Human Brain (DLPFC)
- **Platform:** 10X Visium (short-read)
- **Samples:** 12 samples (SR-151507 to SR-151676)
- **Source:** Maynard et al. (spatialLIBD)
- **Isoform quantification:** Sierra

### Human Glioma
- **Platform:** Visium + ONT long-read
- **Samples:** 11 samples (6 GBM, 5 DMG)
- **Source:** Greenwald et al. 2024
- **Isoform quantification:** Direct from long-read

### Mouse Cerebellum
- **Platform:** Visium + ONT, Visium short-read
- **Samples:** 3 samples (ONT-CBS1, ONT-CBS2, SR-CBS)
- **Source:** Fu et al. 2023
- **Isoform quantification:** Direct from long-read / Sierra

---

## Preprocessing Pipeline

### Short-Read Data (Visium)
1. SpaceRanger alignment and quantification
2. Sierra (v0.99.27) for 3' transcript diversity
3. SPLISOSM spatial variability testing

### Long-Read Data (Visium + ONT)
1. Original study pipelines for isoform quantification
2. Direct isoform-level counts
3. SPLISOSM spatial variability testing

### Xenium Data
1. Custom scripts for per-codeword transcript density
2. Spatial binning
3. SPLISOSM analysis

---

## References

**SPLISOSM Paper:**
- Su J, et al. A computational framework for mapping isoform landscape and regulatory mechanisms from spatial transcriptomics data. *bioRxiv* 2025. doi: [10.1101/2025.05.02.651907](https://doi.org/10.1101/2025.05.02.651907)

**GitHub:**
- [https://github.com/JiayuSuPKU/SPLISOSM_paper](https://github.com/JiayuSuPKU/SPLISOSM_paper)

**Original Data Sources:**
- Maynard et al. (2021) - spatialLIBD DLPFC
- Greenwald et al. (2024) - Glioma Visium + ONT
- Fu et al. (2023) - Mouse cerebellum long-read spatial

---

## License

All datasets are available under **Creative Commons Attribution 4.0 International (CC BY 4.0)**.

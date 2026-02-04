# Session Handoff - 2026-02-04

## What Was Done This Session

### 1. Identified Critical Data Issue
- Original 118-study Visium data has **gene-level counts only** (not isoform-level)
- Cannot use for SPLISOSM HSIC-IR analysis which requires isoform quantification

### 2. Established New Workflow
- **Discovery**: Use Zenodo glioma data (`data/zenodo_16905935/`) which has `iso.quant.h5ad`
- **Validation**: Use new scRNA-seq dataset (Zenodo 17055113) with PSI values

### 3. Ran Pilot Analysis
- Tested 10 high-confidence SecAct genes across 13 glioma samples
- **Key findings**:
  - HMGB1: 11/13 samples significant (p < 4e-86)
  - SPP1: 9/11 samples significant (p < 1e-90)
  - VCAN: 8/13 samples significant
  - C3: 7/11 samples significant
- Classic cytokines (IL-*, IFN*, TNF*) NOT detected (too lowly expressed)
- Results saved: `data/zenodo_16905935/pilot_splisosm_results.csv`

### 4. Inventoried Available Targets
- CytoSig: Only 4/43 cytokines have isoform data in glioma
- SecAct: 108/1170 secreted proteins have isoform data
- High-value targets: HMGB1, SPP1, VCAN, C3, CLU, SPARCL1, APP

## Downloading Now

**Zenodo 17055113**: Single-cell RNA alternative splicing in glioma
- 27.8 GB dataset
- SMART-seq2 data with MARVEL PSI values
- 32,304 cells with cell type annotations
- For validating spatial isoform discoveries

## What To Do When Returning

### 1. Check Download Status
```bash
ls -la /data/parks34/projects/3csitme/data/
# Look for new zenodo folder or the downloaded zip file
```

### 2. Extract and Explore Validation Data
```bash
cd /data/parks34/projects/3csitme/data/
unzip zenodo-20250904.zip  # or whatever the file is named
```

### 3. Cross-Reference Genes
Check which of the pilot-validated genes (HMGB1, SPP1, VCAN, C3) have:
- PSI values in the scRNA-seq data
- Sufficient cell coverage for comparison

### 4. Run Validation Analysis
For each gene with spatial isoform variation:
- Identify cell types that correspond to spatial regions (tumor core vs periphery)
- Compare PSI values between cell types
- Check if direction of isoform preference matches spatial pattern

## Key Files

| File | Purpose |
|------|---------|
| `docs/RESEARCH.md` | Full status and workflow (START HERE) |
| `CLAUDE.md` | Quick reference for data locations |
| `data/zenodo_16905935/pilot_splisosm_results.csv` | Pilot discovery results |
| `splisosm/` | Implementation code |

## Validation Logic

```
SPATIAL DISCOVERY                    scRNA-seq VALIDATION
─────────────────                    ────────────────────

HMGB1 spatially variable      →      Check HMGB1 PSI in:
(p < 4e-86 in ZH8811Abulk)           - Tumor cells vs immune cells
                                     - Hypoxic vs normoxic regions

SPP1 at tumor boundary        →      Check SPP1 PSI in:
(p < 1e-90)                          - CAFs vs tumor cells
                                     - Stromal vs epithelial
```

## Questions Answered This Session

1. **Why can't we use original Visium data?** → Gene-level only, need transcript-level
2. **Which genes actually have isoform data?** → 108 SecAct genes, very few cytokines
3. **Do we have spatial isoform variation?** → Yes! Strong signal in HMGB1, SPP1, VCAN, C3
4. **How does SPLISOSM Figure 2 work?** → Simulation benchmarks showing well-calibrated p-values

# CLAUDE.md

## Project Overview

**CSI-TME: Cytokine & Secreted protein spatial Isoform distribution inference in the Tumor MicroEnvironment**

Applies SPLISOSM (Spatial Isoform Statistical Omnibus Model) to detect spatially variable isoform patterns in cancer cytokine/secreted proteins using spatial transcriptomics.

**Current state (2026-02-05)**: Discovery and validation complete. 7 TIER1 genes validated (SPP1, CLU, HMGB1, CD74, APP, SPARCL1, VCAN) + 19 TIER2 spatial-only. **Next step**: RBP regulator identification via conditional HSIC (DU test).

**Documentation**:
- `RESEARCH.md`: Implementation plan, current status, validation workflow
- `MATH.md`: HSIC, kernels, RKHS, theorems
- `BIOBG.md`: Biological rationale
- `docs/CRITICS.md`: Validation assessment and RBP hypotheses

## Data Location

### Discovery Data (Isoform-level)
`data/zenodo_16905935/`
- `human_glioma_sr/`: 13 Visium SR samples with `iso.quant.h5ad`
- `human_glioma_ont/`: 11 Visium-ONT samples with `iso.quant.h5ad`

### Validation Data (scRNA-seq)
`data/GSE182109/`: 32,304 SMART-seq2 cells with PSI values from MARVEL

### Signature Matrices
`/vf/users/parks34/projects/1ridgesig/SecActpy/secactpy/data/`:
- `CytoSig.tsv.gz`: 43 cytokine signatures
- `SecAct.tsv.gz`: 1,170 secreted protein signatures

## SPLISOSM Framework

### Three Test Statistics
| Statistic | Tests | Use Case |
|-----------|-------|----------|
| **HSIC-GC** | Gene expression spatial variation | SVG detection |
| **HSIC-IR** | Isoform ratio spatial variation | **Primary** — isoform switching |
| **HSIC-IC** | Combined expression + composition | Screening |

### Conditional HSIC (DU Test)
Tests whether isoform usage X correlates with RBP expression Z after removing spatial confounding:
1. Spatial residualization: ε_X = X - f(Y), ε_Z = Z - g(Y)
2. Test: HSIC(ε_X, ε_Z)

Implemented in `splisosm/statistics/conditional.py`.

### Platform Constraint
Visium 3' SR detects **APA** (3' UTR length), not internal splice isoforms. ONT captures full-length transcripts.

## RBP Regulator Analysis

### Target Genes
- **TIER1** (validated): SPP1, CLU, HMGB1, CD74, APP, SPARCL1, VCAN
- **Top TIER2** (spatial-only): SPARC, B2M, PTGDS, SERPINE2, NCAM1, C3

### RBP Candidates (19 genes)
RBFOX1/2/3, CELF1/2/4/5, QKI, NOVA1/2, PTBP1/2, ELAVL1/2/3/4, SRSF1/2/3

### Key Hypotheses (from `docs/CRITICS.md`)
| Gene | Candidate RBPs | Rationale |
|------|---------------|-----------|
| APP | RBFOX1/2, NOVA1/2 | Regulate neuronal exon 7/8 |
| VCAN | CELF1/2, QKI | ECM gene APA regulation |
| HMGB1 | ELAVL1, SRSF1/3 | 3' end processing |

## Key Parameters

| Parameter | Default | Tumor-specific |
|-----------|---------|----------------|
| k_neighbors | 20 | 10-15 |
| ρ | 0.99 | 0.95-0.99 |
| min_spots | 50 | 20-30 |
| FDR | 0.05 | 0.01 |

## Technology Stack

Python 3.9+, NumPy, SciPy, AnnData, scanpy, SPLISOSM

## Git Commit Guidelines

- Do NOT include `Co-Authored-By` lines in commit messages
- Keep commit messages concise and descriptive

## Key References

- Su J, et al. SPLISOSM. *Nat Biotechnol* 2026. doi: 10.1038/s41587-025-02965-6
- Gretton et al. (2007) — HSIC theory
- Zhu et al. (2021) — SPARK-X
- Liu et al. (2009) — Chi-square approximation

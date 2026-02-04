# SPLISOSM Adoption Plan: Characterizing Cytokine/Secreted Protein Isoforms in Tumor Microenvironment with Single-Cell RNA-seq Validation

## Executive Summary

This document provides a comprehensive adoption plan for applying SPLISOSM to analyze cytokine and secreted protein isoforms in cancer spatial transcriptomics, followed by validation using single-cell RNA-seq methods like SCALPEL. The SPLISOSM framework offers robust kernel-based statistical methods for detecting spatially variable isoform patterns, with significant potential for analyzing cytokine diversity in tumors.

**Critical Platform Constraint**: While standard Visium's 3' bias fundamentally limits internal splice variant detection, alternative polyadenylation events and certain C-terminal isoform differences remain detectable, creating viable pathways for cancer cytokine spatial isoformics.

---

## CURRENT STATUS (Updated 2026-02-04)

### Critical Data Availability Issue

**Problem**: The original 118-study Visium dataset (`/data/parks34/projects/0sigdiscov/moran_i/datasets/visium_large/`) contains **gene-level counts only**, not isoform/transcript-level quantification. SPLISOSM requires isoform-level data for HSIC-IR analysis.

### Revised Workflow

```
DISCOVERY (Zenodo Glioma Data)          VALIDATION (scRNA-seq)
─────────────────────────────          ────────────────────────

Zenodo 16905935:                        Zenodo 17055113:
- 13 Visium-SR glioma samples           - 32,304 cells (SMART-seq2)
- 11 Visium-ONT glioma samples          - GBM + IDH-mutant glioma
- Isoform quantification (iso.quant.h5ad) - PSI values (MARVEL)
                                        - Cell type annotations
        │                                       │
        ▼                                       ▼
   SPLISOSM HSIC-IR                      Check same splicing
   → Spatially variable                  patterns in target
   isoform usage                         cell types
```

### Available Data

| Dataset | Location | Type | Isoform Data? |
|---------|----------|------|---------------|
| Zenodo glioma (SR) | `data/zenodo_16905935/human_glioma_sr/` | Spatial | ✅ Yes |
| Zenodo glioma (ONT) | `data/zenodo_16905935/human_glioma_ont/` | Spatial | ✅ Yes |
| Zenodo DLPFC | `data/zenodo_16905935/human_dlpfc/` | Spatial | ✅ Yes |
| scRNA-seq validation | `data/GSE182109/` (downloading) | Single-cell | ✅ Yes (PSI) |
| Original Visium 118 studies | `/data/parks34/.../visium_large/` | Spatial | ❌ Gene-level only |

### Pilot Results (2026-02-04)

Ran SPLISOSM HSIC-IR on 13 glioma SR samples with high-confidence SecAct targets:

| Gene | Samples Tested | Significant (p<0.05) | Min p-value | Biological Role |
|------|----------------|---------------------|-------------|-----------------|
| **HMGB1** | 13 | **11 (85%)** | 4.0e-86 | Alarmin/DAMP |
| **SPP1** | 11 | **9 (82%)** | 1.1e-90 | Osteopontin, tumor-stroma |
| **VCAN** | 13 | **8 (62%)** | 3.7e-19 | Versican, ECM |
| **C3** | 11 | **7 (64%)** | 1.2e-11 | Complement |
| **IGFBP5** | 6 | **5 (83%)** | 1.3e-05 | IGF signaling |
| **APP** | 7 | **4 (57%)** | 5.7e-28 | Amyloid precursor |

**Key finding**: Classic cytokines (IL-*, IFN*, TNF*) are NOT detected at isoform level due to low expression in glioma. However, **TME-associated secreted proteins** (HMGB1, SPP1, VCAN, C3) show strong spatially variable isoform patterns.

Results saved: `data/zenodo_16905935/pilot_splisosm_results.csv`

### Revised Target Gene List

Based on actual data availability (SecAct genes with ≥2 isoforms, >500 UMIs):

```python
# HIGH-CONFIDENCE TARGETS (pilot validated)
pilot_validated = ['HMGB1', 'SPP1', 'VCAN', 'C3', 'IGFBP5', 'APP']

# ADDITIONAL SECACT GENES WITH ISOFORMS (108 total)
high_expression_targets = [
    'CLU',      # Clusterin - 17,661 UMIs
    'SPARCL1',  # ECM - 8,219 UMIs
    'ITM2B',    # Membrane protein - 5,032 UMIs
    'AGT',      # Angiotensinogen - 2,508 UMIs
    'EDIL3',    # ECM - 1,939 UMIs
    'DKK3',     # Wnt inhibitor - 1,853 UMIs
    'NCAM1',    # Neural adhesion - 1,823 UMIs
    'CTSB',     # Cathepsin B - 1,819 UMIs
    'SPARC',    # ECM - 1,584 UMIs
    'CCL2',     # Chemokine - 654 UMIs (only cytokine with coverage)
]
```

### Next Steps (When Returning)

1. **Check download completion**: `ls -la data/GSE182109/` or new validation dataset
2. **Cross-reference genes**: Check which pilot-validated genes have PSI values in scRNA-seq
3. **Run validation analysis**: Compare isoform usage in cell types that map to spatial regions
4. **Expand discovery**: Run SPLISOSM on full 108 SecAct genes with isoforms

### Validation Logic

```
Spatial Discovery                    scRNA-seq Validation
───────────────                      ────────────────────

HMGB1 isoform A enriched      →      Do tumor cells (core-enriched)
at tumor core                        show higher PSI for same exon
                                     vs immune cells (periphery)?

SPP1 spatial gradient         →      Do CAFs show different SPP1
at tumor-stroma boundary             isoform usage vs tumor cells?
```

---

## 1. Core SPLISOSM Statistical Architecture

### 1.1 HSIC Kernel Independence Testing Framework

SPLISOSM reformulates spatial isoform analysis as multivariate independence testing using the Hilbert-Schmidt Independence Criterion (HSIC). For *n* spatial locations with isoform compositions X ∈ [0,1]^q and spatial coordinates Y ∈ ℝ², the framework tests whether X and Y are statistically independent through kernel embeddings.

The fundamental HSIC statistic measures dependence via the Hilbert-Schmidt norm of the cross-covariance operator between feature spaces:

$$\text{HSIC}(X,Y) := \|C_{XY}\|_{\text{HS}}^2$$

The empirical estimator takes the elegant form:

$$T_{\text{HSIC}}(X,Y) = \frac{1}{n}\text{tr}(K_X H L_Y H)$$

where **K_X** is the isoform kernel matrix, **L_Y** is the spatial kernel, and H = I − (1/n)11^T is the centering matrix. Under the null hypothesis of independence, this statistic follows an asymptotic chi-square mixture distribution with weights derived from the eigenvalues {λᵢ} and {μⱼ} of the centered kernel matrices:

$$T_{\text{HSIC}} \xrightarrow{d} \sum_{i,j=1}^{n} \lambda_i \mu_j z_{ij}^2$$

This eliminates the need for computationally expensive permutation testing—a critical advantage for large-scale analysis.

### 1.2 ICAR Spatial Kernel: Theoretically Optimal Spatial Dependency Modeling

SPLISOSM's key innovation is the Intrinsic Conditional Autoregressive (ICAR) spatial kernel based on Graph Fourier Transform principles:

$$L_Y = (D - \rho W)^{-1}$$

where W is the binary adjacency matrix from a k-nearest-neighbor spatial graph, D is the degree matrix with entries D_{ii} = Σⱼ Wᵢⱼ, and ρ ∈ [0,1) controls spatial autocorrelation strength.

The connection to Graph Fourier Transform provides crucial interpretability: low eigenvalues correspond to smooth, large-scale spatial patterns and high eigenvalues represent local fluctuations. SPLISOSM's Theorem 1 proves that any low-rank approximation (as used in SPARK-X) inevitably sacrifices statistical power because truncating high-frequency eigenvectors may discard genuine biological signal.

### 1.3 Three Complementary Test Statistics

**HSIC-GC (Gene Counts)**: Tests whether total gene expression varies spatially—equivalent to spatially variable expression (SVE) detection.

**HSIC-IR (Isoform Ratios)**: Tests whether isoform usage proportions vary spatially, independent of total expression changes. This specifically identifies spatially variable processing (SVP) genes where alternative splicing or polyadenylation patterns change across tissue regions. **This is optimal for cytokine analysis.**

**HSIC-IC (Isoform Counts)**: Detects either isoform usage or expression changes, providing a combined measure of spatial isoform diversity.

**Key biological insight**: HSIC-IR ignores total expression changes, while HSIC-IC captures both types of variation. For cytokine analysis, HSIC-IR is optimal for detecting microenvironment-driven isoform switching independent of overall cytokine upregulation.

### 1.4 Compositional Data Handling

Spatial transcriptomics exhibits extreme sparsity at the isoform level—many spots have zero UMI counts, making isoform ratios undefined. SPLISOSM's Theorem 2 provides a mathematically rigorous solution through a zero-padded centered kernel with mean replacement:

$$K_{X_n} = \begin{pmatrix} H_m K_{X_m} H_m & 0 \\ 0 & 0 \end{pmatrix}$$

This is equivalent to replacing undefined ratios with global average isoform proportions before kernel computation—a strategy that preserves test validity while maintaining computational efficiency.

---

## 2. Platform Constraints: Understanding Visium's 3' Capture Bias

The critical challenge for cytokine isoform analysis is Visium's 3' capture bias. Standard Visium uses poly-dT priming that captures only ~500-1000 bp from the 3' end, fundamentally limiting detection of internal splice variants.

### 2.1 VEGF-A Isoform Detection Feasibility

| Isoform | Exon Structure | 3' Visium Detectable? | Rationale |
|---------|---------------|----------------------|-----------|
| VEGF121 vs VEGF165 vs VEGF189 | Differ in exons 6/7 (internal) | **No** | Same terminal exon 8a; differences are internal |
| VEGF165 vs VEGF165b | Exon 8a vs 8b (terminal) | **Potentially yes** | Different 3' terminal exon usage |
| CXCL12α vs CXCL12β vs CXCL12γ | C-terminal differences | **Potentially yes** | Isoforms differ near 3' end |

### 2.2 Recommended Cytokine Targets for 3' Visium

Based on platform constraints, prioritize:

1. **Alternative polyadenylation** in IL-6 family receptors (gp130/IL6ST has documented APA)
2. **TGF-β1/β2/β3 discrimination** (different genes, fully distinguishable)
3. **CXCL12 isoforms** (C-terminal differences detectable)
4. **CCL21 length variants** (extended C-terminus)
5. **3' UTR length variation** via TREND analysis

For comprehensive VEGF isoform analysis, long-read spatial transcriptomics (Nanopore + Visium, SiT protocol) is required to resolve exon 6/7 alternative splicing.

---

## 3. Target Gene Selection for Tumor Microenvironment Analysis

### 3.1 High-Priority Cytokine Targets

```python
cytokine_targets = {
    # Cytokines with splicing-regulated activity
    'interleukins': ['IL33', 'IL15', 'IL4', 'IL18', 'IL37', 'IL16'],
    'growth_factors': ['VEGFA', 'TGFB1', 'TGFB2', 'TGFB3', 'EGF', 'FGF2'],
    'chemokines': ['CXCL12', 'CCL2', 'CCL5', 'CXCL8'],
    'death_ligands': ['TNFSF10', 'FASLG', 'TNF'],  # TRAIL, FasL, TNF-alpha
    
    # Secreted inhibitors
    'cytokine_inhibitors': ['IL18BP', 'IL1RN'],  # IL-18BP, IL-1Ra
    
    # Proteases and regulators
    'proteases': ['MMP2', 'MMP9', 'ADAMTS17', 'CTSB'],
    'convertases': ['PCSK5', 'PCSK6', 'FURIN'],
}
```

### 3.2 RBP Regulators for Differential Association Testing

```python
rbp_regulators = [
    'RBFOX1', 'RBFOX2', 'RBFOX3',  # RBFOX family
    'CELF1', 'CELF2', 'CELF4', 'CELF5',  # CELF family
    'QKI',  # Quaking
    'NOVA1', 'NOVA2',  # NOVA family
    'PTBP1', 'PTBP2',  # PTB proteins
    'ELAVL1', 'ELAVL2', 'ELAVL3', 'ELAVL4',  # HuR family
    'SRSF1', 'SRSF2', 'SRSF3',  # SR proteins
]
```

---

## 4. Technical Implementation for Pan-Cancer Analysis

### 4.1 Computational Requirements

| Component | Memory | Computation Time |
|-----------|--------|------------------|
| Spatial kernel per sample | ~60 MB (sparse) | 2-5 seconds |
| HSIC-IR test per gene | ~100 KB working | 0.1 seconds |
| Full analysis (10K genes) | 8-16 GB per sample | 20-40 minutes |
| Batch integration (1000 samples) | 100+ GB | 12-48 hours |

### 4.2 Analysis Pipeline

```python
import scanpy as sc
import pandas as pd
from splisosm import SplisosmNP
from statsmodels.stats.multitest import multipletests

# === STEP 1: Load and preprocess data ===
adata = sc.read_h5ad('tumor_visium_isoforms.h5ad')

# Filter for cytokine targets
cytokine_genes = [g for category in cytokine_targets.values() for g in category]
cytokine_isoforms = adata.var[adata.var['gene_name'].isin(cytokine_genes)].index
adata_cytokines = adata[:, cytokine_isoforms].copy()

# === STEP 2: Spatial variability testing ===
model = SplisosmNP(adata_cytokines, spatial_key='spatial')

sv_results = model.test_spatial_variability(
    test_type='HSIC-IR',  # Focus on isoform ratios
    kernel_type='icar',
    rho=0.99,
    k=4  # 4 mutual nearest neighbors
)

# Identify SVP cytokines
sv_results['padj'] = multipletests(sv_results['pvalue'], method='fdr_bh')[1]
svp_cytokines = sv_results[sv_results['padj'] < 0.05]['gene_name'].unique()
print(f"Spatially variable cytokines: {svp_cytokines}")

# === STEP 3: Differential usage with RBPs ===
# Add RBP expression as covariates
for rbp in rbp_regulators:
    if rbp in adata.var['gene_name'].values:
        rbp_idx = adata.var[adata.var['gene_name'] == rbp].index[0]
        adata_cytokines.obs[f'{rbp}_expr'] = adata[:, rbp_idx].X.toarray().flatten()

du_results = model.test_differential_usage(
    genes=svp_cytokines,
    covariates=[f'{rbp}_expr' for rbp in rbp_regulators 
                if f'{rbp}_expr' in adata_cytokines.obs],
    conditional=True,  # Use conditional HSIC to control spatial confounding
    method='hsic'
)

# === STEP 4: Identify regulatory relationships ===
du_results['padj'] = multipletests(du_results['pvalue'], method='fdr_bh')[1]
significant_pairs = du_results[du_results['padj'] < 0.05]

print("Significant RBP-cytokine isoform associations:")
print(significant_pairs[['gene', 'covariate', 'pvalue', 'padj']])
```

### 4.3 Parameter Modifications for Tumor Microenvironment

**Spatial kernel bandwidth**: Tumor microenvironments exhibit high cellular heterogeneity at fine scales compared to structured neural tissue:
- **Tumor core/stroma**: k = 10-15 neighbors (narrow bandwidth for niche-specific patterns)
- **Leading edge**: k = 5-10 (very fine scale for invasion zone gradients)
- **Neural tissue comparison**: k = 20-50 (broader for layer-level patterns)

**Sparsity thresholds**: Cytokines are often lowly expressed:
- Minimum spots with nonzero counts: m ≥ 20 (vs 50+ for abundant genes)
- Consider gene-level aggregation before isoform ratio computation
- Apply KNN imputation for sparse APA signals

**Multiple testing correction**: For 1000 samples × ~100 cytokine genes × 4 spatial regions:
- ~400,000 tests requiring stringent FDR control
- Hierarchical correction: within-sample BH followed by meta-analysis
- Consider q-value approach with empirical null estimation

---

## 5. Single-Cell Validation Using SCALPEL

### 5.1 SCALPEL Overview

SCALPEL (Single-Cell ALternative PoLyadenylation and Expression quantification using Long-reads) provides a Nextflow-based pipeline for transcript isoform quantification using 3'-tagged single-cell RNA-seq data. It enables validation of spatial isoform patterns detected by SPLISOSM.

**Full Citation**: Ake F, Schilling M, Fernández-Moya SM, et al. Quantification of transcript isoforms at the single-cell level using SCALPEL. *Nature Communications*. 2025;16:6402. doi: 10.1038/s41467-025-61118-0

**Key Resources**:
- GitHub: https://github.com/plasslab/SCALPEL
- Zenodo container: doi: 10.5281/zenodo.15717636

### 5.2 SCALPEL Methodological Approach

SCALPEL uses a hybrid annotation approach combining:
1. **Reference transcriptome** (Gencode/Ensembl annotations)
2. **De novo 3' end detection** from clustered read terminations
3. **EM algorithm** for isoform-level abundance estimation

The truncation parameter (default: 600 nucleotides from 3' end) defines the analysis window, accounting for 3' bias in typical scRNA-seq protocols.

### 5.3 Validation Strategy: Integrating SCALPEL with SPLISOSM

**Rationale**: SCALPEL provides cell-type-level isoform resolution, while SPLISOSM detects spatial patterns. If a gene shows spatially variable isoform usage (SVP) in SPLISOSM, we expect:
- Cell types enriched in different spatial regions should show differential isoform usage in SCALPEL
- The direction of isoform preference should match the spatial pattern

```python
# Workflow: Validate spatial isoform patterns with single-cell data

import pandas as pd
import scanpy as sc
import numpy as np
from itertools import combinations
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# 1. Load SPLISOSM results
splisosm_svp = pd.read_csv('splisosm_sv_results.csv')
svp_genes = splisosm_svp[splisosm_svp['padj'] < 0.05]['gene'].unique()
print(f"Found {len(svp_genes)} SVP genes from SPLISOSM")

# 2. Load SCALPEL results
# Option A: Direct from APADGE matrix
isoform_counts = pd.read_csv('sample_APADGE.txt', sep='\t', index_col=0)
scalpel_adata = sc.AnnData(isoform_counts.T)

# Option B: From converted Seurat object
# scalpel_adata = sc.read_h5ad('scalpel_isoforms.h5ad')

# 3. Add gene-isoform mapping
scalpel_adata.var['gene'] = scalpel_adata.var_names.str.split('-').str[0]

# 4. Define isoform ratio computation
def compute_isoform_ratios(adata, gene):
    """Compute isoform ratios for a gene across cells."""
    gene_isoforms = adata.var[adata.var['gene'] == gene].index.tolist()
    if len(gene_isoforms) < 2:
        return None
    
    counts = adata[:, gene_isoforms].X
    if hasattr(counts, 'toarray'):
        counts = counts.toarray()
    
    total = counts.sum(axis=1)
    valid = total > 0
    ratios = np.zeros(len(total))
    ratios[valid] = counts[valid, 0] / total[valid]
    return ratios, valid

# 5. Test for cell-type-specific isoform usage
validation_results = []
for gene in svp_genes:
    result = compute_isoform_ratios(scalpel_adata, gene)
    if result is None:
        continue
    ratios, valid = result
    
    # Compare isoform ratios between cell types
    cell_types = scalpel_adata.obs['cell_type'].unique()
    for ct1, ct2 in combinations(cell_types, 2):
        mask1 = (scalpel_adata.obs['cell_type'] == ct1) & valid
        mask2 = (scalpel_adata.obs['cell_type'] == ct2) & valid
        
        if mask1.sum() < 10 or mask2.sum() < 10:
            continue
        
        stat, pval = mannwhitneyu(ratios[mask1], ratios[mask2], 
                                   alternative='two-sided')
        
        validation_results.append({
            'gene': gene,
            'cell_type_1': ct1,
            'cell_type_2': ct2,
            'mean_ratio_ct1': ratios[mask1].mean(),
            'mean_ratio_ct2': ratios[mask2].mean(),
            'pvalue': pval,
            'statistic': stat
        })

validation_df = pd.DataFrame(validation_results)

# 6. Multiple testing correction
validation_df['padj'] = multipletests(validation_df['pvalue'], method='fdr_bh')[1]

# 7. Assess concordance
significant_validations = validation_df[validation_df['padj'] < 0.05]
validated_genes = significant_validations['gene'].unique()
concordance_rate = len(validated_genes) / len(svp_genes) * 100
print(f"Validated {len(validated_genes)}/{len(svp_genes)} SVP genes ({concordance_rate:.1f}%)")

# 8. Save validation report
validation_df.to_csv('splisosm_scalpel_validation.csv', index=False)
```

### 5.4 Comparison: SPLISOSM vs SCALPEL

| Aspect | SPLISOSM | SCALPEL |
|--------|----------|---------|
| **Data type** | Spatial transcriptomics (Visium, Slide-seq, Xenium) | 3' scRNA-seq (10x Chromium, Drop-seq) |
| **Resolution** | Spot-level (55μm Visium) or near-single-cell | True single-cell |
| **Isoform coverage** | Full-length (ONT) or 3' end (Sierra TREND) | 3' end only (last 600nt) |
| **Statistical test** | HSIC-based kernel independence test | Chi-squared test for isoform distribution |
| **Spatial information** | ✅ Preserved and tested directly | ❌ Lost (requires spatial inference) |
| **Cell type resolution** | Requires deconvolution | ✅ Direct single-cell resolution |
| **RBP association** | ✅ Conditional HSIC with spatial confounding | Manual correlation analysis |
| **Complementary use** | Detects spatial patterns | Validates cell-type specificity |

---

## 6. Testable Biological Hypotheses

### 6.1 Hypothesis 1: VEGF Isoform Gradients Correlate with Tumor Hypoxia

Test using TREND (3' UTR length) changes as proxy for transcriptional regulation. Correlate with hypoxia markers (HIF1A, LDHA, CA9) expression. Expect shorter 3' UTRs in hypoxic core (rapid mRNA turnover).

### 6.2 Hypothesis 2: IL-6 Trans-Signaling Components at Leading Edge

Soluble IL-6R (sIL-6R) is partially generated by alternative splicing (~10%). Test for differential gp130 (IL6ST) APA at tumor-stroma interface. Correlate with ADAM17 expression (proteolytic sIL-6R source).

### 6.3 Hypothesis 3: TGF-β Isoform Ratios Distinguish Tumor-Educated Stroma

TGF-β1/β2/β3 have distinct receptor affinities (β2 requires betaglycan). Test spatial co-localization of TGF-β2 with endoglin (CD105) negative regions. Expect TGF-β2 enrichment in EMT zones.

### 6.4 Hypothesis 4: Chemokine Isoform Patterns Predict Immune Infiltration

CXCL12 isoforms (α vs β vs γ) have different ECM retention. Test whether CXCL12γ (high ECM affinity) localizes to stroma. Correlate with T cell and B cell infiltration patterns.

---

## 7. Interpretation Framework for Cytokine Isoforms

| Pattern | Biological Interpretation | Example |
|---------|---------------------------|---------|
| Agonist isoform in tumor core | Pro-tumor signaling | VEGF165 (angiogenic) |
| Antagonist isoform at margin | Tumor suppression at boundary | VEGF165b (anti-angiogenic) |
| Matrix-bound isoform localized | Gradient formation | CXCL12γ (high GAG binding) |
| Secreted isoform diffuse | Systemic signaling | CXCL12α (diffusible) |
| Intron-retaining isoform | Non-functional/NMD target | CD74 variants |

---

## 8. Validation Checklist

Before biological interpretation, validate findings with:

- [ ] **Cross-platform**: Do ONT and Sierra detect same genes?
- [ ] **Cross-sample**: Is pattern reproducible across replicates?
- [ ] **CLIP data**: Do associated RBPs bind near variable regions?
- [ ] **SCALPEL**: Do scRNA-seq show cell-type-specific usage?
- [ ] **Literature**: Are isoforms previously characterized?
- [ ] **Protein evidence**: Do isoforms produce distinct proteins?

---

## 9. Expected Outputs from Pan-Cancer Analysis

1. **Spatially variable cytokine gene catalog**: List of cytokines/secreted proteins with significant HSIC-IR (p < 0.05, FDR corrected) across multiple cancer types

2. **Tumor boundary-enriched isoform atlas**: Systematic identification of isoforms preferentially expressed at leading edge vs tumor core vs stroma

3. **Pan-cancer spatial isoform signatures**: Conserved vs cancer-type-specific patterns through hierarchical clustering of HSIC statistics

4. **Candidate therapeutic targets**: Isoform-specific cytokines with tumor-restricted spatial patterns amenable to therapeutic targeting (e.g., anti-angiogenic VEGF165b rescue strategies, IL-6 trans-signaling blockade at invasion zones)

---

## 10. Key References

### SPLISOSM
Su J, Qu Y, Schertzer M, et al. Mapping isoforms and regulatory mechanisms from spatial transcriptomics data with SPLISOSM. *Nature Biotechnology* 2026. doi: 10.1038/s41587-025-02965-6
- GitHub: https://github.com/JiayuSuPKU/SPLISOSM
- Documentation: https://splisosm.readthedocs.io/

### SCALPEL
Ake F, Schilling M, Fernández-Moya SM, et al. Quantification of transcript isoforms at the single-cell level using SCALPEL. *Nature Communications* 2025;16:6402. doi: 10.1038/s41467-025-61118-0
- GitHub: https://github.com/plasslab/SCALPEL
- Zenodo container: doi: 10.5281/zenodo.15717636

### Spatial Transcriptomics Validation
Fu Y, Kim H, Roy S, et al. Single cell and spatial alternative splicing analysis with Nanopore long read sequencing. *Nature Communications* 2025;16:6654. doi: 10.1038/s41467-025-60902-2

Lebrigand K, et al. The spatial landscape of gene expression isoforms in tissue sections. *Nucleic Acids Research* 2023;51:e47. doi: 10.1093/nar/gkac994

### Supporting Methods
Zhu J, Sun S, Zhou X. SPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns. *Genome Biology* 2021;22:184. doi: 10.1186/s13059-021-02404-0

---

*Document Version: 2.0 | Last Updated: February 2026*
*Covers: SPLISOSM adoption plan, platform constraints, target gene selection, SCALPEL validation workflow, biological hypotheses, interpretation framework*

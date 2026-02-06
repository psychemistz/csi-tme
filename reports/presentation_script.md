# CSI-TME Presentation Script

## Slide 1: Discovery — Spatial Isoform Variation in Glioma Secreted Proteins

**Title**: Can we detect alternative polyadenylation and splicing of secreted proteins in the glioma tumor microenvironment using spatial transcriptomics?

**Visual**: Table 1 (13 genes, HSIC-GC/IR/IC statistics across 24 samples)

### Script

We applied SPLISOSM — a spatial isoform statistical framework based on the Hilbert-Schmidt Independence Criterion — to detect spatially variable isoform patterns in glioma spatial transcriptomics data.

We tested secreted proteins and cytokines across 24 glioma samples: 13 Visium short-read and 11 Visium ONT long-read. The key test statistic is HSIC-IR, which measures whether isoform ratios vary spatially — indicating cell-type-dependent alternative polyadenylation or splicing.

The table shows our top 13 hits. HSIC-GC detects spatially variable gene expression — nearly universal for these genes. HSIC-IR, our primary statistic, detects spatially variable isoform ratios.

Seven genes reach TIER 1, meaning significant HSIC-IR in the majority of tested samples:
- SPP1 (osteopontin): significant in 14 of 20 samples, best p-value 5.1e-124 — detected on both SR and ONT platforms
- CLU (clusterin): 12/12 samples significant — ONT-dominant
- HMGB1 (alarmin): 11/16 — SR platform
- CD74: 8/10 — ONT
- APP: 4/10 — SR
- SPARCL1 and VCAN: ECM genes with consistent spatial isoform variation

Six additional genes reach TIER 2 with spatial-only evidence — notably SPARC (12/20 IR-significant, best p=2.2e-121) and B2M (11/11).

The SR platform primarily detects 3' UTR length variants from alternative polyadenylation, while ONT captures full-length isoforms including internal splicing events. Both platforms converge on the same genes, supporting genuine biological isoform switching.

---

## Slide 2: Validation — scRNA-seq Differential Splicing & RBP Regulators

**Title**: Independent validation by scRNA-seq and identification of candidate RNA-binding protein regulators

**Visual**: Table 2 (16 genes, scRNA-seq validation) + Table 3 (top RBP-gene pairs from conditional HSIC)

### Script

To validate the spatial isoform patterns, we performed differential splicing analysis on an independent scRNA-seq dataset — 32,000 SMART-seq2 cells from glioblastoma — using the MARVEL pipeline. For each gene, we tested all annotated splice events across pairwise cell-type comparisons. We defined significant events as those with adjusted p-value below 0.05 AND absolute delta-PSI above 10%.

The table is sorted by maximum delta-PSI. The top validated genes show dramatic isoform switching:
- VCAN: 99% delta-PSI between macrophages and tumor cells — near-complete isoform switching. All 6 tested splice events are significant.
- HMGB1: 96% delta-PSI — alternative last exon choice between tumor and immune cells
- CD74: 83% delta-PSI — immune-specific exon inclusion
- APP: 67% delta-PSI — the classic neuronal APP695 versus APP751 exon 7/8 skipping

These first four show the spatial isoform signal we detected with SPLISOSM is driven by genuine cell-type-specific splicing events.

SPP1 and CLU show moderate but significant delta-PSI at 28% and 19% respectively — real biological effects confirmed at single-cell level.

The TIER 2 genes with scRNA data mostly show sub-threshold effects — SERPINE2 has 18% max delta-PSI but doesn't pass significance, consistent with spatial-only evidence.

For the regulator analysis, we used conditional HSIC — a test that asks: does the isoform ratio correlate with an RNA-binding protein's expression after removing spatial confounding? We tested 13 target genes against 19 candidate RBPs across all 24 samples, totaling 971 individual tests.

One hit reaches FDR significance: HMGB1 and RBFOX1 at adjusted p-value 0.049. RBFOX1 is a neuronal splicing factor — this is a novel finding suggesting neuronal-compartment-specific regulation of HMGB1 alternative last exon choice.

The near-significant hit is SPP1 and ELAVL1/HuR at 0.062. This is the most literature-supported result: HuR binds AU-rich elements in SPP1's 3' UTR, and the longer 3' UTR variant has more HuR binding sites. The convergent signal from both ELAVL1 (ubiquitous) and ELAVL3 (neuronal paralog) on SPP1 strengthens the case for ELAV-family regulation of SPP1 APA.

---

## Slide 3: Spatial Validation — VCAN Isoform Switching (ZH1019T1)

**Title**: VCAN shows tumor-specific isoform 1 usage with near-complete spatial switching

**Visual**: 5-panel figure (vcan_validation_ZH1019T1.png)

### Script

This figure shows the spatial validation of VCAN isoform switching in sample ZH1019T1.

Panel A: VCAN is broadly expressed — 1,028 out of 1,109 spots, nearly ubiquitous across the tissue.

Panel B: SpaCET deconvolution shows the tumor (malignant) cell proportion. This sample has a clear tumor core with high malignant fraction, surrounded by regions with lower tumor content.

Panel C: Macrophage proportion from SpaCET — macrophages are present at low levels throughout, with some enrichment at the tissue periphery and in scattered foci.

Panel D: The VCAN isoform 1 fraction — this is the key panel. Red spots have high isoform 1 usage, blue spots have low. You can see the spatial pattern mirrors the tumor proportion in panel B: the tumor core is red, the periphery is blue.

Panel E: To quantify this, we computed Spearman correlations between each cell type's proportion and the VCAN isoform 1 fraction across all expressing spots. Malignant is the only cell type with a positive correlation — r equals +0.22, highly significant. Every immune and stromal cell type is negatively correlated. The strongest negative correlations are Endothelial at -0.21, CAF at -0.20, and pDC at -0.19.

VCAN isoform 1 maps to the last exon at chr5:83,581,043 — the canonical 3' UTR polyadenylation site. Tumor cells preferentially use this proximal APA site, while immune and stromal cells use alternative upstream peaks mapping to exon 8, the GAG-beta domain. This is consistent with the scRNA-seq validation showing 99% delta-PSI — near-complete isoform switching between tumor and non-tumor cell types.

---

## Slide 4: Spatial Validation — HMGB1 Isoform Switching (ZH8811Abulk)

**Title**: HMGB1 alternative last exon usage is spatially resolved in the largest glioma sample

**Visual**: 5-panel figure (hmgb1_validation_ZH8811Abulk.png)

### Script

This is the same analysis for HMGB1, shown in sample ZH8811Abulk — the sample with the strongest HSIC-IR significance for HMGB1 at adjusted p-value 4.2e-80, and also the sample that drove the novel RBFOX1 conditional HSIC hit.

Panel A: HMGB1 is expressed in 2,278 of 2,383 spots — this is the largest sample with the most cellular diversity.

Panel B: Tumor proportion shows a complex spatial architecture with distinct tumor regions.

Panel C: Macrophage proportion — higher macrophage enrichment compared to ZH1019T1, with clear spatial segregation from tumor-dominant regions.

Panel D: HMGB1 isoform 1 fraction shows a spatial gradient that tracks with tumor proportion. Red regions correspond to tumor-rich areas, blue to immune-infiltrated regions. The pattern is more heterogeneous than VCAN, consistent with HMGB1's mechanism being alternative last exon choice rather than simple APA.

Panel E: The correlation bar chart shows the same pattern as VCAN but with stronger effect sizes in this sample. Malignant is the only positive cell type at r=+0.20, highly significant. The negative correlations are broad and strong: Plasma r=-0.17, CAF r=-0.16, Endothelial r=-0.16, Macrophage r=-0.15. Every immune and stromal cell type is negative, all significant at p<0.001 except T CD8 and cDC.

HMGB1 is an alarmin — a damage-associated molecular pattern protein. Its alternative last exon choice produces different C-terminal domains affecting protein secretion and immune signaling. The tumor compartment uses one last exon while immune cells use another, with 96% delta-PSI in scRNA-seq validation. This is the gene where we identified RBFOX1 as a candidate regulator via conditional HSIC — suggesting the neuronal splicing factor RBFOX1 influences HMGB1 3' end processing in the tumor-neuron interface of the glioma microenvironment.

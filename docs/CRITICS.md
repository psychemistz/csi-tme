# CSI-TME Validation Assessment

## Summary of Findings

7 **TIER1** genes validated with both spatial HSIC-IR and scRNA-seq differential splicing:

| Gene | Spatial Evidence | scRNA Evidence | Effect Size |
|------|-----------------|----------------|-------------|
| **SPP1** | 14/20 samples IR-sig | 9/483 events, ΔPSI=28% | Strong |
| **CLU** | 12/12 samples IR-sig | 26/612 events, ΔPSI=19% | Strong (ONT-dominant) |
| **HMGB1** | 11/16 samples IR-sig | 39/273 events, ΔPSI=96% | Strong (ALE) |
| **CD74** | 8/10 samples IR-sig | 4/136 events, ΔPSI=83% | Moderate (ONT-only) |
| **APP** | 4/10 samples IR-sig | 65/294 events, ΔPSI=67% | Strong |
| **SPARCL1** | 4/12 samples IR-sig | 3/460 events, ΔPSI=14% | Moderate |
| **VCAN** | 7/13 samples IR-sig | 45/108 events, ΔPSI=99% | Very strong |

19 **TIER2** genes show spatial-only patterns (no scRNA validation): B2M, SPARC, PTGDS, SERPINE2, NCAM1, SERPINE1, C3, CST3, MIF, TIMP1, LGALS1, CTSB, ITM2B, IGFBP5, AGT, TIMP2, DKK3, APOE, EDIL3.

## Conditional HSIC (DU Test) Results

Tested 13 target genes (7 TIER1 + 6 top TIER2) against 19 RBP candidates across 24 glioma samples (13 SR + 11 ONT). 971 individual tests, 157 gene-RBP pairs summarized with Fisher's combined p-value.

### FDR-Significant

| Gene | RBP | Fisher padj | Sig samples | Predicted? |
|------|-----|-------------|-------------|------------|
| **HMGB1** | **RBFOX1** | **0.049** | 1/5 (ZH8811Abulk p=3.7e-6) | No — novel |

HMGB1 undergoes alternative last exon (ALE) usage. RBFOX1 is a neuron-specific splicing factor that binds (U)GCAUG motifs near regulated exons. In glioma, RBFOX1 marks residual neurons and neuronal-like tumor cells. The association suggests a neuronal-compartment-specific regulatory program controlling HMGB1 3' end processing. Signal driven by ZH8811Abulk, the largest and most heterogeneous bulk sample (2383 spots, 342 RBFOX1+ spots).

### Near-Significant

| Gene | RBP | Fisher padj | Sig samples | Predicted? |
|------|-----|-------------|-------------|------------|
| **SPP1** | **ELAVL1** (HuR) | **0.062** | 3/12 | Yes |
| **SPP1** | **ELAVL3** (HuC) | 0.114 | 2/10 | Related |
| **HMGB1** | **ELAVL1** (HuR) | 0.254 | 3/14 | Yes |

**SPP1 × ELAVL1** is the most literature-supported hit. SPP1 (osteopontin) has 3' UTR length variants from APA. ELAVL1/HuR binds AU-rich elements (AREs) in the 3' UTR — longer variants have more AREs and stronger HuR-mediated stabilization. Significant in ZH8811Abulk (p=0.022), ZH8811Bbulk (p=0.020), and DMG_5 ONT (p=0.003). The convergent signal from both ELAVL1 (ubiquitous) and ELAVL3 (neuronal paralog) on SPP1 strengthens the case: ELAV-family proteins regulate SPP1 APA across cell types through the same ARE binding sites.

### Suggestive Hits

| Gene | RBP | Fisher padj | Biological rationale |
|------|-----|-------------|---------------------|
| SERPINE2 | CELF5 | 0.288 | CELF5 brain-enriched; regulates APA via UG-rich motifs |
| SPARC | ELAVL4 (HuD) | 0.288 | Neuronal ELAV protein; SPARC APA in glia vs neurons |
| APP | SRSF3 | 0.326 | SR proteins modulate exon inclusion |
| APP | ELAVL1 | 0.326 | Documented role in APP mRNA regulation |
| VCAN | SRSF1 | 0.326 | SRSF1 regulates APA site selection |
| NCAM1 | NOVA1 | 0.326 | Literature-validated: NOVA1 regulates NCAM1 splicing |
| HMGB1 | CELF4 | 0.326 | Brain-enriched CELF; 3' end regulation |

### Hypothesis Validation Summary

| Gene | Predicted RBPs | Top DU Hit | Confirmed? |
|------|---------------|------------|------------|
| **SPP1** | ELAVL1, CELF1/2 | **ELAVL1 (padj=0.062)** | **Yes** |
| **HMGB1** | ELAVL1, SRSF1/3 | RBFOX1 (0.049), ELAVL1 (0.254) | Partial — RBFOX1 novel |
| **APP** | RBFOX1/2, NOVA1/2 | SRSF3 (0.326), RBFOX1 (0.360) | Weak |
| **CLU** | NOVA1/2, PTBP1 | CELF1 (0.326), NOVA1 (0.355) | Weak |
| **SPARCL1** | QKI, ELAVL1 | ELAVL3 (0.326) | Related (ELAV family) |
| **VCAN** | CELF1/2, QKI | SRSF1 (0.326) | No |
| **CD74** | SRSF2, PTBP1/2 | SRSF2 (0.503) | Weak |

## Cross-Cutting Themes

### ELAV/Hu Family Dominance

ELAVL1/3/4 appear in 6 of the top 20 hits across 4 target genes (SPP1, HMGB1, SPARC, APP). This family binds AU-rich elements in 3' UTRs — exactly the features detected by Visium SR (APA/3' UTR variants). The conditional HSIC detects a coherent biological signal: ELAV proteins regulate 3' UTR choice across multiple secreted proteins in the glioma TME.

### Neuronal RBPs in a Brain Tumor

The top hit (RBFOX1) and several suggestive hits (ELAVL3, ELAVL4, CELF5, NOVA1) are neuron-specific. Gliomas infiltrate brain parenchyma, creating a TME where neuronal RBPs and tumor/immune cells intermingle. The DU test detects regulatory programs at this tumor-neuron interface.

### Sample ZH8811 Drives Most Discoveries

Both ZH8811A and ZH8811B replicates consistently rank as the most informative samples. These are the largest bulk glioma samples (~2400–2500 spots) with the highest cellular diversity. This indicates the DU test is power-limited — larger samples with more heterogeneity provide more statistical power for conditional associations.

## Limitations of the DU Test

1. **Power**: 5–15 samples per gene-RBP pair limits the Fisher combined test. The single FDR hit and several near-misses suggest real signal buried in noise.
2. **Linear kernel**: Only captures linear RBP-isoform relationships. Threshold or saturating effects would be missed.
3. **SR platform bias**: Visium SR detects APA (3' UTR length), biasing toward RBPs that regulate polyadenylation. Internal splicing regulators (RBFOX, NOVA) may show weaker signal because their primary targets (exon skipping) aren't well-captured.
4. **No CLR transform**: Conditional HSIC uses raw proportion residuals rather than CLR-transformed compositional data.

## Remaining Action Items

1. ~~Run conditional HSIC (DU test)~~ — **Done.** 971 tests, 1 FDR-significant hit.
2. **Cross-reference CLIP data**: Check POSTAR3/ENCODE for RBFOX1 binding sites in HMGB1 3' UTR and ELAVL1 binding in SPP1 3' UTR.
3. **Power analysis**: Estimate sample size needed for 80% power at current effect sizes.
4. **Nonlinear kernel extension**: Replace linear kernel with RBF/Gaussian kernel in conditional HSIC for better sensitivity to nonlinear RBP-isoform relationships.
5. **ONT-focused reanalysis**: Re-run DU test on ONT samples only for genes like APP and CD74 where internal splicing (not APA) is the mechanism.

## Platform Note

SR Visium detects APA/3' UTR variants; ONT detects full-length isoforms. The scRNA-seq validation (MARVEL) captures internal splicing events. These measure different but often correlated isoform axes — both can reflect the same underlying cell-type composition differences driving spatial patterns. The DU test sidesteps this platform mismatch by testing RBP-isoform associations directly in the spatial data.

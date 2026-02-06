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

## Key Gap: No Regulator Identification

The SPLISOSM paper (Fig 3-5) uses **conditional HSIC** to identify RBP regulators driving isoform patterns. CSI-TME has only run unconditional SV tests. Without regulator identification:

- We know *where* isoform switching occurs but not *why*
- We cannot distinguish cell-intrinsic splicing programs from microenvironment-driven regulation
- We cannot prioritize mechanistic follow-up experiments

This is the single most important gap to address.

## RBP-Isoform Hypotheses per TIER1 Gene

| Gene | APA/Isoform Mechanism | Candidate RBPs | Rationale |
|------|----------------------|-----------------|-----------|
| **SPP1** | 3' UTR length (APA) | ELAVL1 (HuR), CELF1/2 | HuR stabilizes AU-rich cytokine mRNAs; CELF binds UG-rich 3' UTR motifs |
| **CLU** | Alternative first exon (AFE) | NOVA1/2, PTBP1 | NOVA regulates neuronal exon selection; PTBP1 represses in non-neuronal cells |
| **HMGB1** | Alternative last exon (ALE) | ELAVL1, SRSF1/3 | ELAVL1 regulates 3' end processing; SR proteins modulate ALE choice |
| **CD74** | Intron retention (ONT) | SRSF2, PTBP1/2 | SRSF2 promotes intron inclusion; PTB regulates in immune cells |
| **APP** | Exon 7/8 skipping | RBFOX1/2, NOVA1/2 | RBFOX/NOVA regulate neuronal APP695 vs APP751 |
| **SPARCL1** | 3' UTR variants | QKI, ELAVL1 | QKI regulates ECM gene mRNA processing in glia |
| **VCAN** | 3' UTR length (APA) | CELF1/2, QKI | CELF family regulates ECM gene APA; QKI expressed in glia |

## Action Items

1. **Run conditional HSIC (DU test)**: Test 19 RBPs × 7 TIER1 genes + top TIER2 genes → `scripts/run_du_test.py`
2. **Interpret significant hits**: RBP-gene pairs with padj < 0.05 after spatial residualization indicate non-spatial regulatory associations
3. **Cross-reference CLIP data**: For significant hits, check POSTAR3/ENCODE for RBP binding evidence in target gene 3' UTRs
4. **Expand to TIER2**: Genes like SPARC, B2M, SERPINE2 with strong spatial patterns may reveal additional RBP associations

## Platform Note

SR Visium detects APA/3' UTR variants; ONT detects full-length isoforms. The scRNA-seq validation (MARVEL) captures internal splicing events. These measure different but often correlated isoform axes — both can reflect the same underlying cell-type composition differences driving spatial patterns. The DU test sidesteps this platform mismatch by testing RBP-isoform associations directly in the spatial data.

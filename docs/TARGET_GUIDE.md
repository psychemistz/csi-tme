# Target Selection Guidelines for Cytokine and Secreted Protein Isoform Analysis in Visium 3′ Short-Read Spatial Transcriptomics

## Executive Summary

This document provides systematic criteria for selecting cytokines and secreted proteins from CytoSig and SecAct gene lists for spatial isoform analysis using 10x Genomics Visium 3′ short-read (SR) chemistry. The fundamental constraint is that **Visium 3′ SR data only reliably detect isoform differences arising from 3′ architecture variations** (alternative polyadenylation, terminal exon usage, and unique 3′ exonic sequences). Internal splicing events, alternative promoters, and 5′ differences are not detectable.

---

## Part I: Technical Constraints of Visium 3′ Short-Read Chemistry

### What Visium 3′ SR Can Detect

| Isoform Feature | Detectability | Biological Relevance |
|----------------|---------------|---------------------|
| Alternative polyadenylation (APA) | **High** | 3′ UTR length affects mRNA stability, localization, translation efficiency |
| 3′ UTR length variation | **High** | miRNA binding sites, RBP regulatory elements |
| Terminal exon choice | **Moderate** | C-terminal protein variation, functional antagonism |
| Unique 3′ exonic sequences | **Low-Moderate** | Context-dependent, requires sufficient sequence uniqueness |

### What Visium 3′ SR Cannot Detect

| Isoform Feature | Reason for Limitation |
|----------------|----------------------|
| Internal exon skipping | Requires junction-spanning reads across full transcript |
| Alternative splice donor/acceptor sites | Internal to transcript, not captured |
| Intron retention | Poor internal coverage |
| Alternative TSS / 5′ UTR variation | 5′ regions not captured by poly(A) priming |
| Full isoform reconstruction | Short reads + fragmentation eliminate connectivity |

### Interpretation Framework

Any "isoform" signal from Visium 3′ SR data represents:
- Indirect inference from 3′-end read distributions
- Annotation-dependent quantification
- Relative enrichment at annotated 3′ features, not direct isoform measurement

---

## Part II: CytoSig Gene List Overview

CytoSig (Cytokine Signaling Analyzer) covers **43 cytokines, chemokines, and growth factors** with curated target gene signatures. The database includes:

### CytoSig Categories

| Category | Example Members | Count |
|----------|-----------------|-------|
| Interleukins | IL1B, IL2, IL4, IL6, IL10, IL12, IL13, IL15, IL17, IL18, IL21, IL22, IL23, IL33 | ~20 |
| Interferons | IFNA, IFNB, IFNG | 3 |
| TGF-β family | TGFB1, TGFB2, TGFB3 | 3 |
| TNF family | TNF, TRAIL (TNFSF10), LTA | ~5 |
| Growth factors | VEGFA, EGF, HGF, PDGF, FGF | ~8 |
| Chemokines | CXCL12, CCL2, CCL5 | ~5 |

**Note:** Not all CytoSig genes have documented isoforms with 3′ differences detectable by Visium.

---

## Part III: Selection Criteria for Spatial Isoform Analysis

### Tier 1: Priority Targets (3′ Differences Documented)

These genes have well-characterized isoform differences localized to the 3′ region and are amenable to Visium 3′ SR analysis:

#### VEGFA (Vascular Endothelial Growth Factor A)
- **Isoform Mechanism:** Alternative 3′ splice site selection in exon 8 (proximal vs distal)
- **Functional Consequence:** Pro-angiogenic (VEGFxxx) vs anti-angiogenic (VEGFxxxb) isoforms
- **Disease Relevance:** Cancer, diabetic retinopathy, macular degeneration
- **Detectability:** **HIGH** - 3′ terminal exon difference
- **Key Reference:** Bates DO et al. Cancer Res 2002;62:4123-4131

#### IL33 (Interleukin-33)
- **Isoform Mechanism:** Alternative splicing of exons 3, 4, 5 affecting C-terminal activity domain
- **Functional Consequence:** Δexon 3,4 isoform is cytoplasmic (more readily secreted) vs nuclear full-length
- **Disease Relevance:** Asthma, type 2 inflammation
- **Detectability:** **MODERATE-HIGH** - Terminal exon differences affect 3′ read distribution
- **Key Reference:** Gordon ED et al. PNAS 2016;113(31):8765-8770

#### CXCL12 (SDF-1)
- **Isoform Mechanism:** Alternative splicing produces α, β, γ isoforms with different C-termini
- **Functional Consequence:** CXCL12γ has extended C-terminus with GAG-binding domain (retention vs diffusion)
- **Disease Relevance:** Cancer metastasis, stem cell homing
- **Detectability:** **MODERATE** - C-terminal differences visible at 3′ end
- **Key Reference:** Dimberg J et al. Int J Mol Med 2016;37:1197-1203

#### TNFSF10 (TRAIL)
- **Isoform Mechanism:** Alternative splicing produces TRAIL-short (antagonist)
- **Functional Consequence:** TRAIL-short lacks death-inducing domain, antagonizes apoptosis
- **Disease Relevance:** HIV pathogenesis, cancer immune evasion
- **Detectability:** **MODERATE** - Terminal differences
- **Key Reference:** Schnepple DJ et al. J Clin Invest 2011;121:4446-4461

### Tier 2: Moderate Priority Targets (Conditional Detectability)

These genes have isoforms but detection depends on specific conditions:

| Gene | Isoform Type | Limitation for Visium 3′ SR |
|------|-------------|---------------------------|
| **IL15** | Long vs short signal peptide | 5′ difference - **NOT detectable** |
| **IL4** | IL-4δ2 antagonist (exon 2 skip) | Internal splicing - **NOT detectable** |
| **TGFB1/2/3** | Isoforms are separate genes | Compare gene-level expression, not splice variants |
| **IL6** | No major 3′ isoform documented | APA may exist but not well characterized |
| **IFNG** | Single exon structure | Limited isoform diversity |

### Tier 3: Low Priority / Not Recommended

These CytoSig genes are **not suitable** for Visium 3′ SR isoform analysis:

| Gene | Reason |
|------|--------|
| IL1B | Functional regulation by proteolysis, not splicing |
| IL18 | Caspase activation, no major 3′ splice variants |
| TNF | Membrane-bound vs soluble by proteolysis |
| Most chemokines | Short transcripts, minimal 3′ diversity |

---

## Part IV: Decision Framework for Target Selection

### Step-by-Step Selection Process

```
Step 1: Identify genes of interest from CytoSig/SecAct
           ↓
Step 2: Query transcript annotation databases (Ensembl, GENCODE, RefSeq)
        → How many annotated transcripts exist?
        → Do transcripts differ at the 3′ end?
           ↓
Step 3: Literature review for functional isoform differences
        → Is there documented biological consequence?
        → Is isoform ratio disease-associated?
           ↓
Step 4: Assess 3′ architecture
        → APA sites documented? (PolyASite, APADB)
        → Terminal exon alternatives?
        → Unique 3′ sequence elements?
           ↓
Step 5: Platform suitability check
        → 3′ difference: Visium 3′ SR suitable
        → 5′ or internal difference: Requires long-read or probe-based platform
           ↓
Step 6: Validate with existing spatial/scRNA-seq data (if available)
```

### Annotation Database Resources

| Resource | URL | Use Case |
|----------|-----|----------|
| Ensembl | ensembl.org | Transcript isoform annotation |
| GENCODE | gencodegenes.org | Comprehensive transcript catalog |
| PolyASite | polyasite.unibas.ch | APA site database |
| APAdb | 3utr.org | 3′ UTR and APA annotation |
| UniProt | uniprot.org | Protein isoform functional annotation |
| IsoExpresso | github.com/gao-lab/IsoExpresso | Cancer-specific isoform expression |

---

## Part V: APA as Primary Signal for Cytokine Spatial Isoforms

### Why APA Matters for Cytokines

Alternative polyadenylation is **the most reliable isoform signal** from Visium 3′ SR data. For cytokines and secreted proteins, APA can regulate:

1. **mRNA Stability:** Shorter 3′ UTRs often increase stability by evading miRNA regulation
2. **Translation Efficiency:** AU-rich elements in extended 3′ UTRs can destabilize transcripts
3. **Localization:** 3′ UTR elements can direct subcellular localization of mRNA
4. **Protein Output:** Net effect on protein abundance in the tumor microenvironment

### APA in Immune Response

| Context | APA Pattern | Functional Consequence |
|---------|-------------|----------------------|
| T cell activation | Global 3′ UTR shortening | Increased protein output of effector molecules |
| Macrophage activation | Shortening of cytokine transcripts | Enhanced inflammatory response |
| Viral infection (VSV) | Progressive 3′ UTR shortening | JAK-STAT, TLR pathway genes affected |
| Tumor microenvironment | Variable by cell type | Context-dependent immune modulation |

### APA-Regulated Cytokine Candidates

Based on literature evidence for APA regulation during immune responses:

| Gene | Evidence for APA | Spatial Biology Hypothesis |
|------|-----------------|---------------------------|
| IL6 | 3′ UTR lengthening with IL-6 | Spatial variation in IL-6 stability across TME |
| CXCL8 (IL-8) | APA during inflammation | Neutrophil chemotaxis gradient formation |
| IFNG | Limited APA documentation | Less suitable for APA analysis |
| TNF | Documented APA in inflammation | May show spatial 3′ UTR variation |
| IL10 | APA in macrophages | Immunosuppressive niche detection |

---

## Part VI: Practical Implementation

### Recommended Target List for SPLISOSM Analysis

#### High-Confidence Targets (Start Here)

| Gene | Signal Type | Biological Question |
|------|-------------|---------------------|
| VEGFA | Terminal exon (8a vs 8b) | Pro- vs anti-angiogenic balance in tumor |
| IL33 | Exon 3,4,5 alternatives | Alarmin release mechanism in inflammation |
| CXCL12 | C-terminal isoforms | Chemokine retention vs diffusion gradients |

#### APA-Based Targets (Novel Discovery)

| Gene | Approach | Expected Insight |
|------|----------|-----------------|
| IL6 | APA ratio (proximal/distal) | Stability regulation across spatial domains |
| TGFB1 | APA sites if documented | CAF-tumor communication |
| CCL2 | 3′ UTR variation | Monocyte recruitment gradients |

### Validation Considerations

1. **Cross-Platform Validation**
   - Use SCALPEL or similar scRNA-seq isoform tools on matched samples
   - Compare 3′ peaks with bulk RNA-seq from similar tissues
   
2. **Orthogonal Evidence**
   - CLIP-seq data for RBP binding (RBFOX, CELF families)
   - Published APA profiles from PolyA-seq studies
   
3. **Biological Validation**
   - Protein-level confirmation (ELISA, IHC for isoform-specific epitopes)
   - Functional assays for isoform activity

---

## Part VII: Exclusion Criteria

### Genes to Avoid for Visium 3′ SR Isoform Analysis

| Exclusion Criterion | Example Genes | Reason |
|--------------------|---------------|--------|
| 5′-only isoform differences | IL15 (signal peptide) | Not captured by 3′ chemistry |
| Internal exon skipping only | IL4δ2 (exon 2) | No junction coverage |
| Regulation by proteolysis | IL1B, IL18, TNF | Isoform control is post-translational |
| Single-transcript genes | Many chemokines | No isoform diversity to detect |
| Poor annotation quality | Novel genes | Cannot distinguish signal from noise |
| Very short transcripts | Small cytokines | Insufficient 3′ window for analysis |

---

## Part VIII: Summary Decision Matrix

### Quick Reference Table for CytoSig Genes

| Gene | Has 3′ Isoforms | Visium 3′ Suitable | Priority |
|------|-----------------|-------------------|----------|
| VEGFA | ✓ (terminal exon) | ✓✓✓ | **HIGH** |
| IL33 | ✓ (terminal exons) | ✓✓ | **HIGH** |
| CXCL12 | ✓ (C-terminal) | ✓✓ | **HIGH** |
| TNFSF10 | ✓ (TRAIL-short) | ✓✓ | **MODERATE** |
| IL6 | APA only | ✓ | **MODERATE** |
| TGFB1/2/3 | Separate genes | Gene-level only | **LOW** |
| IL15 | ✓ but 5′ | ✗ | **NOT SUITABLE** |
| IL4 | ✓ but internal | ✗ | **NOT SUITABLE** |
| IL1B | No splice isoforms | ✗ | **NOT SUITABLE** |

---

## Part IX: References

### Key Literature

1. **Visium 3′ Chemistry Limitations:**
   - 10x Genomics. Visium Spatial Gene Expression. Technical documentation.

2. **CytoSig Database:**
   - Jiang P et al. Systematic investigation of cytokine signaling activity at the tissue and single-cell levels. *Nature Methods* 2021;18:1181-1191. doi: 10.1038/s41592-021-01274-5

3. **Cytokine Isoforms:**
   - Gordon ED et al. Alternative splicing of interleukin-33 and type 2 inflammation in asthma. *PNAS* 2016;113:8765-8770. doi: 10.1073/pnas.1601914113
   - Bates DO et al. VEGF165b, an inhibitory splice variant of vascular endothelial growth factor. *Cancer Res* 2002;62:4123-4131.
   - Dimberg J et al. CXCL12 and CXCR4 in colorectal cancer. *Int J Mol Med* 2016;37:1197-1203.

4. **APA in Immune Response:**
   - Jia X et al. The role of alternative polyadenylation in the antiviral innate immune response. *Nature Communications* 2017;8:14605. doi: 10.1038/ncomms14605
   - Gruber AR et al. Alternative cleavage and polyadenylation in health and disease. *Nature Reviews Genetics* 2019;20:599-614.

5. **Spatial Isoform Methods:**
   - Su J et al. SPLISOSM: Spatial Isoform Analysis in Single-Cell and Spatial Transcriptomics. *Nature Biotechnology* 2026. doi: 10.1038/s41587-025-02965-6
   - Ake F et al. Quantification of transcript isoforms at the single-cell level using SCALPEL. *Nature Communications* 2025;16:6402. doi: 10.1038/s41467-025-61118-0

---

## Appendix A: Detailed Gene Cards for Priority Targets

### VEGFA

**Gene:** Vascular Endothelial Growth Factor A  
**Chromosome:** 6p21.1  
**Exons:** 8  
**Key Isoforms:**
- VEGF165a (pro-angiogenic): PSS selection in exon 8
- VEGF165b (anti-angiogenic): DSS selection in exon 8, 6 aa C-terminal difference

**3′ Architecture:**
- Alternative 3′ splice site in exon 8 creates two families
- Both families share identical 5′ structure
- 3′ UTR also has multiple APA sites

**Spatial Hypothesis:** Tumor core may show VEGF165a dominance; normal stroma may retain VEGF165b. Ratio may predict bevacizumab response.

---

### IL33

**Gene:** Interleukin-33  
**Chromosome:** 9p24.1  
**Exons:** 7  
**Key Isoforms:**
- Full-length (NM_033439): Nuclear localization
- Δexon 3,4 variant: Cytoplasmic, readily secreted
- Δexon 5: Inactive (lacks ST2 binding domain)

**3′ Architecture:**
- Exon 5 is critical for receptor binding (cytokine activity domain)
- Alternative terminal structures affect 3′ read patterns

**Spatial Hypothesis:** Type 2-high asthma epithelium expresses Δexon 3,4 correlating with inflammation. Spatial mapping could identify IL33-active microdomains.

---

### CXCL12

**Gene:** C-X-C Motif Chemokine Ligand 12 (SDF-1)  
**Chromosome:** 10q11.21  
**Exons:** 5  
**Key Isoforms:**
- CXCL12α: Most common, diffusible
- CXCL12β: 4 aa C-terminal extension
- CXCL12γ: 30 aa C-terminal extension, strong GAG binding (retained)

**3′ Architecture:**
- C-terminal differences encoded by alternative terminal exon usage
- γ isoform creates heparan sulfate-binding domain for matrix retention

**Spatial Hypothesis:** CXCL12γ may create immobilized chemokine gradients in stroma; CXCL12α may form diffusible signals. Spatial ratio could reveal chemotaxis guidance mechanisms.

---

## Appendix B: Computational Pipeline Considerations

### For SPLISOSM Analysis

1. **Input Requirements:**
   - Visium count matrix
   - Spatial coordinates
   - Isoform annotation (3′ features only)

2. **Preprocessing:**
   - Focus analysis on genes with documented 3′ isoform differences
   - Aggregate spots to pseudo-bulk for stable isoform ratios
   - Apply ICAR spatial smoothing kernel

3. **Statistical Framework:**
   - HSIC for spatial isoform correlation
   - Permutation testing for significance
   - Multiple testing correction across genes

4. **Interpretation:**
   - Positive spatial correlation: Isoform co-regulation in space
   - Differential isoform ratios: Microenvironment-specific splicing
   - Integration with cell type deconvolution for source attribution

---

*Document Version: 1.0*  
*Last Updated: February 2026*  
*For use with SPLISOSM spatial isoform analysis framework*

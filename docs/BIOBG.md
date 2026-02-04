# Biological Rationale: Alternative Splicing as a Regulatory Mechanism for Cytokines and Secreted Proteins

## Executive Summary

This document establishes the biological foundation for analyzing cytokine and secreted protein isoforms in spatial transcriptomics. While these proteins are traditionally understood to be regulated by proteolytic activation, emerging evidence demonstrates that alternative splicing provides a complementary regulatory layer that can gate secretion competence, modulate proteolytic maturation, alter receptor binding, generate antagonistic molecules, and reshape protease networks.

---

## 1. Motivation: Why Study Isoform-Level Spatial Patterns in Cytokines?

Cytokines and secreted proteins are traditionally understood to be regulated by **proteolytic activation**—many are synthesized as inactive precursors (preproproteins, zymogens, pro-cytokines) and converted to active forms by cleavage events. However, emerging evidence demonstrates that **alternative splicing provides a complementary regulatory layer** that can:

1. **Gate secretion competence** (presence/absence of signal peptides)
2. **Modulate proteolytic maturation** (include/exclude inhibitory propeptides or cleavage sites)
3. **Alter receptor binding and extracellular retention** (C-terminal variations affecting glycosaminoglycan binding)
4. **Generate antagonistic/decoy molecules** (truncated isoforms that bind but do not signal)
5. **Reshape protease networks** (splicing in proteases themselves)

The most compelling evidence shows:
- The existence of **distinct splice isoforms** with demonstrably different secretion/localization, receptor binding, or bioactivity
- **Condition-, tissue-, or disease-associated shifts** in isoform ratios that plausibly modulate physiological effects
- Mechanistic links between splicing and activation, such as splicing that deletes/includes inhibitory propeptides, alters cleavage motifs, or reshapes receptor-binding domains

---

## 2. Canonical Proteolytic Activation Mechanisms

### 2.1 Preproprotein Processing in the Secretory Pathway

For proteins entering the classical secretory pathway, an N-terminal **signal peptide** targets the nascent chain to the endoplasmic reticulum (ER), after which **signal peptidase** cleaves the signal peptide to release the mature region into the ER lumen for folding and trafficking. Signal peptide cleavage is often necessary for maturation and correct localization.

Many secreted proteins then undergo additional processing at specific sites by **proprotein convertases** (PCs; e.g., furin and related enzymes). PCs themselves are synthesized as zymogens and typically auto-activate via regulated cleavage steps; they cleave diverse substrates, including growth factors, receptors, and enzymes, often in the trans-Golgi network or secretory vesicles.

### 2.2 Zymogen Activation for Secreted Enzymes

Many extracellular proteases are synthesized as **inactive zymogens** containing inhibitory pro-domains. Proteolytic activation removes or disrupts these inhibitory segments, enabling catalytic activity; this provides fast, spatially restricted control and helps prevent unwanted proteolysis during biosynthesis and transport.

For **matrix metalloproteinases (MMPs)**, a canonical regulatory step is conversion from proMMP to active MMP by cleavage of the propeptide and disruption of the "cysteine switch," which exposes the catalytic zinc and enables proteolysis. Diverse proteases (including other MMPs, serine proteases, and sometimes furin-like convertases) can participate in activation depending on the MMP and context.

### 2.3 Proteolytic Activation of Select Cytokines

Some cytokines are synthesized as **inactive precursors** and require cleavage for full activity or altered localization. A prominent case is the **IL-1 family**, where the "signal 1 / signal 2" logic of inflammasome biology illustrates proteolytic control: priming induces transcription/accumulation of pro-cytokines, and inflammasome activation leads to **caspase-1** activation, which processes precursors such as **pro–IL‑1β** and **pro–IL‑18** into mature cytokines.

The biological reality can be more nuanced than a single mandatory cleavage step. For example, IL-33 is widely discussed as an alarmin with complex regulation by proteolysis (including activation or inactivation depending on the protease context) and by non-classical secretion mechanisms.

### 2.4 Conceptual Schematic of Proteolytic Control

```
CANONICAL PROTEOLYTIC CONTROL:
Gene → mRNA → Translation → Precursor protein (prepro/pro/membrane-tethered)
                  ↓
         Signal peptide cleavage (secretory entry)
                  ↓
         Proprotein / zymogen / pro-cytokine
                  ↓
         Activating cleavage by convertase/caspase/protease
                  ↓
         Mature bioactive secreted protein
```

This "precursor → cleavage → active product" architecture is a recurring design in extracellular biology because it supports rapid activation, compartmentalization, and strong on/off behavior.

---

## 3. Alternative Splicing as a Regulatory Mechanism

### 3.1 How Alternative Splicing Can Regulate Extracellular Function

**Alternative splicing** controls which exons (or splice sites) are used to build mature mRNA. It can therefore change the coding sequence of secreted proteins in ways that affect:

- **Secretion competence** (e.g., gain/loss of a signal peptide; altered trafficking motifs)
- **Proteolytic maturation potential** (e.g., deletion of inhibitory propeptides or removal of a cleavage site)
- **Receptor binding and extracellular retention** (e.g., altered C-termini affecting co-receptor binding, glycosaminoglycan binding, or dimerization)
- **Production of antagonistic/decoy molecules** (e.g., truncated or altered isoforms that bind but do not signal)
- **Protease-network behavior** (e.g., splicing in the protease itself that changes activity, secretion, or localization)

Because splicing acts *before translation*, it is typically slower than post-translational cleavage and often more "programmed" by cell type, activation state, and splicing-factor availability.

### 3.2 How Splicing Complements Proteolytic Activation

```
SPLICING-BASED CONTROL:
Gene → pre-mRNA → [SPLICING DECISION] → Isoform A or Isoform B
                                              ↓
                            Different precursor structures:
                            - With/without signal peptide
                            - With/without propeptide
                            - Different C-terminal domains
                            - Agonist vs. antagonist function
```

**Key Insight**: Splicing acts *before* translation, implementing a transcription-coupled choice among protein architectures. This is slower than post-translational cleavage but more "programmable" by cell type, activation state, and splicing factor availability.

A compact way to see complementarity:
```
Splicing control:     pre-mRNA → isoform choice → (which precursor exists?)
Proteolytic control:  precursor protein → cleavage → (which mature form exists now?)
```

This means splicing can **gate** whether a protein is even eligible for canonical activation (e.g., by including/excluding a signal peptide or a key domain), while proteolysis can **gate** the timing and location of activity for whichever isoform is present.

---

## 4. Well-Documented Cytokine Examples with Functional Isoform Differences

| Cytokine | Splice Isoforms | Functional Consequence | Disease Relevance |
|----------|-----------------|------------------------|-------------------|
| **IL-33** | Full-length vs Δexon3,4 | Δexon3,4 is cytoplasmic, more readily secreted, retains signaling | Associated with chronic stable asthma type 2 inflammation |
| **IL-15** | Long signal peptide (LSP) vs Short signal peptide (SSP) | SSP can be secreted and bioactive with IL-15Rα; tunes magnitude/duration of IL-15 activity | NK and T-cell signaling modulation |
| **VEGF-A** | VEGF165 vs VEGF165b | VEGF165b binds VEGFR2 but antagonizes VEGF165-driven angiogenesis | Tumor angiogenesis; affects anti-VEGF therapy |
| **CXCL12** | α, β, γ isoforms | CXCL12γ has high heparan sulfate affinity, strong surface retention | Gradient formation and bioactivity |
| **TRAIL** | TRAIL-α vs TRAIL-β/γ vs TRAIL-short | TRAIL-β/γ lose pro-apoptotic potential; TRAIL-short is dominant negative | HIV immune evasion; cancer therapy resistance |
| **IL-4** | Canonical vs IL-4δ2 | IL-4δ2 binds receptor but inhibits IL-4 signaling | Built-in RNA-level brake on IL-4 activity |
| **IL-18BP** | IL-18BPa/c vs IL-18BPb/d | Only a/c isoforms bind and neutralize IL-18 | Cytokine inhibition capacity |

### 4.1 Secreted Protein and Enzyme Examples

| Protein | Splice Isoforms | Functional Consequence |
|---------|-----------------|------------------------|
| **Cathepsin B** | Canonical vs CB(-2,3) missing exons 2-3 | CB(-2,3) bypasses secretory pathway, enables internal translation start, yields active enzyme with altered localization |
| **PCSK5/6** | PC5/6A vs PC5/6B | Splicing determines secreted vs membrane-associated forms, reshaping proteolytic maturation landscape |
| **ADAMTS17** | Spacer-domain variants | Affects secretion and autoproteolytic activity |

---

## 5. Mechanistic Comparison: Proteolytic vs. Splicing-Based Regulation

| Aspect | Proteolytic Activation | Alternative Splicing |
|--------|----------------------|---------------------|
| **Speed** | Fast (minutes to hours) | Slow (transcription/translation-dependent) |
| **Reversibility** | Typically irreversible | Reversible via splicing factor changes |
| **Spatial control** | Requires protease/substrate co-localization | Cell-autonomous, programmed by splicing factors |
| **Control point** | Post-translational (protein level) | Pre-translational (mRNA level) |
| **Outcome** | Activates existing precursor pool | Changes which protein architectures are produced |

### 5.1 Functional Motifs Frequently Targeted by Splicing

Across the examples, several recurring "splicing targets" stand out:

1. **Signal peptides / trafficking determinants**: IL-15 LSP vs SSP; Cathepsin B variants affecting secretory routing
2. **C-terminal tails controlling diffusion vs retention**: CXCL12γ heparan sulfate retention; VEGF-A terminal exon switching
3. **Domains required for binding/neutralization**: IL-18BP Ig-like domain integrity
4. **Loss-of-function or dominant-negative variants**: IL-4δ2 antagonism; TRAIL splice variants
5. **Protease activation logic**: PCSK family splicing; ADAMTS17 spacer domain variants

---

## 6. Disease Relevance of Cytokine Splicing

### 6.1 Airway Disease and Epithelial Cytokines (IL-33)

A notably direct disease-linked example is **IL-33 in asthma**: alternative splicing in airway epithelial cells produces an IL-33 isoform (Δexon3,4) that is cytoplasmic, more readily secreted, and associated with type 2 inflammation in chronic stable asthma, whereas full-length IL-33 is not similarly associated in that setting. This implies that chronic disease may rely on an isoform repertoire that enables extracellular signaling without widespread necrosis—conceptually different from an "alarmin released only upon cell death" model.

**Reference**: Gordon ED, Simpson LJ, Rios CL, et al. Alternative splicing of interleukin-33 and type 2 inflammation in asthma. *Proc Natl Acad Sci USA*. 2016;113(31):8765-8770. doi: 10.1073/pnas.1601914113

### 6.2 Cancer Angiogenesis and Therapy Relevance (VEGF-A)

**VEGF-A** provides one of the clearest examples where splicing generates isoforms with opposing angiogenic outputs. The VEGF165b isoform is an endogenous inhibitory/anti-angiogenic splice variant.

Beyond biology, isoform balance may influence therapy: work in colorectal carcinoma suggests VEGF165b binds bevacizumab similarly to VEGF165 and that shifts in VEGFxxx vs VEGFxxxb balance could affect predicted responses because the antibody does not inherently distinguish "pro-angiogenic" vs "anti-angiogenic" families.

**Reference**: Bates DO, Cui TG, Doughty JM, et al. VEGF165b, an inhibitory splice variant of vascular endothelial growth factor, is down-regulated in renal cell carcinoma. *Cancer Res*. 2002;62(14):4123-4131. PMID: 12124351

### 6.3 Immune Evasion and Pathogen-Associated Mechanisms (TRAIL)

TRAIL splice variants illustrate how splicing can modulate "death cytokine" function. Canonical TRAIL can induce apoptosis via its death receptors, but alternative splice variants (TRAIL-β, TRAIL-γ) lack apoptotic potential.

A particularly strong mechanistic disease-context connection is **TRAIL-short**, reported as an HIV-induced splice variant that antagonizes TRAIL receptor signaling and correlates with viral load, providing an example of splicing-mediated evasion of cytokine-driven cell death. TRAIL-short has also been implicated in cancer therapy resistance.

**References**:
- Schnepple DJ, Shepard B, Bren GD, et al. Isolation of a TRAIL antagonist from the serum of HIV-infected patients. *J Biol Chem*. 2011;286(41):35742-35754. doi: 10.1074/jbc.M111.274639
- Aboulnasr F, Krogman A, Graham RP, et al. Human cancers express TRAILshort, a dominant negative TRAIL splice variant, which impairs immune effector cell killing of tumor cells. *Clin Cancer Res*. 2020;26(21):5759-5771. doi: 10.1158/1078-0432.CCR-20-0251

### 6.4 Immune Homeostasis: Tuning Cytokine Magnitude (IL-15, IL-4)

**IL-15** shows a splicing-driven control point that does not merely "turn IL-15 on/off," but can tune *how much* bioactive IL-15 is produced and how stable the IL-15/IL-15Rα complex is—allowing fine-grained modulation of NK and T-cell–relevant signaling.

**Reference**: Tagaya Y, Kurys G, Thies TA, et al. Generation of secretable and nonsecretable interleukin 15 isoforms through alternate usage of signal peptides. *Proc Natl Acad Sci USA*. 1997;94(26):14444-14449. doi: 10.1073/pnas.94.26.14444

**IL-4δ2** demonstrates that splicing can generate a **ligand antagonist** (or functional inhibitor) that competes at the receptor level and is differentially expressed across tissues, implying a built-in RNA-level "brake" on IL-4 activity in certain contexts.

**Reference**: Atamas SP, Choi J, Yurovsky VV, White B. An alternative splice variant of human IL-4, IL-4 delta 2, inhibits IL-4-stimulated T cell proliferation. *J Immunol*. 1996;156(2):435-441. PMID: 8543791

### 6.5 Chemokine Gradients and Cancer Metastasis (CXCL12)

CXCL12/SDF-1 isoforms demonstrate how alternative splicing creates chemokines with different extracellular matrix retention properties. CXCL12γ has an extended C-terminus with high heparan sulfate binding affinity, leading to strong surface retention and steep local gradients, while CXCL12α is more diffusible.

**Reference**: Dimberg J, Skarstedt M, Löfgren S, Zar N, Matussek A. Analysis of the expression of SDF-1 splicing variants in human colorectal cancer and normal mucosa tissues. *Oncol Lett*. 2016;11(3):1679-1684. doi: 10.3892/ol.2016.4139

---

## 7. Implications for Spatial Transcriptomics Analysis

The existence of functionally distinct cytokine isoforms with **condition-, tissue-, or disease-associated expression shifts** motivates spatial isoform analysis:

1. **Tumor Microenvironment**: Different cytokine isoforms may predominate at tumor core vs. invasive margin (e.g., pro-angiogenic VEGF165 in tumor core, anti-angiogenic VEGF165b at margin)

2. **Immune Infiltration Patterns**: Antagonist vs. agonist isoform ratios may correlate with immune cell distribution (e.g., TRAIL-short in immunosuppressive regions)

3. **Hypoxia Gradients**: Metabolic stress may drive isoform switching via hypoxia-responsive splicing factors

4. **Cell-Cell Communication**: Matrix-bound vs. diffusible isoforms (e.g., CXCL12γ vs CXCL12α) create different signaling gradients

5. **Therapeutic Response Prediction**: Spatial mapping of agonist/antagonist isoform ratios may predict responses to targeted therapies

**SPLISOSM enables detection of these spatially variable processing (SVP) patterns without requiring prior knowledge of cell type boundaries.**

---

## 8. Key References

### Cytokine and Secreted Protein Splicing

1. Gordon ED, et al. Alternative splicing of interleukin-33 and type 2 inflammation in asthma. *PNAS* 2016;113:8765-8770. doi: 10.1073/pnas.1601914113

2. Tagaya Y, et al. Generation of secretable and nonsecretable interleukin 15 isoforms through alternate usage of signal peptides. *PNAS* 1997;94:14444-14449. doi: 10.1073/pnas.94.26.14444

3. Bates DO, et al. VEGF165b, an inhibitory splice variant of vascular endothelial growth factor, is down-regulated in renal cell carcinoma. *Cancer Res* 2002;62:4123-4131. PMID: 12124351

4. Atamas SP, et al. An alternative splice variant of human IL-4, IL-4 delta 2, inhibits IL-4-stimulated T cell proliferation. *J Immunol* 1996;156:435-441. PMID: 8543791

5. Dimberg J, et al. Analysis of the expression of SDF-1 splicing variants in human colorectal cancer and normal mucosa tissues. *Oncol Lett* 2016;11:1679-1684. doi: 10.3892/ol.2016.4139

6. Schnepple DJ, et al. Isolation of a TRAIL antagonist from the serum of HIV-infected patients. *J Biol Chem* 2011;286:35742-35754. doi: 10.1074/jbc.M111.274639

7. Aboulnasr F, et al. Human cancers express TRAILshort, a dominant negative TRAIL splice variant. *Clin Cancer Res* 2020;26:5759-5771. doi: 10.1158/1078-0432.CCR-20-0251

8. Kim HS, et al. Identification of novel IL-18 binding protein isoforms. *J Immunol* 2000;164:5905-5910.

9. Seidah NG, Prat A. The biology and therapeutic targeting of the proprotein convertases. *Nat Rev Drug Discov* 2012;11:367-383.

### Spatial Transcriptomics Methods

10. Su J, et al. Mapping isoforms and regulatory mechanisms from spatial transcriptomics data with SPLISOSM. *Nature Biotechnology* 2026. doi: 10.1038/s41587-025-02965-6

---

*Document Version: 2.0 | Last Updated: February 2026*
*Covers: Cytokine splicing biology, proteolytic activation mechanisms, disease relevance, spatial transcriptomics implications*

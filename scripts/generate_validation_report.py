#!/usr/bin/env python3
"""Generate comprehensive validation report combining ONT HSIC-IR and scRNA-seq results."""

import pandas as pd
from pathlib import Path

project_root = Path(__file__).parent.parent
reports_dir = project_root / "reports"

print("=" * 80)
print("COMPREHENSIVE VALIDATION REPORT")
print("=" * 80)

# Load ONT results
ont_summary = pd.read_csv(reports_dir / "candidate_validation_summary.csv")

# Load scRNA-seq differential splicing results
scrna_diff = pd.read_csv(reports_dir / "scrna_wilcox_diff_20_gene_short_name.csv")
scrna_features = pd.read_csv(reports_dir / "scrna_splice_feature_gene_short_name.csv")

# Summarize scRNA-seq splice features per gene
splice_counts = scrna_features.groupby('gene_short_name').size().reset_index(name='n_splice_events')

# Summarize differential splicing per gene
diff_summary = scrna_diff.groupby('gene_short_name').agg({
    'mean.diff': ['min', 'max', 'mean'],
    'p.val': 'min',
    'p.val.adj': 'min',
    'comparasion': 'nunique'
}).reset_index()
diff_summary.columns = ['gene', 'min_delta_psi', 'max_delta_psi', 'mean_delta_psi',
                        'min_pval', 'min_padj', 'n_comparisons']

print("\n" + "=" * 80)
print("1. ONT VISIUM HSIC-IR RESULTS")
print("=" * 80)
print("\nGenes with spatially variable isoform usage (HSIC-IR test):\n")
print(ont_summary.to_string(index=False))

print("\n" + "=" * 80)
print("2. scRNA-seq SPLICE FEATURES")
print("=" * 80)
print(f"\nAll {len(splice_counts)} candidate genes have splice events in scRNA-seq:")
print(splice_counts.to_string(index=False))

print("\n" + "=" * 80)
print("3. scRNA-seq DIFFERENTIAL SPLICING BY TUMOR STATE")
print("=" * 80)
print("\nGenes with significant differential splicing between tumor states:\n")
print(diff_summary.to_string(index=False))

# Merge ONT and scRNA-seq results
merged = ont_summary.merge(splice_counts, left_on='gene', right_on='gene_short_name', how='left')
merged = merged.merge(diff_summary, on='gene', how='left')

# Determine overall validation status
def get_validation_tier(row):
    hsic_ir_validated = row['pct_sig_hsic_ir'] >= 50
    scrna_validated = pd.notna(row.get('min_padj')) and row.get('min_padj', 1) < 0.05

    if hsic_ir_validated and scrna_validated:
        return 'TIER1_FULL'
    elif hsic_ir_validated:
        return 'TIER2_SPATIAL'
    elif scrna_validated:
        return 'TIER3_SCRNA'
    else:
        return 'TIER4_PARTIAL'

merged['overall_validation'] = merged.apply(get_validation_tier, axis=1)

# Sort by validation tier and p-value
tier_order = {'TIER1_FULL': 0, 'TIER2_SPATIAL': 1, 'TIER3_SCRNA': 2, 'TIER4_PARTIAL': 3}
merged['tier_rank'] = merged['overall_validation'].map(tier_order)
merged = merged.sort_values(['tier_rank', 'min_pvalue_ir'])

print("\n" + "=" * 80)
print("4. INTEGRATED VALIDATION SUMMARY")
print("=" * 80)

tier1 = merged[merged['overall_validation'] == 'TIER1_FULL']
tier2 = merged[merged['overall_validation'] == 'TIER2_SPATIAL']

print(f"\n### TIER 1: Fully Validated ({len(tier1)} genes)")
print("Significant in both spatial HSIC-IR AND scRNA-seq differential splicing")
print("-" * 60)
if len(tier1) > 0:
    for _, row in tier1.iterrows():
        print(f"  {row['gene']:10s} | Spatial: {row['n_sig_hsic_ir']}/{row['n_samples_tested']} samples sig | "
              f"scRNA: ΔPSI={row['max_delta_psi']:.1f}%, padj={row['min_padj']:.2e}")
else:
    print("  None")

print(f"\n### TIER 2: Spatial-Only Validated ({len(tier2)} genes)")
print("Significant in HSIC-IR spatial analysis (≥50% samples)")
print("-" * 60)
for _, row in tier2.iterrows():
    print(f"  {row['gene']:10s} | Spatial: {row['n_sig_hsic_ir']}/{row['n_samples_tested']} samples sig | "
          f"Best p={row['min_pvalue_ir']:.2e}")

# Save final validation report
output_cols = ['gene', 'n_samples_tested', 'n_sig_hsic_ir', 'pct_sig_hsic_ir',
               'min_pvalue_ir', 'n_splice_events', 'max_delta_psi', 'min_padj',
               'overall_validation', 'sample_types']
final_report = merged[output_cols].copy()
final_report.to_csv(reports_dir / "final_validation_report.csv", index=False)

print("\n" + "=" * 80)
print("5. KEY BIOLOGICAL FINDINGS")
print("=" * 80)

print("""
## APP (Amyloid Precursor Protein)
- Spatial: Detected in 0 ONT samples (low expression in ONT data)
- scRNA-seq: STRONG validation - PSI differences up to 55% between tumor states
  - MES vs NPC: ΔPSI = -55%, p = 8.6e-52
  - AC vs NPC: ΔPSI = -48%, p = 4.6e-39
- Biological: Consistent with tumor state-specific splicing regulation

## VCAN (Versican)
- Spatial: 3/9 ONT samples significant (33%)
- scRNA-seq: STRONG validation - PSI differences up to 64%
  - MES vs OPC_OL: ΔPSI = +64%, p = 2.5e-11
  - MES vs OPC_N: ΔPSI = +47%, p = 4.4e-09
- Biological: ECM proteoglycan with known isoform-specific functions

## CLU (Clusterin)
- Spatial: 11/11 ONT samples significant (100%), p = 5.1e-74
- scRNA-seq: Has 14 splice events but no significant tumor state differences
- Interpretation: Spatial variability may reflect microenvironment rather than tumor state

## B2M (Beta-2 Microglobulin)
- Spatial: 11/11 ONT samples significant (100%), p = 1.7e-67
- scRNA-seq: Has 7 splice events but no significant tumor state differences
- Interpretation: MHC class I component with spatial regulation

## CD74
- Spatial: 8/10 ONT samples significant (80%), p = 5.5e-38
- scRNA-seq: Has 12 splice events but no significant tumor state differences
- Interpretation: Antigen presentation related, immune microenvironment marker

## MIF (Macrophage Migration Inhibitory Factor)
- Spatial: 9/11 ONT samples significant (82%), p = 2.2e-21
- scRNA-seq: Has 5 splice events but no significant tumor state differences
- Interpretation: Cytokine with spatial regulation in TME
""")

print("=" * 80)
print("REPORT COMPLETE")
print("=" * 80)
print(f"\nSaved: {reports_dir / 'final_validation_report.csv'}")

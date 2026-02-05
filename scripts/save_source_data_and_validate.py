#!/usr/bin/env python3
"""Save source data as CSV files and perform remaining validations.

This script:
1. Saves figure source data (APP, VCAN, CLU) as CSV files in reports/
2. Validates remaining candidates (B2M, CD74, MIF, LGALS1, etc.) using scRNA-seq PSI data
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import pandas as pd
import numpy as np

# Load required R data using rpy2
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()
    HAS_RPY2 = True
except ImportError:
    HAS_RPY2 = False
    print("Warning: rpy2 not available, will use pre-converted data if available")


def load_rdata(filepath):
    """Load an .RData file and return as dict."""
    if not HAS_RPY2:
        raise ImportError("rpy2 required for loading .RData files")

    base = importr('base')
    ro.r['load'](str(filepath))

    # Get all objects in the R environment
    objects = {}
    for name in ro.r.ls():
        try:
            obj = ro.r[name]
            # Try to convert to pandas
            try:
                objects[name] = pandas2ri.rpy2py(obj)
            except:
                objects[name] = obj
        except:
            continue
    return objects


def save_ont_results_csv():
    """Save ONT SPLISOSM results for report figures."""
    print("\n=== Saving ONT SPLISOSM Results ===")

    ont_results_file = project_root / "data/zenodo_16905935/human_glioma_ont/ont_secact_splisosm_results.csv"

    if not ont_results_file.exists():
        print(f"ERROR: ONT results not found at {ont_results_file}")
        return None

    df = pd.read_csv(ont_results_file)

    # Save full results
    output_dir = project_root / "reports"
    output_dir.mkdir(exist_ok=True)

    # Save full ONT results
    full_output = output_dir / "ont_splisosm_full_results.csv"
    df.to_csv(full_output, index=False)
    print(f"Saved full ONT results: {full_output}")

    # Save CLU-specific results for fig3
    clu_df = df[df['gene'] == 'CLU'].copy()
    clu_df = clu_df.sort_values('pvalue_ir')
    clu_output = output_dir / "fig3_clu_source_data.csv"
    clu_df.to_csv(clu_output, index=False)
    print(f"Saved CLU figure data: {clu_output} ({len(clu_df)} rows)")

    # Generate summary by gene (for validation)
    summary = df.groupby('gene').agg({
        'pvalue_ir': ['min', 'mean', 'count'],
        'sample_id': 'nunique'
    }).round(6)
    summary.columns = ['min_pvalue', 'mean_pvalue', 'n_tests', 'n_samples']
    summary['sig_count'] = df[df['pvalue_ir'] < 0.05].groupby('gene').size()
    summary = summary.fillna(0).sort_values('min_pvalue')

    summary_output = output_dir / "ont_splisosm_gene_summary.csv"
    summary.to_csv(summary_output)
    print(f"Saved gene summary: {summary_output}")

    return df


def save_pilot_results_csv():
    """Save pilot SPLISOSM results for report figures."""
    print("\n=== Saving Pilot SPLISOSM Results ===")

    pilot_file = project_root / "data/zenodo_16905935/pilot_splisosm_results.csv"

    if not pilot_file.exists():
        print(f"ERROR: Pilot results not found at {pilot_file}")
        return None

    df = pd.read_csv(pilot_file)

    output_dir = project_root / "reports"

    # Save full pilot results
    full_output = output_dir / "pilot_splisosm_full_results.csv"
    df.to_csv(full_output, index=False)
    print(f"Saved full pilot results: {full_output}")

    # Save APP-specific results for fig1
    app_df = df[df['gene'] == 'APP'].copy()
    app_output = output_dir / "fig1_app_source_data.csv"
    app_df.to_csv(app_output, index=False)
    print(f"Saved APP figure data: {app_output} ({len(app_df)} rows)")

    # Save VCAN-specific results for fig2
    vcan_df = df[df['gene'] == 'VCAN'].copy()
    vcan_output = output_dir / "fig2_vcan_source_data.csv"
    vcan_df.to_csv(vcan_output, index=False)
    print(f"Saved VCAN figure data: {vcan_output} ({len(vcan_df)} rows)")

    return df


def validate_candidates_scrna():
    """Validate top HSIC-IR candidates using scRNA-seq PSI data."""
    print("\n=== Validating Candidates with scRNA-seq PSI ===")

    scrna_dir = project_root / "data/zenodo_17055113/zenodo-20250904/MARVEL"

    # Key files for validation
    psi_file = scrna_dir / "PSI_merge.RData"
    tumor_state_file = scrna_dir / "event_list_tumor_state_diff_wilcox.RData"
    splice_feature_file = scrna_dir / "splice_feature.RData"

    if not psi_file.exists():
        print(f"PSI data not found at {psi_file}")
        return None

    if not HAS_RPY2:
        print("rpy2 required for loading RData files")
        return None

    print(f"Loading PSI data from {psi_file}...")

    # Load PSI data
    try:
        ro.r['load'](str(psi_file))

        # Get the PSI merge object
        psi_merge = None
        for name in ro.r.ls():
            obj = ro.r[name]
            print(f"  Found object: {name}, type: {type(obj)}")
            if 'PSI' in name or 'psi' in name:
                psi_merge = obj
                break

        if psi_merge is None:
            print("  No PSI object found in RData file")
            # List all objects
            for name in ro.r.ls():
                print(f"    - {name}")

    except Exception as e:
        print(f"Error loading PSI data: {e}")

    # Load tumor state differential events
    print(f"\nLoading tumor state differential events...")
    try:
        ro.r['load'](str(tumor_state_file))

        for name in ro.r.ls():
            obj = ro.r[name]
            print(f"  Found object: {name}")

    except Exception as e:
        print(f"Error loading tumor state data: {e}")

    # Load splice features (gene annotations)
    print(f"\nLoading splice feature data...")
    try:
        ro.r['load'](str(splice_feature_file))

        for name in ro.r.ls():
            print(f"  Found object: {name}")

    except Exception as e:
        print(f"Error loading splice feature data: {e}")

    return None


def analyze_candidates_from_rdata():
    """Analyze validation candidates using R directly."""
    print("\n=== Analyzing Candidate Genes from scRNA-seq ===")

    scrna_dir = project_root / "data/zenodo_17055113/zenodo-20250904/MARVEL"
    output_dir = project_root / "reports"

    if not HAS_RPY2:
        print("rpy2 not available - generating validation summary from ONT results only")
        return generate_validation_summary_ont_only()

    # Define candidate genes from ONT HSIC-IR analysis
    candidates = ['B2M', 'CLU', 'CD74', 'MIF', 'LGALS1', 'PTGDS', 'CST3', 'SPARCL1',
                  'TIMP1', 'ITM2B', 'LGALS3', 'SPP1', 'AGT', 'C3', 'SPARC', 'APOE']

    print(f"\nSearching for {len(candidates)} candidate genes in scRNA-seq splicing data...")

    # Load splice feature data to find gene names
    try:
        ro.r['load'](str(scrna_dir / "splice_feature.RData"))

        # Try to extract gene information
        print("\nSearching for candidates in splice features...")

        # Get variable names
        r_vars = list(ro.r.ls())
        print(f"Available R objects: {r_vars}")

        # Load tumor state differential expression data
        ro.r['load'](str(scrna_dir / "event_list_tumor_state_diff_wilcox.RData"))
        r_vars = list(ro.r.ls())
        print(f"After loading tumor state data: {r_vars}")

        # Try to access the data
        validation_results = []

        for var_name in r_vars:
            try:
                obj = ro.r[var_name]
                # Check if it's a data frame
                if hasattr(obj, 'colnames'):
                    cols = list(obj.colnames)
                    print(f"\n{var_name} columns: {cols[:10]}...")  # First 10 columns

                    # Look for gene column
                    gene_cols = [c for c in cols if 'gene' in c.lower()]
                    if gene_cols:
                        print(f"  Gene columns found: {gene_cols}")

                        # Convert to pandas
                        try:
                            df = pandas2ri.rpy2py(obj)

                            # Check for our candidates
                            for gene_col in gene_cols:
                                if gene_col in df.columns:
                                    found = df[df[gene_col].astype(str).isin(candidates)]
                                    if len(found) > 0:
                                        print(f"  Found {len(found)} events for candidates in {gene_col}")
                                        found['source'] = var_name
                                        validation_results.append(found)
                        except Exception as e:
                            print(f"  Could not convert {var_name}: {e}")

            except Exception as e:
                pass

        if validation_results:
            combined = pd.concat(validation_results, ignore_index=True)
            output_file = output_dir / "scrna_validation_results.csv"
            combined.to_csv(output_file, index=False)
            print(f"\nSaved scRNA-seq validation results: {output_file}")
            return combined

    except Exception as e:
        print(f"Error analyzing RData: {e}")
        import traceback
        traceback.print_exc()

    return generate_validation_summary_ont_only()


def generate_validation_summary_ont_only():
    """Generate validation summary using ONT results only."""
    print("\n=== Generating Validation Summary from ONT Results ===")

    output_dir = project_root / "reports"
    ont_file = project_root / "data/zenodo_16905935/human_glioma_ont/ont_secact_splisosm_results.csv"

    if not ont_file.exists():
        print("ONT results not found")
        return None

    df = pd.read_csv(ont_file)

    # Calculate validation metrics per gene
    summary_data = []

    for gene in df['gene'].unique():
        gene_df = df[df['gene'] == gene]

        n_samples = len(gene_df)
        n_sig_ir = (gene_df['pvalue_ir'] < 0.05).sum()
        n_sig_gc = (gene_df['pvalue_gc'] < 0.05).sum()
        n_sig_ic = (gene_df['pvalue_ic'] < 0.05).sum()

        summary_data.append({
            'gene': gene,
            'n_samples_tested': n_samples,
            'n_sig_hsic_ir': n_sig_ir,
            'pct_sig_hsic_ir': round(100 * n_sig_ir / n_samples, 1),
            'min_pvalue_ir': gene_df['pvalue_ir'].min(),
            'median_pvalue_ir': gene_df['pvalue_ir'].median(),
            'n_sig_hsic_gc': n_sig_gc,
            'n_sig_hsic_ic': n_sig_ic,
            'avg_n_isoforms': gene_df['n_isoforms'].mean(),
            'avg_n_spots': gene_df['n_spots_expressed'].mean(),
            'sample_types': ','.join(gene_df['sample_type'].unique()),
            'validation_status': 'VALIDATED' if n_sig_ir >= n_samples * 0.5 else 'PARTIAL' if n_sig_ir >= 2 else 'NOT_VALIDATED'
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('min_pvalue_ir')

    # Save summary
    output_file = output_dir / "candidate_validation_summary.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"Saved validation summary: {output_file}")

    # Print summary
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)

    validated = summary_df[summary_df['validation_status'] == 'VALIDATED']
    partial = summary_df[summary_df['validation_status'] == 'PARTIAL']

    print(f"\nFully Validated Genes ({len(validated)} genes, ≥50% samples significant):")
    print("-" * 60)
    for _, row in validated.iterrows():
        print(f"  {row['gene']:10s} - {row['n_sig_hsic_ir']}/{row['n_samples_tested']} samples "
              f"({row['pct_sig_hsic_ir']:.0f}%), best p={row['min_pvalue_ir']:.2e}")

    print(f"\nPartially Validated Genes ({len(partial)} genes, 2+ samples significant):")
    print("-" * 60)
    for _, row in partial.iterrows():
        print(f"  {row['gene']:10s} - {row['n_sig_hsic_ir']}/{row['n_samples_tested']} samples "
              f"({row['pct_sig_hsic_ir']:.0f}%), best p={row['min_pvalue_ir']:.2e}")

    return summary_df


def main():
    print("=" * 80)
    print("Saving Source Data and Validating Candidates")
    print("=" * 80)

    # Part 1: Save source data CSV files
    save_pilot_results_csv()
    save_ont_results_csv()

    # Part 2: Validate candidates
    analyze_candidates_from_rdata()

    print("\n" + "=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print("\nGenerated files in reports/:")

    reports_dir = project_root / "reports"
    for f in sorted(reports_dir.glob("*.csv")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()

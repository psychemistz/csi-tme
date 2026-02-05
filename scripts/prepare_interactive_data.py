#!/usr/bin/env python3
"""
Prepare data for interactive HTML dashboard.

Reads all CSV files and images, converts to JSON format with base64-encoded images.
Output: reports/dashboard_data.json
"""

import json
import base64
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
REPORTS_DIR = PROJECT_ROOT / "reports"
OUTPUT_FILE = REPORTS_DIR / "dashboard_data.json"


def read_csv_to_records(filepath):
    """Read CSV file and return as list of dicts."""
    if not filepath.exists():
        print(f"  Warning: {filepath} not found")
        return []
    df = pd.read_csv(filepath)
    # Convert NaN to None for JSON serialization
    df = df.where(pd.notnull(df), None)
    return df.to_dict(orient='records')


def encode_image_base64(filepath):
    """Read image and encode as base64 string."""
    if not filepath.exists():
        print(f"  Warning: {filepath} not found")
        return None
    with open(filepath, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode('utf-8')


def main():
    print("Preparing interactive dashboard data...")

    dashboard_data = {
        'metadata': {
            'generated': pd.Timestamp.now().isoformat(),
            'project': 'CSI-TME: Spatial Isoform Validation'
        },
        'csv_data': {},
        'images': {},
        'summary_stats': {}
    }

    # === CSV FILES ===
    csv_files = {
        # HIGH PRIORITY - Main results
        'ont_splisosm_full_results': REPORTS_DIR / 'ont_splisosm_full_results.csv',
        'sr_splisosm_full_results': REPORTS_DIR / 'sr_splisosm_full_results.csv',
        'final_validation_report': REPORTS_DIR / 'final_validation_report.csv',
        'scrna_full_psi_differential': REPORTS_DIR / 'scrna_full_psi_differential.csv',

        # MEDIUM PRIORITY - Validation
        'spacet_isoform_correlations': REPORTS_DIR / 'spacet_isoform_correlations.csv',
        'spacet_correlation_summary': REPORTS_DIR / 'spacet_correlation_summary.csv',
        'scrna_validation_corrected': REPORTS_DIR / 'scrna_validation_corrected.csv',
        'scrna_validation_all_genes': REPORTS_DIR / 'scrna_validation_all_genes.csv',
        'candidate_validation_summary': REPORTS_DIR / 'candidate_validation_summary.csv',
        'ont_splisosm_gene_summary': REPORTS_DIR / 'ont_splisosm_gene_summary.csv',
        'sr_splisosm_gene_summary': REPORTS_DIR / 'sr_splisosm_gene_summary.csv',

        # LOW PRIORITY - Figure source data
        'fig1_app_source_data': REPORTS_DIR / 'fig1_app_source_data.csv',
        'fig2_vcan_source_data': REPORTS_DIR / 'fig2_vcan_source_data.csv',
        'fig3_clu_source_data': REPORTS_DIR / 'fig3_clu_source_data.csv',
    }

    print("\nLoading CSV files...")
    for name, filepath in csv_files.items():
        print(f"  Loading {name}...")
        records = read_csv_to_records(filepath)
        dashboard_data['csv_data'][name] = records
        print(f"    -> {len(records)} records")

    # === IMAGES ===
    print("\nEncoding images as base64...")

    # Main figures
    main_figures = [
        ('fig1_app_psi_by_tumor_state', 'fig1_app_psi_by_tumor_state.png'),
        ('fig1_app_psi_boxplot', 'fig1_app_psi_boxplot.png'),
        ('fig2_vcan_psi_by_tumor_state', 'fig2_vcan_psi_by_tumor_state.png'),
        ('fig2_vcan_psi_boxplot', 'fig2_vcan_psi_boxplot.png'),
        ('fig3_clu_visium_hsic_ir', 'fig3_clu_visium_hsic_ir.png'),
        ('fig4_spatial_isoform_switching', 'fig4_spatial_isoform_switching.png'),
        ('fig5_spatial_isoform_ZH1019T1', 'fig5_spatial_isoform_ZH1019T1.png'),
    ]

    dashboard_data['images']['main_figures'] = {}
    for name, filename in main_figures:
        filepath = REPORTS_DIR / filename
        print(f"  Encoding {filename}...")
        b64 = encode_image_base64(filepath)
        if b64:
            dashboard_data['images']['main_figures'][name] = {
                'filename': filename,
                'data': b64,
                'type': 'image/png'
            }

    # SpaCET deconvolution plots
    spacet_files = list(REPORTS_DIR.glob('spacet_*_macro.png'))
    dashboard_data['images']['spacet_plots'] = {}
    for filepath in spacet_files:
        sample = filepath.stem.replace('spacet_', '').replace('_macro', '')
        print(f"  Encoding SpaCET {sample}...")
        b64 = encode_image_base64(filepath)
        if b64:
            dashboard_data['images']['spacet_plots'][sample] = {
                'filename': filepath.name,
                'data': b64,
                'type': 'image/png'
            }

    # Spatial isoform plots
    spatial_plots_dir = REPORTS_DIR / 'spatial_plots'
    if spatial_plots_dir.exists():
        spatial_files = list(spatial_plots_dir.glob('spatial_isoform_*.png'))
        dashboard_data['images']['spatial_plots'] = {}
        for filepath in spatial_files:
            sample = filepath.stem.replace('spatial_isoform_', '')
            print(f"  Encoding spatial {sample}...")
            b64 = encode_image_base64(filepath)
            if b64:
                dashboard_data['images']['spatial_plots'][sample] = {
                    'filename': filepath.name,
                    'data': b64,
                    'type': 'image/png'
                }

    # === COMPUTE SUMMARY STATISTICS ===
    print("\nComputing summary statistics...")

    # ONT results summary
    ont_data = dashboard_data['csv_data'].get('ont_splisosm_full_results', [])
    if ont_data:
        ont_df = pd.DataFrame(ont_data)
        genes_ont = ont_df['gene'].nunique()
        samples_ont = ont_df['sample_id'].nunique()
        sig_ont = len(ont_df[ont_df['pvalue_ir'] < 0.05])
        dashboard_data['summary_stats']['ont'] = {
            'n_genes': genes_ont,
            'n_samples': samples_ont,
            'n_significant': sig_ont,
            'best_pvalue': ont_df['pvalue_ir'].min() if 'pvalue_ir' in ont_df else None
        }

    # SR results summary
    sr_data = dashboard_data['csv_data'].get('sr_splisosm_full_results', [])
    if sr_data:
        sr_df = pd.DataFrame(sr_data)
        genes_sr = sr_df['gene'].nunique()
        samples_sr = sr_df['sample'].nunique()
        # Use pvalue_ir for the new comprehensive SR results
        pval_col = 'pvalue_ir' if 'pvalue_ir' in sr_df.columns else 'pvalue'
        sig_sr = len(sr_df[sr_df[pval_col] < 0.05])
        dashboard_data['summary_stats']['sr'] = {
            'n_genes': genes_sr,
            'n_samples': samples_sr,
            'n_significant': sig_sr,
            'best_pvalue': sr_df[pval_col].min() if pval_col in sr_df else None
        }

    # Validation summary
    val_data = dashboard_data['csv_data'].get('final_validation_report', [])
    if val_data:
        val_df = pd.DataFrame(val_data)
        # Handle both old and new column names
        tier_col = 'validation_tier' if 'validation_tier' in val_df.columns else 'overall_validation'
        tier_counts = val_df[tier_col].value_counts().to_dict()
        dashboard_data['summary_stats']['validation'] = {
            'n_genes_tested': len(val_df),
            'tier_counts': tier_counts
        }

    # scRNA validation summary
    scrna_data = dashboard_data['csv_data'].get('scrna_full_psi_differential', [])
    if scrna_data:
        scrna_df = pd.DataFrame(scrna_data)
        n_events = len(scrna_df)
        sig_events = len(scrna_df[(scrna_df['padj'] < 0.05) & (scrna_df['delta_psi'].abs() > 0.1)])
        dashboard_data['summary_stats']['scrna'] = {
            'n_comparisons': n_events,
            'n_significant': sig_events,
            'genes_with_sig': scrna_df[scrna_df['padj'] < 0.05]['gene'].nunique()
        }

    # Gene-level aggregation for search table
    print("\nAggregating gene-level summary...")
    gene_summary = aggregate_gene_summary(dashboard_data['csv_data'])
    dashboard_data['gene_summary'] = gene_summary

    # Sample-level aggregation
    print("Aggregating sample-level summary...")
    sample_summary = aggregate_sample_summary(dashboard_data['csv_data'])
    dashboard_data['sample_summary'] = sample_summary

    # === WRITE OUTPUT ===
    print(f"\nWriting output to {OUTPUT_FILE}...")
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(dashboard_data, f)

    file_size_mb = OUTPUT_FILE.stat().st_size / (1024 * 1024)
    print(f"Done! Output size: {file_size_mb:.2f} MB")

    # Summary
    print("\n=== Data Summary ===")
    print(f"CSV files loaded: {len(dashboard_data['csv_data'])}")
    print(f"Main figures: {len(dashboard_data['images'].get('main_figures', {}))}")
    print(f"SpaCET plots: {len(dashboard_data['images'].get('spacet_plots', {}))}")
    print(f"Spatial plots: {len(dashboard_data['images'].get('spatial_plots', {}))}")
    print(f"Gene summaries: {len(dashboard_data.get('gene_summary', []))}")
    print(f"Sample summaries: {len(dashboard_data.get('sample_summary', []))}")


def aggregate_gene_summary(csv_data):
    """Create gene-level summary table from validation report."""
    gene_info = []

    # Use final_validation_report as the primary source (has all IR, GC, IC data)
    for row in csv_data.get('final_validation_report', []):
        gene = row.get('gene')
        if not gene:
            continue

        info = {
            'gene': gene,
            'validation_tier': row.get('validation_tier'),
            'platforms': row.get('platforms'),
            'total_samples': row.get('total_samples_tested', 0),
            # IR stats
            'ir_sig': row.get('total_ir_sig', 0),
            'ir_best_pval': row.get('best_ir_pval'),
            # GC stats
            'gc_sig': row.get('total_gc_sig', 0),
            'gc_best_pval': row.get('best_gc_pval'),
            # IC stats
            'ic_sig': row.get('total_ic_sig', 0),
            'ic_best_pval': row.get('best_ic_pval'),
            # ONT breakdown
            'ont_samples_total': row.get('ont_samples_tested', 0),
            'ont_ir_sig': row.get('ont_ir_sig', 0),
            'ont_gc_sig': row.get('ont_gc_sig', 0),
            'ont_ic_sig': row.get('ont_ic_sig', 0),
            # SR breakdown
            'sr_samples_total': row.get('sr_samples_tested', 0),
            'sr_ir_sig': row.get('sr_ir_sig', 0),
            'sr_gc_sig': row.get('sr_gc_sig', 0),
            'sr_ic_sig': row.get('sr_ic_sig', 0),
            # scRNA
            'scrna_n_events': row.get('scrna_total_events', 0),
            'scrna_sig_events': row.get('scrna_sig_events', 0),
            'scrna_max_delta_psi': row.get('scrna_max_delta_psi'),
            'scrna_best_padj': row.get('scrna_min_padj'),
        }
        gene_info.append(info)

    # Add any genes from scRNA that might not be in validation report
    validation_genes = {g['gene'] for g in gene_info}
    for row in csv_data.get('scrna_full_psi_differential', []):
        gene = row.get('gene')
        if gene and gene not in validation_genes:
            validation_genes.add(gene)
            # Add minimal entry for scRNA-only genes
            gene_info.append({
                'gene': gene,
                'validation_tier': 'TIER3_SCRNA',
                'platforms': 'scRNA',
                'total_samples': 0,
                'ir_sig': 0, 'ir_best_pval': None,
                'gc_sig': 0, 'gc_best_pval': None,
                'ic_sig': 0, 'ic_best_pval': None,
                'ont_samples_total': 0, 'ont_ir_sig': 0, 'ont_gc_sig': 0, 'ont_ic_sig': 0,
                'sr_samples_total': 0, 'sr_ir_sig': 0, 'sr_gc_sig': 0, 'sr_ic_sig': 0,
                'scrna_n_events': 0, 'scrna_sig_events': 0,
                'scrna_max_delta_psi': None, 'scrna_best_padj': None,
            })

    return gene_info


def aggregate_sample_summary(csv_data):
    """Create sample-level summary table."""
    sample_info = defaultdict(lambda: {
        'sample_id': None,
        'platform': None,
        'n_spots': None,
        'genes_tested': 0,
        'genes_significant': 0,
        'best_gene': None,
        'best_pval': None
    })

    # ONT samples
    for row in csv_data.get('ont_splisosm_full_results', []):
        sample = row.get('sample_id')
        if not sample:
            continue
        sample_info[sample]['sample_id'] = sample
        sample_info[sample]['platform'] = 'ONT'
        sample_info[sample]['genes_tested'] += 1

        spots = row.get('n_spots_expressed')
        if spots:
            sample_info[sample]['n_spots'] = spots

        pval = row.get('pvalue_ir')
        if pval is not None and pval < 0.05:
            sample_info[sample]['genes_significant'] += 1
        if pval is not None:
            if sample_info[sample]['best_pval'] is None or pval < sample_info[sample]['best_pval']:
                sample_info[sample]['best_pval'] = pval
                sample_info[sample]['best_gene'] = row.get('gene')

    # SR samples
    for row in csv_data.get('sr_splisosm_full_results', []):
        sample = row.get('sample')
        if not sample:
            continue
        sample_info[sample]['sample_id'] = sample
        sample_info[sample]['platform'] = 'SR'
        sample_info[sample]['genes_tested'] += 1

        spots = row.get('n_spots')
        if spots:
            sample_info[sample]['n_spots'] = spots

        # Use pvalue_ir for the new comprehensive SR results
        pval = row.get('pvalue_ir') or row.get('pvalue')
        if pval is not None and pval < 0.05:
            sample_info[sample]['genes_significant'] += 1
        if pval is not None:
            if sample_info[sample]['best_pval'] is None or pval < sample_info[sample]['best_pval']:
                sample_info[sample]['best_pval'] = pval
                sample_info[sample]['best_gene'] = row.get('gene')

    return list(sample_info.values())


if __name__ == '__main__':
    main()

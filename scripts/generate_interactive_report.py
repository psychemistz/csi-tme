#!/usr/bin/env python3
"""
Generate interactive HTML dashboard for CSI-TME analysis results.

Reads dashboard_data.json and generates a single-file HTML with:
- DataTables for searchable/sortable tables
- Embedded base64 images
- Tab-based navigation
- Interactive filtering

Output: reports/interactive_dashboard.html
"""

import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
REPORTS_DIR = PROJECT_ROOT / "reports"
DATA_FILE = REPORTS_DIR / "dashboard_data.json"
OUTPUT_FILE = REPORTS_DIR / "interactive_dashboard.html"


def generate_html(data):
    """Generate the complete HTML document."""

    # Pre-compute JSON strings for embedding
    gene_summary_json = json.dumps(data.get('gene_summary', []))
    sample_summary_json = json.dumps(data.get('sample_summary', []))
    ont_results_json = json.dumps(data.get('csv_data', {}).get('ont_splisosm_full_results', []))
    sr_results_json = json.dumps(data.get('csv_data', {}).get('sr_splisosm_full_results', []))
    scrna_results_json = json.dumps(data.get('csv_data', {}).get('scrna_full_psi_differential', []))
    spacet_results_json = json.dumps(data.get('csv_data', {}).get('spacet_isoform_correlations', []))
    validation_json = json.dumps(data.get('csv_data', {}).get('final_validation_report', []))

    # Prepare image data for interactive viewer
    images = data.get('images', {})

    def prepare_image_data(img_dict, category):
        """Convert image dict to list format for JS."""
        result = []
        for name, img_data in img_dict.items():
            result.append({
                'name': name,
                'label': name.replace('_', ' ').replace('fig', 'Fig').title(),
                'data': f"data:{img_data.get('type', 'image/png')};base64,{img_data.get('data', '')}",
                'filename': img_data.get('filename', name)
            })
        return result

    main_figures_list = prepare_image_data(images.get('main_figures', {}), 'main')
    spacet_plots_list = prepare_image_data(images.get('spacet_plots', {}), 'spacet')
    spatial_plots_list = prepare_image_data(images.get('spatial_plots', {}), 'spatial')

    images_json = json.dumps({
        'main': main_figures_list,
        'spacet': spacet_plots_list,
        'spatial': spatial_plots_list
    })

    # Summary stats
    stats = data.get('summary_stats', {})
    ont_stats = stats.get('ont', {})
    sr_stats = stats.get('sr', {})
    val_stats = stats.get('validation', {})
    scrna_stats = stats.get('scrna', {})

    # Images
    images = data.get('images', {})
    main_figures = images.get('main_figures', {})
    spacet_plots = images.get('spacet_plots', {})
    spatial_plots = images.get('spatial_plots', {})

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CSI-TME Interactive Dashboard</title>

    <!-- DataTables CSS -->
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">

    <style>
        :root {{
            --primary-color: #3498db;
            --secondary-color: #2c3e50;
            --success-color: #27ae60;
            --warning-color: #f39c12;
            --danger-color: #e74c3c;
            --light-bg: #f5f5f5;
            --card-bg: #ffffff;
        }}

        * {{
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 0;
            background-color: var(--light-bg);
            color: #333;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}

        .header h1 {{
            margin: 0 0 10px 0;
            font-size: 2em;
        }}

        .header p {{
            margin: 0;
            opacity: 0.9;
        }}

        .nav-tabs {{
            display: flex;
            background-color: var(--secondary-color);
            padding: 0 20px;
            overflow-x: auto;
        }}

        .nav-tab {{
            padding: 15px 25px;
            color: rgba(255,255,255,0.7);
            cursor: pointer;
            border: none;
            background: none;
            font-size: 1em;
            white-space: nowrap;
            transition: all 0.3s;
        }}

        .nav-tab:hover {{
            color: white;
            background-color: rgba(255,255,255,0.1);
        }}

        .nav-tab.active {{
            color: white;
            background-color: var(--primary-color);
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}

        .tab-content {{
            display: none;
        }}

        .tab-content.active {{
            display: block;
        }}

        .card {{
            background-color: var(--card-bg);
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            padding: 20px;
            margin-bottom: 20px;
        }}

        .card h2 {{
            color: var(--secondary-color);
            border-left: 4px solid var(--primary-color);
            padding-left: 15px;
            margin-top: 0;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}

        .stat-box {{
            text-align: center;
            padding: 25px;
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            border-radius: 10px;
        }}

        .stat-number {{
            font-size: 2.5em;
            font-weight: bold;
            color: var(--primary-color);
        }}

        .stat-label {{
            color: #666;
            font-size: 0.9em;
            margin-top: 5px;
        }}

        .tier-badge {{
            display: inline-block;
            padding: 3px 10px;
            border-radius: 15px;
            font-size: 0.85em;
            font-weight: bold;
        }}

        .tier-1 {{ background-color: #27ae60; color: white; }}
        .tier-2 {{ background-color: #3498db; color: white; }}
        .tier-3 {{ background-color: #f39c12; color: white; }}
        .tier-4 {{ background-color: #95a5a6; color: white; }}

        .platform-badge {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 4px;
            font-size: 0.8em;
            font-weight: bold;
        }}

        .platform-ont {{ background-color: #9b59b6; color: white; }}
        .platform-sr {{ background-color: #16a085; color: white; }}

        table.dataTable {{
            width: 100% !important;
        }}

        .dataTables_wrapper {{
            padding: 10px 0;
        }}

        .filter-row {{
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            margin-bottom: 20px;
            align-items: center;
        }}

        .filter-group {{
            display: flex;
            flex-direction: column;
            gap: 5px;
        }}

        .filter-group label {{
            font-size: 0.85em;
            color: #666;
            font-weight: 500;
        }}

        .filter-group select, .filter-group input {{
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 5px;
            font-size: 0.9em;
        }}

        .image-gallery {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 20px;
        }}

        .gallery-item {{
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            cursor: pointer;
            transition: transform 0.3s;
        }}

        .gallery-item:hover {{
            transform: translateY(-5px);
        }}

        .gallery-item img {{
            width: 100%;
            height: 200px;
            object-fit: cover;
        }}

        .gallery-item .caption {{
            padding: 15px;
            font-size: 0.9em;
        }}

        .modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.9);
        }}

        .modal-content {{
            max-width: 90%;
            max-height: 90%;
            margin: auto;
            display: block;
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
        }}

        .modal-close {{
            position: absolute;
            top: 20px;
            right: 30px;
            color: white;
            font-size: 40px;
            cursor: pointer;
        }}

        .detail-panel {{
            background: #f8f9fa;
            border-radius: 10px;
            padding: 20px;
            margin-top: 20px;
            display: none;
        }}

        .detail-panel.active {{
            display: block;
        }}

        .detail-panel h3 {{
            margin-top: 0;
            color: var(--secondary-color);
        }}

        .key-finding {{
            background-color: #e8f6ff;
            border-left: 4px solid var(--primary-color);
            padding: 15px;
            margin: 20px 0;
        }}

        .download-section {{
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
        }}

        .download-btn {{
            display: inline-block;
            padding: 12px 25px;
            background-color: var(--primary-color);
            color: white;
            text-decoration: none;
            border-radius: 5px;
            font-weight: 500;
            transition: background-color 0.3s;
            border: none;
            cursor: pointer;
            font-size: 1em;
        }}

        .download-btn:hover {{
            background-color: #2980b9;
        }}

        .pval-significant {{
            color: var(--success-color);
            font-weight: bold;
        }}

        .pval-marginal {{
            color: var(--warning-color);
        }}

        @media (max-width: 768px) {{
            .nav-tabs {{
                flex-wrap: nowrap;
                overflow-x: auto;
            }}

            .nav-tab {{
                padding: 12px 15px;
                font-size: 0.9em;
            }}

            .stats-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>CSI-TME: Spatial Isoform Analysis Dashboard</h1>
        <p>Cytokine & Secreted protein spatial Isoform distribution in the Tumor MicroEnvironment</p>
    </div>

    <div class="nav-tabs">
        <button class="nav-tab active" data-tab="overview">Overview</button>
        <button class="nav-tab" data-tab="genes">Gene Results</button>
        <button class="nav-tab" data-tab="samples">Sample Results</button>
        <button class="nav-tab" data-tab="scrna">scRNA Validation</button>
        <button class="nav-tab" data-tab="visualizations">Visualizations</button>
        <button class="nav-tab" data-tab="downloads">Downloads</button>
    </div>

    <div class="container">
        <!-- OVERVIEW TAB -->
        <div id="overview" class="tab-content active">
            <div class="card">
                <h2>Analysis Summary</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-number">{ont_stats.get('n_samples', 0) + sr_stats.get('n_samples', 0)}</div>
                        <div class="stat-label">Total Visium Samples<br>(ONT: {ont_stats.get('n_samples', 0)}, SR: {sr_stats.get('n_samples', 0)})</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-number">{val_stats.get('n_genes_tested', 0)}</div>
                        <div class="stat-label">Genes Analyzed</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-number">{scrna_stats.get('genes_with_sig', 0)}</div>
                        <div class="stat-label">Genes Validated (scRNA)</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-number">{scrna_stats.get('n_significant', 0)}</div>
                        <div class="stat-label">Significant PSI Events</div>
                    </div>
                </div>
            </div>

            <div class="card">
                <h2>Validation Tier Summary</h2>
                <div class="key-finding">
                    <strong>Validation Tiers:</strong><br>
                    <span class="tier-badge tier-1">TIER1_FULL</span> Strong spatial + scRNA validation<br>
                    <span class="tier-badge tier-2">TIER2_SPATIAL</span> Reproducible spatial patterns only<br>
                    <span class="tier-badge tier-3">TIER3_SCRNA</span> scRNA validation but limited spatial<br>
                    <span class="tier-badge tier-4">TIER4_PARTIAL</span> Partial evidence in either modality
                </div>
                <div id="tier-chart-container">
                    <canvas id="tierChart" width="400" height="200"></canvas>
                </div>
            </div>

            <div class="card">
                <h2>Key Findings</h2>
                <div class="key-finding">
                    <strong>VCAN (Versican):</strong> Shows near-complete isoform switching (99% ΔPSI) between macrophages and tumor cells. The largest splicing effect in the dataset, validated both spatially (8/13 SR samples significant) and in scRNA-seq (p=1.4e-180).
                </div>
                <div class="key-finding">
                    <strong>APP (Amyloid Precursor Protein):</strong> Exon 7/8 inclusion varies dramatically between tumor states (MES: 82% vs NPC: 27%). Spatial HSIC-IR significant in 4/7 samples.
                </div>
                <div class="key-finding">
                    <strong>CLU (Clusterin):</strong> Most reproducible spatial isoform variation, significant in ALL 11 ONT samples (100%) and validated in scRNA-seq with 26 differential events.
                </div>
            </div>
        </div>

        <!-- GENES TAB -->
        <div id="genes" class="tab-content">
            <div class="card">
                <h2>Gene-Level Results</h2>
                <div class="filter-row">
                    <div class="filter-group">
                        <label>Validation Tier</label>
                        <select id="filter-tier">
                            <option value="">All Tiers</option>
                            <option value="TIER1">TIER1 (Full)</option>
                            <option value="TIER2">TIER2 (Spatial)</option>
                            <option value="TIER3">TIER3 (scRNA)</option>
                            <option value="TIER4">TIER4 (Partial)</option>
                        </select>
                    </div>
                    <div class="filter-group">
                        <label>Min Samples Significant</label>
                        <input type="number" id="filter-min-samples" value="0" min="0" max="13">
                    </div>
                </div>
                <div style="background:#f8f9fa; padding:10px; border-radius:8px; margin-bottom:15px; font-size:0.9em;">
                    <strong>Test Statistics:</strong>
                    <code>IR</code>=Isoform Ratio (primary),
                    <code>GC</code>=Gene Counts (expression),
                    <code>IC</code>=Isoform Counts (combined)
                </div>
                <table id="geneTable" class="display" style="width:100%">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Tier</th>
                            <th>IR Sig</th>
                            <th>IR p-val</th>
                            <th>GC Sig</th>
                            <th>GC p-val</th>
                            <th>IC Sig</th>
                            <th>IC p-val</th>
                            <th>scRNA ΔPSI</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
            </div>

            <div id="geneDetailPanel" class="detail-panel">
                <h3>Gene Details: <span id="selectedGeneName"></span></h3>
                <div id="geneDetailContent"></div>
            </div>
        </div>

        <!-- SAMPLES TAB -->
        <div id="samples" class="tab-content">
            <div class="card">
                <h2>Sample-Level Results</h2>
                <div class="filter-row">
                    <div class="filter-group">
                        <label>Platform</label>
                        <select id="filter-platform">
                            <option value="">All Platforms</option>
                            <option value="ONT">ONT (Long-read)</option>
                            <option value="SR">SR (Short-read)</option>
                        </select>
                    </div>
                </div>
                <table id="sampleTable" class="display" style="width:100%">
                    <thead>
                        <tr>
                            <th>Sample ID</th>
                            <th>Platform</th>
                            <th>Spots</th>
                            <th>Genes Tested</th>
                            <th>Genes Significant</th>
                            <th>Best Gene</th>
                            <th>Best p-value</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
            </div>

            <div id="sampleDetailPanel" class="detail-panel">
                <h3>Sample Details: <span id="selectedSampleName"></span></h3>
                <div id="sampleDetailContent"></div>
            </div>
        </div>

        <!-- SCRNA TAB -->
        <div id="scrna" class="tab-content">
            <div class="card">
                <h2>scRNA-seq Differential Splicing</h2>
                <p>32,304 cells from glioma scRNA-seq dataset. PSI values quantified by MARVEL.</p>

                <div class="event-type-legend" style="background:#f8f9fa; padding:15px; border-radius:8px; margin-bottom:20px;">
                    <strong>Event Types:</strong>
                    <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(200px, 1fr)); gap:10px; margin-top:10px;">
                        <div><code>SE</code> - Skipped Exon (exon included or excluded)</div>
                        <div><code>MXE</code> - Mutually Exclusive Exons (one or the other)</div>
                        <div><code>A5SS</code> - Alternative 5' Splice Site</div>
                        <div><code>A3SS</code> - Alternative 3' Splice Site</div>
                        <div><code>AFE</code> - Alternative First Exon</div>
                        <div><code>ALE</code> - Alternative Last Exon</div>
                    </div>
                </div>

                <div class="filter-row">
                    <div class="filter-group">
                        <label>Gene</label>
                        <select id="filter-scrna-gene">
                            <option value="">All Genes</option>
                        </select>
                    </div>
                    <div class="filter-group">
                        <label>Min |ΔPSI|</label>
                        <input type="number" id="filter-min-dpsi" value="0.1" min="0" max="1" step="0.05">
                    </div>
                    <div class="filter-group">
                        <label>Max padj</label>
                        <input type="number" id="filter-max-padj" value="0.05" min="0" max="1" step="0.01">
                    </div>
                </div>
                <table id="scrnaTable" class="display" style="width:100%">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Event Type</th>
                            <th>Exon/Region</th>
                            <th>State 1</th>
                            <th>State 2</th>
                            <th>PSI 1</th>
                            <th>PSI 2</th>
                            <th>ΔPSI</th>
                            <th>padj</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
            </div>
        </div>

        <!-- VISUALIZATIONS TAB -->
        <div id="visualizations" class="tab-content">
            <div class="card">
                <h2>Interactive Figure Viewer</h2>
                <p>Select a category and figure to view. Click the image to enlarge.</p>

                <div class="filter-row" style="margin-bottom:20px;">
                    <div class="filter-group">
                        <label>Category</label>
                        <select id="viz-category" onchange="updateVizOptions()">
                            <option value="main">Main Figures</option>
                            <option value="spacet">SpaCET Deconvolution</option>
                            <option value="spatial">Spatial Isoform Maps</option>
                        </select>
                    </div>
                    <div class="filter-group">
                        <label>Figure</label>
                        <select id="viz-figure" onchange="showSelectedFigure()">
                        </select>
                    </div>
                </div>

                <div id="viz-description" style="background:#e8f6ff; padding:15px; border-radius:8px; margin-bottom:20px;">
                    <strong>Main Figures:</strong> Key results including PSI boxplots and spatial isoform switching patterns.
                </div>

                <div id="viz-container" style="text-align:center; background:#f8f9fa; padding:20px; border-radius:10px; min-height:400px;">
                    <img id="viz-image" src="" alt="Select a figure" style="max-width:100%; max-height:600px; cursor:pointer; border-radius:5px; box-shadow:0 2px 10px rgba(0,0,0,0.1);" onclick="openModal(this.src)">
                </div>

                <div id="viz-caption" style="text-align:center; margin-top:15px; font-style:italic; color:#666;">
                </div>
            </div>
        </div>

        <!-- DOWNLOADS TAB -->
        <div id="downloads" class="tab-content">
            <div class="card">
                <h2>Download Data</h2>
                <p>Export analysis results and source data.</p>
                <div class="download-section">
                    <button class="download-btn" onclick="downloadCSV('gene_summary', 'gene_summary.csv')">Gene Summary (CSV)</button>
                    <button class="download-btn" onclick="downloadCSV('sample_summary', 'sample_summary.csv')">Sample Summary (CSV)</button>
                    <button class="download-btn" onclick="downloadCSV('scrna', 'scrna_differential.csv')">scRNA Results (CSV)</button>
                    <button class="download-btn" onclick="downloadCSV('ont', 'ont_results.csv')">ONT Results (CSV)</button>
                    <button class="download-btn" onclick="downloadCSV('sr', 'sr_results.csv')">SR Results (CSV)</button>
                </div>
            </div>

            <div class="card">
                <h2>Source Files</h2>
                <p>Original data files in reports/ directory:</p>
                <ul>
                    <li><code>ont_splisosm_full_results.csv</code> - ONT HSIC-IR results</li>
                    <li><code>pilot_splisosm_full_results.csv</code> - SR HSIC-IR results</li>
                    <li><code>scrna_full_psi_differential.csv</code> - All PSI comparisons</li>
                    <li><code>spacet_isoform_correlations.csv</code> - Cell type correlations</li>
                    <li><code>final_validation_report.csv</code> - Validation tier summary</li>
                </ul>
            </div>
        </div>
    </div>

    <!-- Image Modal -->
    <div id="imageModal" class="modal">
        <span class="modal-close" onclick="closeModal()">&times;</span>
        <img class="modal-content" id="modalImage">
    </div>

    <!-- JavaScript Libraries -->
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>

    <script>
        // Embedded data
        const geneSummary = {gene_summary_json};
        const sampleSummary = {sample_summary_json};
        const ontResults = {ont_results_json};
        const srResults = {sr_results_json};
        const scrnaResults = {scrna_results_json};
        const spacetResults = {spacet_results_json};
        const validationData = {validation_json};
        const imagesData = {images_json};

        // Category descriptions
        const categoryDescriptions = {{
            'main': '<strong>Main Figures:</strong> Key results including PSI boxplots by tumor state and spatial isoform switching patterns.',
            'spacet': '<strong>SpaCET Deconvolution:</strong> Cell type proportion correlations with VCAN isoform ratios. Shows tumor vs macrophage isoform preferences.',
            'spatial': '<strong>Spatial Isoform Maps:</strong> Spatial distribution of isoform ratios across tissue sections for each sample.'
        }};

        // Initialize visualization viewer
        function updateVizOptions() {{
            const category = document.getElementById('viz-category').value;
            const figSelect = document.getElementById('viz-figure');
            const images = imagesData[category] || [];

            // Clear and populate figure dropdown
            figSelect.innerHTML = '';
            images.forEach((img, idx) => {{
                const option = document.createElement('option');
                option.value = idx;
                option.textContent = img.label;
                figSelect.appendChild(option);
            }});

            // Update description
            document.getElementById('viz-description').innerHTML = categoryDescriptions[category] || '';

            // Show first figure
            showSelectedFigure();
        }}

        function showSelectedFigure() {{
            const category = document.getElementById('viz-category').value;
            const figIdx = document.getElementById('viz-figure').value;
            const images = imagesData[category] || [];

            if (images.length > 0 && figIdx < images.length) {{
                const img = images[figIdx];
                document.getElementById('viz-image').src = img.data;
                document.getElementById('viz-caption').textContent = img.label + ' (' + img.filename + ')';
            }} else {{
                document.getElementById('viz-image').src = '';
                document.getElementById('viz-caption').textContent = 'No figures available';
            }}
        }}

        // Initialize on page load
        document.addEventListener('DOMContentLoaded', function() {{
            updateVizOptions();
        }});

        // Tab switching
        document.querySelectorAll('.nav-tab').forEach(tab => {{
            tab.addEventListener('click', function() {{
                document.querySelectorAll('.nav-tab').forEach(t => t.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
                this.classList.add('active');
                document.getElementById(this.dataset.tab).classList.add('active');
            }});
        }});

        // Format p-value for display
        function formatPval(pval) {{
            if (pval === null || pval === undefined) return '-';
            if (pval < 0.001) return pval.toExponential(1);
            return pval.toFixed(4);
        }}

        // Format p-value with color
        function formatPvalHtml(pval) {{
            if (pval === null || pval === undefined) return '-';
            const formatted = formatPval(pval);
            if (pval < 0.001) return '<span class="pval-significant">' + formatted + '</span>';
            if (pval < 0.05) return '<span class="pval-marginal">' + formatted + '</span>';
            return formatted;
        }}

        // Get tier badge HTML
        function getTierBadge(tier) {{
            if (!tier) return '-';
            const tierNum = tier.includes('1') ? 1 : tier.includes('2') ? 2 : tier.includes('3') ? 3 : 4;
            return '<span class="tier-badge tier-' + tierNum + '">' + tier + '</span>';
        }}

        // Get platform badge
        function getPlatformBadge(platform) {{
            const cls = platform === 'ONT' ? 'platform-ont' : 'platform-sr';
            return '<span class="platform-badge ' + cls + '">' + platform + '</span>';
        }}

        // Initialize Gene Table
        let geneTable;
        $(document).ready(function() {{
            geneTable = $('#geneTable').DataTable({{
                data: geneSummary,
                columns: [
                    {{ data: 'gene', render: (data, type, row) => '<a href="#" onclick="showGeneDetail(\\'' + data + '\\'); return false;">' + data + '</a>' }},
                    {{ data: 'validation_tier', render: (data, type) => type === 'display' ? getTierBadge(data) : data }},
                    {{ data: null, render: (data, type, row) => row.ir_sig + '/' + row.total_samples }},
                    {{ data: 'ir_best_pval', render: (data, type) => type === 'display' ? formatPvalHtml(data) : data }},
                    {{ data: null, render: (data, type, row) => row.gc_sig + '/' + row.total_samples }},
                    {{ data: 'gc_best_pval', render: (data, type) => type === 'display' ? formatPvalHtml(data) : data }},
                    {{ data: null, render: (data, type, row) => row.ic_sig + '/' + row.total_samples }},
                    {{ data: 'ic_best_pval', render: (data, type) => type === 'display' ? formatPvalHtml(data) : data }},
                    {{ data: 'scrna_max_delta_psi', render: (data, type) => data ? (data * 100).toFixed(1) + '%' : '-' }}
                ],
                order: [[3, 'asc']],
                pageLength: 25
            }});

            // Gene table filters
            $('#filter-tier').on('change', function() {{
                const val = this.value;
                geneTable.column(1).search(val).draw();
            }});

            $('#filter-min-samples').on('change', function() {{
                geneTable.draw();
            }});

            $.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {{
                if (settings.nTable.id !== 'geneTable') return true;
                const minSamples = parseInt($('#filter-min-samples').val()) || 0;
                const irSig = parseInt(data[2].split('/')[0]) || 0;
                return irSig >= minSamples;
            }});
        }});

        // Initialize Sample Table
        let sampleTable;
        $(document).ready(function() {{
            sampleTable = $('#sampleTable').DataTable({{
                data: sampleSummary,
                columns: [
                    {{ data: 'sample_id', render: (data, type, row) => '<a href="#" onclick="showSampleDetail(\\'' + data + '\\'); return false;">' + data + '</a>' }},
                    {{ data: 'platform', render: (data, type) => type === 'display' ? getPlatformBadge(data) : data }},
                    {{ data: 'n_spots' }},
                    {{ data: 'genes_tested' }},
                    {{ data: 'genes_significant' }},
                    {{ data: 'best_gene' }},
                    {{ data: 'best_pval', render: (data, type) => type === 'display' ? formatPvalHtml(data) : data }}
                ],
                order: [[6, 'asc']],
                pageLength: 25
            }});

            $('#filter-platform').on('change', function() {{
                sampleTable.column(1).search(this.value).draw();
            }});
        }});

        // Initialize scRNA Table
        let scrnaTable;
        $(document).ready(function() {{
            // Populate gene filter
            const genes = [...new Set(scrnaResults.map(r => r.gene))].sort();
            genes.forEach(g => {{
                $('#filter-scrna-gene').append('<option value="' + g + '">' + g + '</option>');
            }});

            // Function to shorten event_id to readable exon coordinates
            function shortenEventId(eventId) {{
                if (!eventId) return '-';
                // Extract coordinates: chr5:83519349:83522309:+@chr5:83537007:83542268
                // Show as: 83519349-83522309
                const match = eventId.match(/:(\d+):(\d+):/);
                if (match) {{
                    return match[1] + '-' + match[2];
                }}
                return eventId.substring(0, 20) + '...';
            }}

            scrnaTable = $('#scrnaTable').DataTable({{
                data: scrnaResults,
                columns: [
                    {{ data: 'gene' }},
                    {{ data: 'event_type' }},
                    {{ data: 'event_id', render: (data, type) => type === 'display' ? '<span title=\"' + data + '\">' + shortenEventId(data) + '</span>' : data }},
                    {{ data: 'state1' }},
                    {{ data: 'state2' }},
                    {{ data: 'mean_psi_state1', render: (data) => data ? (data * 100).toFixed(1) + '%' : '-' }},
                    {{ data: 'mean_psi_state2', render: (data) => data ? (data * 100).toFixed(1) + '%' : '-' }},
                    {{ data: 'delta_psi', render: (data) => data ? (data * 100).toFixed(1) + '%' : '-' }},
                    {{ data: 'padj', render: (data, type) => type === 'display' ? formatPvalHtml(data) : data }}
                ],
                order: [[8, 'asc']],
                pageLength: 25
            }});

            // scRNA filters
            $('#filter-scrna-gene').on('change', function() {{
                scrnaTable.column(0).search(this.value).draw();
            }});

            $.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {{
                if (settings.nTable.id !== 'scrnaTable') return true;
                const minDpsi = parseFloat($('#filter-min-dpsi').val()) || 0;
                const maxPadj = parseFloat($('#filter-max-padj').val()) || 1;
                const dpsi = Math.abs(parseFloat(data[7])) / 100 || 0;  // ΔPSI is now column 7
                const padj = parseFloat(scrnaResults[dataIndex].padj) || 1;
                return dpsi >= minDpsi && padj <= maxPadj;
            }});

            $('#filter-min-dpsi, #filter-max-padj').on('change', function() {{
                scrnaTable.draw();
            }});
        }});

        // Tier distribution chart
        $(document).ready(function() {{
            const tierCounts = {{}};
            validationData.forEach(row => {{
                const tier = row.validation_tier || row.overall_validation || 'Unknown';
                tierCounts[tier] = (tierCounts[tier] || 0) + 1;
            }});

            const ctx = document.getElementById('tierChart');
            if (ctx) {{
                new Chart(ctx, {{
                    type: 'bar',
                    data: {{
                        labels: Object.keys(tierCounts),
                        datasets: [{{
                            label: 'Number of Genes',
                            data: Object.values(tierCounts),
                            backgroundColor: ['#27ae60', '#3498db', '#f39c12', '#95a5a6'],
                            borderWidth: 0
                        }}]
                    }},
                    options: {{
                        responsive: true,
                        plugins: {{
                            legend: {{ display: false }}
                        }},
                        scales: {{
                            y: {{ beginAtZero: true }}
                        }}
                    }}
                }});
            }}
        }});

        // Gene detail panel
        function showGeneDetail(gene) {{
            document.getElementById('geneDetailPanel').classList.add('active');
            document.getElementById('selectedGeneName').textContent = gene;

            // Find gene data
            const geneData = geneSummary.find(g => g.gene === gene);
            const ontGeneResults = ontResults.filter(r => r.gene === gene);
            const srGeneResults = srResults.filter(r => r.gene === gene);
            const scrnaGeneResults = scrnaResults.filter(r => r.gene === gene).slice(0, 10);

            let html = '<h4>Summary</h4>';
            html += '<p><strong>Validation Tier:</strong> ' + getTierBadge(geneData?.validation_tier) + '</p>';
            html += '<p><strong>Platforms:</strong> ' + (geneData?.platforms || '-') + '</p>';

            // Statistics summary table
            html += '<h4>Spatial Statistics (ONT + SR combined)</h4>';
            html += '<table class="display compact"><thead><tr><th>Test</th><th>Significant</th><th>Best p-value</th></tr></thead><tbody>';
            html += '<tr><td><strong>HSIC-IR</strong> (Isoform Ratio)</td><td>' + (geneData?.ir_sig || 0) + '/' + (geneData?.total_samples || 0) + '</td><td>' + formatPvalHtml(geneData?.ir_best_pval) + '</td></tr>';
            html += '<tr><td><strong>HSIC-GC</strong> (Gene Counts)</td><td>' + (geneData?.gc_sig || 0) + '/' + (geneData?.total_samples || 0) + '</td><td>' + formatPvalHtml(geneData?.gc_best_pval) + '</td></tr>';
            html += '<tr><td><strong>HSIC-IC</strong> (Isoform Counts)</td><td>' + (geneData?.ic_sig || 0) + '/' + (geneData?.total_samples || 0) + '</td><td>' + formatPvalHtml(geneData?.ic_best_pval) + '</td></tr>';
            html += '</tbody></table>';

            // Platform breakdown
            if (geneData?.ont_samples_total > 0 || geneData?.sr_samples_total > 0) {{
                html += '<h4>Platform Breakdown</h4>';
                html += '<table class="display compact"><thead><tr><th>Platform</th><th>Samples</th><th>IR Sig</th><th>GC Sig</th><th>IC Sig</th></tr></thead><tbody>';
                if (geneData?.ont_samples_total > 0) {{
                    html += '<tr><td>ONT</td><td>' + geneData.ont_samples_total + '</td><td>' + (geneData.ont_ir_sig || 0) + '</td><td>' + (geneData.ont_gc_sig || 0) + '</td><td>' + (geneData.ont_ic_sig || 0) + '</td></tr>';
                }}
                if (geneData?.sr_samples_total > 0) {{
                    html += '<tr><td>SR</td><td>' + geneData.sr_samples_total + '</td><td>' + (geneData.sr_ir_sig || 0) + '</td><td>' + (geneData.sr_gc_sig || 0) + '</td><td>' + (geneData.sr_ic_sig || 0) + '</td></tr>';
                }}
                html += '</tbody></table>';
            }}

            // Per-sample ONT results
            if (ontGeneResults.length > 0) {{
                html += '<h4>ONT Per-Sample Results</h4><table class="display compact"><thead><tr><th>Sample</th><th>Isoforms</th><th>IR p</th><th>GC p</th><th>IC p</th></tr></thead><tbody>';
                ontGeneResults.sort((a,b) => (a.pvalue_ir || 1) - (b.pvalue_ir || 1));
                ontGeneResults.forEach(r => {{
                    html += '<tr><td>' + r.sample_id + '</td><td>' + r.n_isoforms + '</td><td>' + formatPvalHtml(r.pvalue_ir) + '</td><td>' + formatPvalHtml(r.pvalue_gc) + '</td><td>' + formatPvalHtml(r.pvalue_ic) + '</td></tr>';
                }});
                html += '</tbody></table>';
            }}

            // Per-sample SR results
            if (srGeneResults.length > 0) {{
                html += '<h4>SR Per-Sample Results</h4><table class="display compact"><thead><tr><th>Sample</th><th>Isoforms</th><th>IR p</th><th>GC p</th><th>IC p</th></tr></thead><tbody>';
                srGeneResults.sort((a,b) => (a.pvalue_ir || 1) - (b.pvalue_ir || 1));
                srGeneResults.forEach(r => {{
                    html += '<tr><td>' + r.sample + '</td><td>' + r.n_isoforms + '</td><td>' + formatPvalHtml(r.pvalue_ir) + '</td><td>' + formatPvalHtml(r.pvalue_gc) + '</td><td>' + formatPvalHtml(r.pvalue_ic) + '</td></tr>';
                }});
                html += '</tbody></table>';
            }}

            // scRNA results
            if (scrnaGeneResults.length > 0) {{
                html += '<h4>Top scRNA-seq Events</h4><table class="display compact"><thead><tr><th>States</th><th>ΔPSI</th><th>padj</th></tr></thead><tbody>';
                scrnaGeneResults.forEach(r => {{
                    html += '<tr><td>' + r.state1 + ' vs ' + r.state2 + '</td><td>' + (r.delta_psi*100).toFixed(1) + '%</td><td>' + formatPvalHtml(r.padj) + '</td></tr>';
                }});
                html += '</tbody></table>';
            }}

            document.getElementById('geneDetailContent').innerHTML = html;
        }}

        // Sample detail panel
        function showSampleDetail(sample) {{
            document.getElementById('sampleDetailPanel').classList.add('active');
            document.getElementById('selectedSampleName').textContent = sample;

            const sampleData = sampleSummary.find(s => s.sample_id === sample);
            const isONT = sampleData?.platform === 'ONT';
            const results = isONT ? ontResults.filter(r => r.sample_id === sample) : srResults.filter(r => r.sample === sample);
            const spacetData = spacetResults.filter(r => r.sample === sample);

            let html = '<h4>Summary</h4>';
            html += '<p><strong>Platform:</strong> ' + getPlatformBadge(sampleData?.platform || '-') + '</p>';
            html += '<p><strong>Spots:</strong> ' + (sampleData?.n_spots || '-') + '</p>';
            html += '<p><strong>Genes Significant:</strong> ' + (sampleData?.genes_significant || 0) + '/' + (sampleData?.genes_tested || 0) + '</p>';

            if (results.length > 0) {{
                html += '<h4>Gene Results</h4><table class="display compact"><thead><tr><th>Gene</th><th>Isoforms</th><th>HSIC-IR p</th><th>HSIC-GC p</th></tr></thead><tbody>';
                results.sort((a, b) => (a.pvalue_ir || a.pvalue || 1) - (b.pvalue_ir || b.pvalue || 1));
                results.forEach(r => {{
                    const pval_ir = r.pvalue_ir || r.pvalue;
                    const pval_gc = r.pvalue_gc;
                    const niso = r.n_isoforms;
                    html += '<tr><td>' + r.gene + '</td><td>' + niso + '</td><td>' + formatPvalHtml(pval_ir) + '</td><td>' + formatPvalHtml(pval_gc) + '</td></tr>';
                }});
                html += '</tbody></table>';
            }}

            if (spacetData.length > 0) {{
                html += '<h4>SpaCET Cell Type Correlations</h4><table class="display compact"><thead><tr><th>Cell Type</th><th>Correlation</th><th>p-value</th></tr></thead><tbody>';
                spacetData.sort((a, b) => Math.abs(b.correlation) - Math.abs(a.correlation));
                spacetData.forEach(r => {{
                    html += '<tr><td>' + r.cell_type + '</td><td>' + (r.correlation?.toFixed(3) || '-') + '</td><td>' + formatPvalHtml(r.pvalue) + '</td></tr>';
                }});
                html += '</tbody></table>';
            }}

            document.getElementById('sampleDetailContent').innerHTML = html;
        }}

        // Image modal
        function openModal(imgSrc) {{
            document.getElementById('imageModal').style.display = 'block';
            document.getElementById('modalImage').src = imgSrc;
        }}

        function closeModal() {{
            document.getElementById('imageModal').style.display = 'none';
        }}

        document.getElementById('imageModal').addEventListener('click', function(e) {{
            if (e.target === this) closeModal();
        }});

        document.addEventListener('keydown', function(e) {{
            if (e.key === 'Escape') closeModal();
        }});

        // CSV download
        function downloadCSV(type, filename) {{
            let data;
            switch(type) {{
                case 'gene_summary': data = geneSummary; break;
                case 'sample_summary': data = sampleSummary; break;
                case 'scrna': data = scrnaResults; break;
                case 'ont': data = ontResults; break;
                case 'sr': data = srResults; break;
                default: return;
            }}

            if (!data || data.length === 0) {{
                alert('No data available');
                return;
            }}

            const headers = Object.keys(data[0]);
            let csv = headers.join(',') + '\\n';
            data.forEach(row => {{
                csv += headers.map(h => {{
                    let val = row[h];
                    if (val === null || val === undefined) val = '';
                    if (typeof val === 'string' && val.includes(',')) val = '"' + val + '"';
                    return val;
                }}).join(',') + '\\n';
            }});

            const blob = new Blob([csv], {{ type: 'text/csv' }});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            a.click();
            URL.revokeObjectURL(url);
        }}
    </script>
</body>
</html>
'''
    return html


def generate_image_gallery(images_dict, prefix):
    """Generate HTML for image gallery items."""
    if not images_dict:
        return '<p>No images available.</p>'

    items = []
    for name, img_data in images_dict.items():
        b64 = img_data.get('data', '')
        filename = img_data.get('filename', name)
        img_type = img_data.get('type', 'image/png')

        # Create a friendly caption
        caption = filename.replace('.png', '').replace('_', ' ').title()

        items.append(f'''
            <div class="gallery-item" onclick="openModal('data:{img_type};base64,{b64}')">
                <img src="data:{img_type};base64,{b64}" alt="{caption}">
                <div class="caption">{caption}</div>
            </div>
        ''')

    return '\n'.join(items)


def main():
    print("Generating interactive HTML dashboard...")

    # Load prepared data
    if not DATA_FILE.exists():
        print(f"Error: {DATA_FILE} not found. Run prepare_interactive_data.py first.")
        return

    print(f"Loading data from {DATA_FILE}...")
    with open(DATA_FILE, 'r') as f:
        data = json.load(f)

    print("Generating HTML...")
    html = generate_html(data)

    print(f"Writing output to {OUTPUT_FILE}...")
    with open(OUTPUT_FILE, 'w') as f:
        f.write(html)

    file_size_mb = OUTPUT_FILE.stat().st_size / (1024 * 1024)
    print(f"Done! Output size: {file_size_mb:.2f} MB")

    if file_size_mb > 50:
        print("Warning: File size exceeds 50 MB. Consider reducing image quality.")


if __name__ == '__main__':
    main()

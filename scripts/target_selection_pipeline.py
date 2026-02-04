#!/usr/bin/env python3
"""
Systematic Target Selection Pipeline for Cytokine/Secreted Protein Isoform Analysis

This script implements the TARGET_GUIDE.md criteria to systematically evaluate
all CytoSig and SecAct genes for suitability in Visium 3' SR isoform analysis.

Criteria evaluated:
1. Transcript 3' architecture (Ensembl REST API)
2. APA site documentation (PolyASite 2.0)
3. Data availability (Zenodo isoform data)
4. Literature evidence (PubMed query templates)

Output: Tiered target list with confidence scores

Usage:
    python scripts/target_selection_pipeline.py --download-polyasite
    python scripts/target_selection_pipeline.py --run-full
    python scripts/target_selection_pipeline.py --quick-check HMGB1 SPP1 VEGFA
"""

import argparse
import gzip
import json
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.request import urlretrieve
from urllib.error import URLError

import numpy as np
import pandas as pd
import requests

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))


# =============================================================================
# Configuration
# =============================================================================

POLYASITE_BED_URL = "https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz"
ENSEMBL_REST_URL = "https://rest.ensembl.org"
DATA_DIR = Path("data/target_selection")
ZENODO_ISO_DIR = Path("data/zenodo_16905935/human_glioma_sr")

# Signature data paths
CYTOSIG_PATH = Path("/vf/users/parks34/projects/1ridgesig/SecActpy/secactpy/data/CytoSig.tsv.gz")
SECACT_PATH = Path("/vf/users/parks34/projects/1ridgesig/SecActpy/secactpy/data/SecAct.tsv.gz")
GENE_MAPPING_PATH = Path("/data/parks34/projects/2secactpy/cytoatlas-api/static/data/signature_gene_mapping.json")


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class GeneAssessment:
    """Assessment results for a single gene."""
    gene_symbol: str
    source: str  # 'cytosig' or 'secact'

    # Step 1: Ensembl transcript architecture
    n_transcripts: int = 0
    has_3prime_differences: bool = False
    last_exon_variants: int = 0
    utr3_length_range: Tuple[int, int] = (0, 0)

    # Step 2: PolyASite APA documentation
    n_apa_clusters: int = 0
    apa_cluster_positions: List[int] = field(default_factory=list)
    has_brain_expression: bool = False
    has_immune_expression: bool = False

    # Step 3: Data availability
    in_zenodo_data: bool = False
    n_peaks_detected: int = 0
    total_umis: int = 0
    peak_positions_match_apa: bool = False

    # Step 4: Literature evidence
    pubmed_apa_hits: int = 0
    pubmed_isoform_hits: int = 0
    known_functional_isoforms: bool = False
    isoform_type: str = ""  # '3prime', '5prime', 'internal', 'unknown'

    # Known issues from TARGET_GUIDE
    known_non_3prime_isoforms: bool = False
    exclusion_reason: str = ""

    # Final score
    score: int = 0
    tier: str = ""  # 'A', 'B', 'C', 'Exclude'

    def calculate_score(self):
        """Calculate confidence score (0-4 points)."""
        self.score = 0

        # +1: Has 3' transcript differences
        if self.has_3prime_differences or self.last_exon_variants > 1:
            self.score += 1

        # +1: Has documented APA sites
        if self.n_apa_clusters >= 2:
            self.score += 1

        # +1: Detectable in our data
        if self.in_zenodo_data and self.n_peaks_detected >= 2 and self.total_umis >= 100:
            self.score += 1

        # +1: Has literature support
        if self.pubmed_apa_hits > 0 or self.known_functional_isoforms:
            self.score += 1

        # Penalties
        if self.known_non_3prime_isoforms:
            self.score -= 1  # Known isoforms are not 3'-detectable

        if self.exclusion_reason:
            self.score = 0

        # Tier assignment
        if self.exclusion_reason:
            self.tier = "Exclude"
        elif self.score >= 4:
            self.tier = "A"
        elif self.score >= 3:
            self.tier = "B"
        elif self.score >= 2:
            self.tier = "C"
        else:
            self.tier = "Exclude"

        return self.score


# =============================================================================
# Step 1: Ensembl Transcript Architecture
# =============================================================================

def query_ensembl_transcripts(gene_symbol: str, max_retries: int = 3) -> Dict:
    """Query Ensembl REST API for transcript information."""

    # First, get gene ID from symbol
    url = f"{ENSEMBL_REST_URL}/xrefs/symbol/homo_sapiens/{gene_symbol}"
    headers = {"Content-Type": "application/json"}

    for attempt in range(max_retries):
        try:
            response = requests.get(url, headers=headers, timeout=10)
            if response.status_code == 429:  # Rate limited
                time.sleep(1)
                continue
            response.raise_for_status()
            data = response.json()
            break
        except (requests.RequestException, json.JSONDecodeError):
            if attempt == max_retries - 1:
                return {}
            time.sleep(0.5)
    else:
        return {}

    # Find Ensembl gene ID
    gene_id = None
    for entry in data:
        if entry.get('type') == 'gene':
            gene_id = entry.get('id')
            break

    if not gene_id:
        return {}

    # Get transcript information
    url = f"{ENSEMBL_REST_URL}/lookup/id/{gene_id}?expand=1"

    for attempt in range(max_retries):
        try:
            response = requests.get(url, headers=headers, timeout=10)
            if response.status_code == 429:
                time.sleep(1)
                continue
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, json.JSONDecodeError):
            if attempt == max_retries - 1:
                return {}
            time.sleep(0.5)

    return {}


def assess_3prime_architecture(gene_data: Dict) -> Tuple[int, bool, int, Tuple[int, int]]:
    """Assess 3' architecture from Ensembl transcript data."""

    transcripts = gene_data.get('Transcript', [])
    n_transcripts = len(transcripts)

    if n_transcripts == 0:
        return 0, False, 0, (0, 0)

    # Collect 3' end positions and UTR lengths
    end_positions = set()
    utr3_lengths = []

    for tx in transcripts:
        # Get transcript end (3' for + strand, 5' for - strand)
        strand = tx.get('strand', 1)
        if strand == 1:
            end_pos = tx.get('end', 0)
        else:
            end_pos = tx.get('start', 0)
        end_positions.add(end_pos)

        # Try to get 3' UTR info from exons
        exons = tx.get('Exon', [])
        if exons:
            # Approximate 3' UTR from last exon
            last_exon = exons[-1] if strand == 1 else exons[0]
            utr3_lengths.append(last_exon.get('end', 0) - last_exon.get('start', 0))

    # Assess 3' differences
    last_exon_variants = len(end_positions)
    has_3prime_diff = last_exon_variants > 1

    utr3_range = (min(utr3_lengths), max(utr3_lengths)) if utr3_lengths else (0, 0)

    return n_transcripts, has_3prime_diff, last_exon_variants, utr3_range


# =============================================================================
# Step 2: PolyASite APA Database
# =============================================================================

def download_polyasite(output_dir: Path) -> Path:
    """Download PolyASite human APA atlas."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "polyasite_human_hg38.bed.gz"

    if output_file.exists():
        print(f"PolyASite data already exists: {output_file}")
        return output_file

    print(f"Downloading PolyASite data from {POLYASITE_BED_URL}...")
    try:
        urlretrieve(POLYASITE_BED_URL, output_file)
        print(f"Downloaded to {output_file}")
    except URLError as e:
        print(f"Failed to download PolyASite data: {e}")
        return None

    return output_file


def load_polyasite_data(bed_file: Path) -> pd.DataFrame:
    """Load PolyASite BED file into DataFrame."""
    if not bed_file or not bed_file.exists():
        return pd.DataFrame()

    # PolyASite BED format columns
    columns = [
        'chrom', 'start', 'end', 'cluster_id', 'score', 'strand',
        'rep_pos', 'n_protocols', 'avg_tpm', 'gene_symbol', 'gene_id',
        'pas_signal', 'annotation'
    ]

    df = pd.read_csv(bed_file, sep='\t', names=columns, comment='#')
    return df


def query_polyasite(gene_symbol: str, polyasite_df: pd.DataFrame) -> Tuple[int, List[int], bool, bool]:
    """Query PolyASite data for a gene."""
    if polyasite_df.empty:
        return 0, [], False, False

    gene_data = polyasite_df[polyasite_df['gene_symbol'] == gene_symbol]

    n_clusters = len(gene_data)
    positions = gene_data['rep_pos'].tolist() if n_clusters > 0 else []

    # Check tissue expression (would need sample metadata for full check)
    # For now, assume present if detected
    has_brain = n_clusters > 0
    has_immune = n_clusters > 0

    return n_clusters, positions, has_brain, has_immune


# =============================================================================
# Step 3: Zenodo Data Availability
# =============================================================================

def check_zenodo_availability(gene_symbol: str, zenodo_dir: Path) -> Tuple[bool, int, int, List[int]]:
    """Check if gene has isoform data in Zenodo glioma samples."""
    import anndata as ad

    sample_dirs = list(zenodo_dir.glob("*"))

    total_peaks = 0
    total_umis = 0
    all_positions = []
    found = False

    for sample_dir in sample_dirs:
        iso_file = sample_dir / "iso.quant.h5ad"
        if not iso_file.exists():
            continue

        try:
            adata = ad.read_h5ad(iso_file)
            gene_mask = adata.var['gene_symbol'] == gene_symbol

            if gene_mask.sum() > 0:
                found = True
                total_peaks = max(total_peaks, gene_mask.sum())

                # Get UMI counts
                counts = adata[:, gene_mask].layers['counts']
                if hasattr(counts, 'toarray'):
                    counts = counts.toarray()
                total_umis = max(total_umis, int(counts.sum()))

                # Get positions
                gene_data = adata.var[gene_mask]
                if 'Fit.start' in gene_data.columns:
                    positions = gene_data['Fit.start'].tolist()
                    all_positions.extend(positions)

            # Only check first available sample for efficiency
            if found:
                break

        except Exception as e:
            continue

    return found, total_peaks, total_umis, list(set(all_positions))


# =============================================================================
# Step 4: Known Issues from TARGET_GUIDE
# =============================================================================

# Genes with known NON-3' isoform differences
KNOWN_NON_3PRIME_ISOFORMS = {
    'IL15': "5' signal peptide difference",
    'IL4': "Internal exon 2 skipping (IL-4delta2)",
    'VCAN': "Internal GAG domain exons (V0/V1/V2/V3)",
    'APP': "Internal exon 7/8 skipping (APP695/751/770)",
    'CLU': "5' TSS difference (nuclear vs secreted)",
    'SPP1': "Internal exon 4/5 skipping (SPP1a/b/c)",
}

# Genes regulated by proteolysis, not splicing
PROTEOLYSIS_REGULATED = {
    'IL1B': "Caspase-1 cleavage",
    'IL18': "Caspase-1 cleavage",
    'TNF': "ADAM17 shedding (membrane vs soluble)",
}

# Tier 1 targets from TARGET_GUIDE (documented 3' differences)
TIER1_TARGETS = {
    'VEGFA': "Terminal exon 8a vs 8b",
    'IL33': "Terminal exon differences affecting 3' pattern",
    'CXCL12': "C-terminal isoforms (alpha/beta/gamma)",
    'TNFSF10': "TRAIL-short terminal difference",
}

# APA-documented targets
APA_DOCUMENTED = {
    'CCL2': "3' UTR APA variation",
    'IL6': "APA in inflammation",
    'CXCL8': "APA during inflammation",
    'IL10': "APA in macrophages",
}


def apply_known_issues(assessment: GeneAssessment) -> GeneAssessment:
    """Apply known issues from TARGET_GUIDE to assessment."""
    gene = assessment.gene_symbol

    # Check for known non-3' isoforms
    if gene in KNOWN_NON_3PRIME_ISOFORMS:
        assessment.known_non_3prime_isoforms = True
        assessment.isoform_type = "non_3prime"
        assessment.exclusion_reason = f"Known isoforms are NOT 3'-detectable: {KNOWN_NON_3PRIME_ISOFORMS[gene]}"

    # Check for proteolysis regulation
    if gene in PROTEOLYSIS_REGULATED:
        assessment.exclusion_reason = f"Regulated by proteolysis: {PROTEOLYSIS_REGULATED[gene]}"

    # Check for Tier 1 targets
    if gene in TIER1_TARGETS:
        assessment.known_functional_isoforms = True
        assessment.isoform_type = "3prime"
        assessment.has_3prime_differences = True

    # Check for APA-documented targets
    if gene in APA_DOCUMENTED:
        assessment.pubmed_apa_hits = 1  # Documented
        assessment.isoform_type = "3prime"

    return assessment


# =============================================================================
# Main Assessment Pipeline
# =============================================================================

def assess_gene(gene_symbol: str, source: str, polyasite_df: pd.DataFrame,
                check_ensembl: bool = True, check_zenodo: bool = True) -> GeneAssessment:
    """Run full assessment pipeline for a single gene."""

    assessment = GeneAssessment(gene_symbol=gene_symbol, source=source)

    # Step 1: Ensembl transcript architecture
    if check_ensembl:
        try:
            gene_data = query_ensembl_transcripts(gene_symbol)
            if gene_data:
                n_tx, has_3p, n_last_exon, utr3_range = assess_3prime_architecture(gene_data)
                assessment.n_transcripts = n_tx
                assessment.has_3prime_differences = has_3p
                assessment.last_exon_variants = n_last_exon
                assessment.utr3_length_range = utr3_range
        except Exception as e:
            pass

    # Step 2: PolyASite APA
    n_apa, apa_pos, has_brain, has_immune = query_polyasite(gene_symbol, polyasite_df)
    assessment.n_apa_clusters = n_apa
    assessment.apa_cluster_positions = apa_pos
    assessment.has_brain_expression = has_brain
    assessment.has_immune_expression = has_immune

    # Step 3: Zenodo data availability
    if check_zenodo and ZENODO_ISO_DIR.exists():
        try:
            found, n_peaks, umis, positions = check_zenodo_availability(gene_symbol, ZENODO_ISO_DIR)
            assessment.in_zenodo_data = found
            assessment.n_peaks_detected = n_peaks
            assessment.total_umis = umis

            # Check if peaks match APA sites (within 1kb)
            if positions and apa_pos:
                for pos in positions:
                    for apa in apa_pos:
                        if abs(pos - apa) < 1000:
                            assessment.peak_positions_match_apa = True
                            break
        except Exception as e:
            pass

    # Step 4: Apply known issues
    assessment = apply_known_issues(assessment)

    # Calculate final score
    assessment.calculate_score()

    return assessment


def load_gene_lists() -> Tuple[List[str], List[str]]:
    """Load CytoSig and SecAct gene lists."""

    cytosig_genes = []
    secact_genes = []

    # Load CytoSig
    if CYTOSIG_PATH.exists():
        with gzip.open(CYTOSIG_PATH, 'rt') as f:
            df = pd.read_csv(f, sep='\t', index_col=0, nrows=0)
            cytosig_genes = list(df.columns)

    # Load gene mapping for CytoSig -> HGNC
    if GENE_MAPPING_PATH.exists():
        with open(GENE_MAPPING_PATH) as f:
            mapping = json.load(f)
            # Convert signature names to HGNC symbols
            cytosig_genes = [mapping['signature_to_hgnc'].get(g, g) for g in cytosig_genes]

    # Load SecAct
    if SECACT_PATH.exists():
        with gzip.open(SECACT_PATH, 'rt') as f:
            df = pd.read_csv(f, sep='\t', index_col=0, nrows=0)
            secact_genes = list(df.columns)

    return cytosig_genes, secact_genes


def run_full_pipeline(output_dir: Path, skip_ensembl: bool = False) -> pd.DataFrame:
    """Run full target selection pipeline."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load gene lists
    print("Loading gene lists...")
    cytosig_genes, secact_genes = load_gene_lists()
    print(f"  CytoSig: {len(cytosig_genes)} genes")
    print(f"  SecAct: {len(secact_genes)} genes")

    # Download/load PolyASite
    print("\nLoading PolyASite data...")
    polyasite_file = download_polyasite(output_dir)
    polyasite_df = load_polyasite_data(polyasite_file) if polyasite_file else pd.DataFrame()
    print(f"  Loaded {len(polyasite_df)} APA clusters")

    # Assess all genes
    results = []

    print("\nAssessing CytoSig genes...")
    for i, gene in enumerate(cytosig_genes):
        print(f"  [{i+1}/{len(cytosig_genes)}] {gene}", end='\r')
        assessment = assess_gene(gene, 'cytosig', polyasite_df,
                                 check_ensembl=not skip_ensembl)
        results.append(assessment)
        if not skip_ensembl:
            time.sleep(0.1)  # Rate limiting for Ensembl API

    print(f"\nAssessing SecAct genes (this may take a while)...")
    for i, gene in enumerate(secact_genes):
        if (i + 1) % 50 == 0:
            print(f"  [{i+1}/{len(secact_genes)}] {gene}")
        assessment = assess_gene(gene, 'secact', polyasite_df,
                                 check_ensembl=not skip_ensembl,
                                 check_zenodo=True)
        results.append(assessment)
        if not skip_ensembl:
            time.sleep(0.05)  # Rate limiting

    # Convert to DataFrame
    df = pd.DataFrame([vars(a) for a in results])

    # Save results
    output_file = output_dir / "target_assessment_results.csv"
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    tier_counts = df['tier'].value_counts()
    for tier in ['A', 'B', 'C', 'Exclude']:
        count = tier_counts.get(tier, 0)
        print(f"  Tier {tier}: {count} genes")

    # Show Tier A genes
    tier_a = df[df['tier'] == 'A']
    if len(tier_a) > 0:
        print(f"\nTier A targets ({len(tier_a)}):")
        for _, row in tier_a.iterrows():
            print(f"  {row['gene_symbol']}: score={row['score']}, "
                  f"APA={row['n_apa_clusters']}, peaks={row['n_peaks_detected']}")

    return df


def quick_check(genes: List[str]) -> pd.DataFrame:
    """Quick check for specific genes without full pipeline."""

    print("Loading PolyASite data...")
    polyasite_file = DATA_DIR / "polyasite_human_hg38.bed.gz"
    polyasite_df = load_polyasite_data(polyasite_file) if polyasite_file.exists() else pd.DataFrame()

    results = []
    for gene in genes:
        print(f"\nAssessing {gene}...")
        assessment = assess_gene(gene, 'query', polyasite_df,
                                 check_ensembl=True, check_zenodo=True)
        results.append(assessment)

        print(f"  Tier: {assessment.tier} (score: {assessment.score})")
        print(f"  Transcripts: {assessment.n_transcripts}")
        print(f"  3' differences: {assessment.has_3prime_differences}")
        print(f"  APA clusters: {assessment.n_apa_clusters}")
        print(f"  In Zenodo: {assessment.in_zenodo_data} ({assessment.n_peaks_detected} peaks)")
        if assessment.exclusion_reason:
            print(f"  ⚠️ Exclusion: {assessment.exclusion_reason}")

    return pd.DataFrame([vars(a) for a in results])


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Target selection pipeline")
    parser.add_argument("--download-polyasite", action="store_true",
                        help="Download PolyASite data only")
    parser.add_argument("--run-full", action="store_true",
                        help="Run full pipeline on all genes")
    parser.add_argument("--skip-ensembl", action="store_true",
                        help="Skip Ensembl API queries (faster)")
    parser.add_argument("--quick-check", nargs="+", metavar="GENE",
                        help="Quick check specific genes")
    parser.add_argument("--output-dir", type=Path, default=DATA_DIR,
                        help="Output directory")

    args = parser.parse_args()

    if args.download_polyasite:
        download_polyasite(args.output_dir)
    elif args.run_full:
        run_full_pipeline(args.output_dir, skip_ensembl=args.skip_ensembl)
    elif args.quick_check:
        quick_check(args.quick_check)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

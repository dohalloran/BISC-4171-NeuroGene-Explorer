#!/usr/bin/env python3
"""
BISC 4171 - NeuroGene Explorer
00_download_gse113754.py

Downloads selected public files for GEO Series GSE113754.

Default: downloads a processed counts table (no FASTQ processing).
Optionally downloads the Series Matrix file for metadata reference.

Usage:
  python src/00_download_gse113754.py --outdir data
  python src/00_download_gse113754.py --outdir data --series-matrix
"""

from __future__ import annotations
import argparse
from pathlib import Path
from urllib.request import urlretrieve

GSE = "GSE113754"
BASE = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/{GSE}"

# NOTE: If NCBI file names change, students can still download from the GEO page.
DEFAULT_SUPPL = "GSE113754_RNASeq_S3_WT_SD5_HC5_Counts.txt.gz"
SERIES_MATRIX = f"{GSE}_series_matrix.txt.gz"

def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading:\n  {url}\n-> {dest}")
    urlretrieve(url, dest)
    print("Done.\n")

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data", help="Output directory (default: data)")
    ap.add_argument("--file", default=DEFAULT_SUPPL, help=f"Supplementary filename (default: {DEFAULT_SUPPL})")
    ap.add_argument("--series-matrix", action="store_true", help="Also download the GEO series matrix file")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    suppl_url = f"{BASE}/suppl/{args.file}"
    download(suppl_url, outdir / args.file)

    if args.series_matrix:
        sm_url = f"{BASE}/matrix/{SERIES_MATRIX}"
        download(sm_url, outdir / SERIES_MATRIX)

if __name__ == "__main__":
    main()

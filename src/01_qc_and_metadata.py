#!/usr/bin/env python3
"""
BISC 4171 - NeuroGene Explorer
01_qc_and_metadata.py

Loads a counts table, builds metadata from sample names (condition + genotype),
and outputs a basic QC plot (library size).

Usage:
  python src/01_qc_and_metadata.py --counts data/counts.txt.gz --outdir results

Outputs:
- results/metadata.csv
- results/qc_library_size.png

IMPORTANT:
You may need to adjust infer_condition() / infer_genotype() after inspecting sample names.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def infer_condition(sample_name: str) -> str:
    s = sample_name.lower()
    if "sd5" in s:
        return "SD5"
    if "hc5" in s:
        return "HC5"
    # fallback patterns
    if "sd" in s and "hc" not in s:
        return "SD5"
    if "hc" in s:
        return "HC5"
    return "UNKNOWN"

def infer_genotype(sample_name: str) -> str:
    s = sample_name.lower()
    # Flexible rules—adjust once you inspect actual sample names.
    if "shank3" in s or "mut" in s or "ko" in s or "delta" in s:
        return "Shank3"
    if "wt" in s:
        return "WT"
    return "UNKNOWN"

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="Path to counts table (tsv or tsv.gz)")
    ap.add_argument("--outdir", default="results", help="Output directory (default: results)")
    args = ap.parse_args()

    counts_path = Path(args.counts)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(counts_path, sep="\t", comment="#")
    gene_col = df.columns[0]
    counts = df.set_index(gene_col)

    samples = list(counts.columns)
    meta = pd.DataFrame({"sample": samples})
    meta["condition"] = meta["sample"].apply(infer_condition)
    meta["genotype"] = meta["sample"].apply(infer_genotype)

    meta.to_csv(outdir / "metadata.csv", index=False)

    libsize = counts.sum(axis=0).sort_values()

    plt.figure()
    plt.plot(range(len(libsize)), libsize.values)
    plt.xticks(range(len(libsize)), libsize.index, rotation=90, fontsize=7)
    plt.ylabel("Total counts per sample")
    plt.title("QC: library size per sample")
    plt.tight_layout()
    plt.savefig(outdir / "qc_library_size.png", dpi=200)

    print(f"Wrote: {outdir / 'metadata.csv'}")
    print(f"Wrote: {outdir / 'qc_library_size.png'}")

    if (meta["condition"] == "UNKNOWN").any():
        print("\nWARNING: condition=UNKNOWN found. Edit infer_condition() to match sample naming.")
        print(meta.loc[meta["condition"] == "UNKNOWN", "sample"].to_string(index=False))

    if (meta["genotype"] == "UNKNOWN").any():
        print("\nWARNING: genotype=UNKNOWN found. Edit infer_genotype() to match sample naming.")
        print(meta.loc[meta["genotype"] == "UNKNOWN", "sample"].to_string(index=False))

if __name__ == "__main__":
    main()

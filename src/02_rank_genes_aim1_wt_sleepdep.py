#!/usr/bin/env python3
"""
BISC 4171 - NeuroGene Explorer
02_rank_genes_aim1_wt_sleepdep.py

Aim 1 (required): WT SD5 vs WT HC5 gene ranking.

Beginner-friendly workflow:
- log2(counts + 1)
- Welch's t-test per gene
- Benjamini-Hochberg FDR correction
- Volcano plot + ranked table

Usage:
  python src/02_rank_genes_aim1_wt_sleepdep.py \
    --counts data/<COUNTS_FILE>.txt.gz \
    --metadata results/metadata.csv \
    --outdir results

Outputs:
- results/aim1_ranked_genes.csv
- results/aim1_volcano.png

Note:
This is a simplified learning workflow (not a full DESeq2/edgeR pipeline).
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

import pandas as pd

def clean_name(x):
    # strip whitespace, then strip surrounding quotes/apostrophes
    s = str(x).strip()
    return s.strip("'").strip('"').strip()

def load_counts_and_metadata(
    counts_path="data/GSE113754_RNASeq_S3_WT_SD5_HC5_Counts.txt.gz",
    metadata_path="results/metadata.csv",
    n_prefix_cols=2,   # counts files often start with GeneID, Length
    sample_col="sample"
):
    # --- read counts ---
    counts_df, meta_df = load_counts_and_metadata()
    counts.columns = [clean_name(c) for c in counts.columns]

    # --- read metadata ---
    counts_df, meta_df = load_counts_and_metadata()
    for col in [sample_col, "condition", "genotype"]:
        if col in meta.columns:
            meta[col] = meta[col].astype(str).map(clean_name)

    # --- sanity checks ---
    count_samples = set(counts.columns[n_prefix_cols:])
    meta_samples  = set(meta[sample_col])

    in_meta_not_counts = sorted(meta_samples - count_samples)
    in_counts_not_meta = sorted(count_samples - meta_samples)

    print("Example cleaned count columns:", counts.columns[:12].tolist())
    print("In metadata but not counts:", in_meta_not_counts[:20])
    print("In counts but not metadata:", in_counts_not_meta[:20])

    # --- align (subset + reorder) ---
    common = [s for s in meta[sample_col].tolist() if s in count_samples]
    if len(common) == 0:
        raise ValueError("No overlapping sample names between counts and metadata after cleaning.")

    # keep only samples present in metadata, ordered to match metadata
    counts_aligned = counts.iloc[:, :n_prefix_cols].copy()
    counts_aligned = pd.concat([counts_aligned, counts[common]], axis=1)

    # optionally filter metadata to only samples actually in counts
    meta_aligned = meta[meta[sample_col].isin(common)].copy()

    return counts_aligned, meta_aligned


# ---- in your script, right where you currently load files ----
counts_df, meta_df = load_counts_and_metadata()

# from here on, use counts_df + meta_df for Aim 1 filtering/ranking
# e.g., select WT samples, SD5 vs HC5, etc.


def load_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    gene_col = df.columns[0]
    return df.set_index(gene_col)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="Path to counts table (tsv or tsv.gz)")
    ap.add_argument("--metadata", required=True, help="Path to metadata.csv from 01_qc_and_metadata.py")
    ap.add_argument("--outdir", default="results", help="Output directory (default: results)")
    ap.add_argument("--min_total_counts", type=float, default=0.0,
                    help="Optional gene filter: keep genes with total counts >= this (default: 0)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts = load_counts(Path(args.counts))
    meta = pd.read_csv(args.metadata)

    # Ensure samples overlap
    meta = meta[meta["sample"].isin(counts.columns)].copy()

    # Aim 1: WT only (if genotype labels are present and WT exists)
    if "genotype" in meta.columns and (meta["genotype"].astype(str).str.upper() == "WT").any():
        meta = meta[meta["genotype"].astype(str).str.upper() == "WT"].copy()

    sd = meta.loc[meta["condition"] == "SD5", "sample"].tolist()
    hc = meta.loc[meta["condition"] == "HC5", "sample"].tolist()

    if len(sd) == 0 or len(hc) == 0:
        raise ValueError(
            "Could not find SD5 and HC5 samples for Aim 1.\n"
            "Check results/metadata.csv and adjust infer_condition() if needed."
        )

    if args.min_total_counts > 0:
        counts = counts.loc[counts.sum(axis=1) >= args.min_total_counts]

    logc = np.log2(counts + 1)
    log2fc = logc[sd].mean(axis=1) - logc[hc].mean(axis=1)

    pvals = []
    for gene_id in logc.index:
        p = ttest_ind(logc.loc[gene_id, sd], logc.loc[gene_id, hc], equal_var=False, nan_policy="omit").pvalue
        pvals.append(p)

    res = pd.DataFrame({
        "gene_id": logc.index,
        "log2FC_SD5_minus_HC5": log2fc.values,
        "pvalue": pvals
    }).set_index("gene_id")

    res["fdr"] = multipletests(res["pvalue"].fillna(1.0), method="fdr_bh")[1]
    res = res.sort_values(["fdr", "pvalue"])

    out_csv = outdir / "aim1_ranked_genes.csv"
    res.to_csv(out_csv)
    print("Wrote:", out_csv)

    plt.figure()
    x = res["log2FC_SD5_minus_HC5"].values
    y = -np.log10(np.clip(res["pvalue"].astype(float).values, 1e-300, 1.0))
    plt.scatter(x, y, s=4)
    plt.xlabel("log2FC (SD5 - HC5)")
    plt.ylabel("-log10(p-value)")
    plt.title("Aim 1: WT sleep deprivation (SD5 vs HC5)")
    plt.tight_layout()
    out_png = outdir / "aim1_volcano.png"
    plt.savefig(out_png, dpi=200)
    print("Wrote:", out_png)

    print("\nTop 10 genes by FDR:")
    print(res.head(10))

if __name__ == "__main__":
    main()

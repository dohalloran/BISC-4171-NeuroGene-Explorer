#!/usr/bin/env python3
"""
BISC 4171 - NeuroGene Explorer
03_rank_genes_aim2_genotype.py

Aim 2 (required): Genotype comparison (Shank3 vs WT) under a matched condition.

Default: compare within HC5. You can run within SD5 via --condition SD5.

Workflow:
- log2(counts + 1)
- Welch's t-test per gene
- BH-FDR correction
- Volcano plot + ranked table

Usage:
  python src/03_rank_genes_aim2_genotype.py \
    --counts data/<COUNTS_WITH_GENOTYPE>.txt.gz \
    --metadata results/metadata.csv \
    --condition HC5 \
    --outdir results

Outputs:
- results/aim2_ranked_genes_<CONDITION>.csv
- results/aim2_volcano_<CONDITION>.png
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def load_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    gene_col = df.columns[0]
    return df.set_index(gene_col)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="Path to counts table (tsv or tsv.gz)")
    ap.add_argument("--metadata", required=True, help="Path to metadata.csv from 01_qc_and_metadata.py")
    ap.add_argument("--condition", default="HC5", choices=["HC5", "SD5"], help="Condition to analyze (default: HC5)")
    ap.add_argument("--outdir", default="results", help="Output directory (default: results)")
    ap.add_argument("--min_total_counts", type=float, default=0.0,
                    help="Optional gene filter: keep genes with total counts >= this (default: 0)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts = load_counts(Path(args.counts))
    meta = pd.read_csv(args.metadata)

    if "genotype" not in meta.columns:
        raise ValueError("metadata.csv must include a 'genotype' column. Re-run 01_qc_and_metadata.py.")

    meta = meta[meta["sample"].isin(counts.columns)].copy()
    meta = meta[meta["condition"] == args.condition].copy()

    wt = meta.loc[meta["genotype"].astype(str).str.upper() == "WT", "sample"].tolist()
    sh = meta.loc[meta["genotype"].astype(str).str.lower().str.contains("shank3"), "sample"].tolist()

    if len(wt) == 0 or len(sh) == 0:
        raise ValueError(
            f"Could not find both WT and Shank3 samples in condition {args.condition}.\n"
            "Check results/metadata.csv and adjust infer_genotype() if needed."
        )

    if args.min_total_counts > 0:
        counts = counts.loc[counts.sum(axis=1) >= args.min_total_counts]

    logc = np.log2(counts + 1)
    log2fc = logc[sh].mean(axis=1) - logc[wt].mean(axis=1)

    pvals = []
    for gene_id in logc.index:
        p = ttest_ind(logc.loc[gene_id, sh], logc.loc[gene_id, wt], equal_var=False, nan_policy="omit").pvalue
        pvals.append(p)

    col_fc = f"log2FC_Shank3_minus_WT_{args.condition}"
    res = pd.DataFrame({
        "gene_id": logc.index,
        col_fc: log2fc.values,
        "pvalue": pvals
    }).set_index("gene_id")

    res["fdr"] = multipletests(res["pvalue"].fillna(1.0), method="fdr_bh")[1]
    res = res.sort_values(["fdr", "pvalue"])

    out_csv = outdir / f"aim2_ranked_genes_{args.condition}.csv"
    res.to_csv(out_csv)
    print("Wrote:", out_csv)

    plt.figure()
    x = res[col_fc].values
    y = -np.log10(np.clip(res["pvalue"].astype(float).values, 1e-300, 1.0))
    plt.scatter(x, y, s=4)
    plt.xlabel(f"log2FC (Shank3 - WT) [{args.condition}]")
    plt.ylabel("-log10(p-value)")
    plt.title(f"Aim 2: Genotype comparison within {args.condition}")
    plt.tight_layout()
    out_png = outdir / f"aim2_volcano_{args.condition}.png"
    plt.savefig(out_png, dpi=200)
    print("Wrote:", out_png)

    print("\nTop 10 genes by FDR:")
    print(res.head(10))

if __name__ == "__main__":
    main()

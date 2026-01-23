#!/usr/bin/env python3
"""
BISC 4171 - NeuroGene Explorer
03_rank_genes_aim2_genotype.py

Aim 2 (required): Genotype comparison (Shank3 vs WT) under a matched condition.

Default: compare Shank3 vs WT within HC5 (homecage).
You can also run within SD5 by setting --condition SD5.

Beginner-friendly workflow:
- log2(counts + 1)
- Welch's t-test per gene
- Benjamini-Hochberg FDR correction
- Volcano plot + ranked table

Why your volcano might look "flat":
Sometimes one gene gets an extremely tiny p-value (underflow to ~0), which makes
-log10(p) ~ 300 and squashes everything else. This script optionally clips the
y-axis so the plot is readable.

Usage:
  python src/03_rank_genes_aim2_genotype.py \
    --counts data/<COUNTS_FILE>.txt.gz \
    --metadata results/metadata.csv \
    --condition HC5 \
    --outdir results \
    --ymax 10 \
    --volcano_metric pvalue

Outputs (in --outdir):
- aim2_ranked_genes_<CONDITION>.csv
- aim2_volcano_<CONDITION>.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def clean_name(x) -> str:
    """Normalize sample/column names: trim whitespace and surrounding quotes/apostrophes."""
    s = str(x).strip()
    return s.strip("'").strip('"').strip()


def load_counts_and_metadata(
    counts_path: str | Path,
    metadata_path: str | Path,
    *,
    n_prefix_cols: int = 2,   # counts files often start with GeneID, Length
    sample_col: str = "sample",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """
    Load counts + metadata, clean sample names, and align counts columns to metadata order.

    Returns counts, meta, and optional prefix columns (e.g., Length).
    """
    counts_path = Path(counts_path)
    metadata_path = Path(metadata_path)

    raw = pd.read_csv(counts_path, sep="\t", compression="infer", comment="#")
    raw.columns = [clean_name(c) for c in raw.columns]

    if raw.shape[1] <= n_prefix_cols:
        raise ValueError(
            f"Counts table has {raw.shape[1]} columns but n_prefix_cols={n_prefix_cols}. "
            "If your counts file has fewer non-sample columns, lower --n_prefix_cols."
        )

    gene_col = raw.columns[0]
    raw[gene_col] = raw[gene_col].astype(str).map(clean_name)

    prefix_cols = list(raw.columns[:n_prefix_cols])

    raw = raw.set_index(gene_col)

    prefix = None
    prefix_cols_wo_gene = [c for c in prefix_cols if c != gene_col and c in raw.columns]
    if prefix_cols_wo_gene:
        prefix = raw[prefix_cols_wo_gene].copy()

    counts = raw.drop(columns=prefix_cols_wo_gene, errors="ignore").copy()
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    meta = pd.read_csv(metadata_path)
    if sample_col not in meta.columns:
        raise ValueError(f"Metadata file is missing required column: '{sample_col}'")

    for col in [sample_col, "condition", "genotype"]:
        if col in meta.columns:
            meta[col] = meta[col].astype(str).map(clean_name)

    # sanity checks
    count_samples = set(counts.columns)
    meta_samples = set(meta[sample_col])

    in_meta_not_counts = sorted(meta_samples - count_samples)
    in_counts_not_meta = sorted(count_samples - meta_samples)

    print("Example cleaned count columns:", list(counts.columns[:12]))
    print("In metadata but not counts:", in_meta_not_counts[:20])
    print("In counts but not metadata:", in_counts_not_meta[:20])

    # align
    common = [s for s in meta[sample_col].tolist() if s in count_samples]
    if len(common) == 0:
        raise ValueError("No overlapping sample names between counts and metadata after cleaning.")

    counts_aligned = counts[common].copy()
    meta_aligned = meta[meta[sample_col].isin(common)].copy()
    meta_aligned["_order"] = meta_aligned[sample_col].map({s: i for i, s in enumerate(common)})
    meta_aligned = meta_aligned.sort_values("_order").drop(columns=["_order"])

    return counts_aligned, meta_aligned, prefix


def volcano_plot(
    res: pd.DataFrame,
    x_col: str,
    metric: str,
    out_png: Path,
    title: str,
    xlabel: str,
    *,
    ymax: float = 10.0,
) -> None:
    x = res[x_col].astype(float).values
    p = res[metric].astype(float).clip(lower=1e-300, upper=1.0).values
    y_raw = -np.log10(p)
    y = np.minimum(y_raw, ymax)
    over = y_raw > ymax

    plt.figure()
    plt.scatter(x[~over], y[~over], s=4, alpha=0.6)
    if over.any():
        plt.scatter(x[over], np.full(over.sum(), ymax), s=10, marker="^", alpha=0.9)
        plt.text(
            0.99, 0.02,
            f"{over.sum()} genes clipped at y={ymax}",
            transform=plt.gca().transAxes,
            ha="right", va="bottom", fontsize=8,
        )

    plt.axvline(0, linewidth=1)
    plt.axhline(-np.log10(0.05), linewidth=1)

    plt.xlabel(xlabel)
    plt.ylabel(f"-log10({metric})")
    plt.title(title)
    plt.ylim(0, ymax * 1.05)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="Path to counts table (tsv or tsv.gz)")
    ap.add_argument("--metadata", required=True, help="Path to metadata.csv from 01_qc_and_metadata.py")
    ap.add_argument("--condition", default="HC5", choices=["HC5", "SD5"], help="Condition to analyze (default: HC5)")
    ap.add_argument("--outdir", default="results", help="Output directory (default: results)")
    ap.add_argument(
        "--min_total_counts",
        type=float,
        default=0.0,
        help="Optional gene filter: keep genes with total counts >= this (default: 0)",
    )
    ap.add_argument(
        "--n_prefix_cols",
        type=int,
        default=2,
        help="Number of non-sample leading columns in counts file (default: 2; e.g., GeneID, Length)",
    )
    ap.add_argument(
        "--volcano_metric",
        choices=["pvalue", "fdr"],
        default="pvalue",
        help="Volcano y-axis metric (default: pvalue). Try 'fdr' for adjusted significance.",
    )
    ap.add_argument(
        "--ymax",
        type=float,
        default=10.0,
        help="Max y-axis for volcano plot (clips extreme points). Default: 10",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts, meta, _prefix = load_counts_and_metadata(
        args.counts,
        args.metadata,
        n_prefix_cols=args.n_prefix_cols,
        sample_col="sample",
    )

    if "genotype" not in meta.columns:
        raise ValueError(
            "metadata.csv is missing a 'genotype' column.\n"
            "Re-run 01_qc_and_metadata.py and ensure infer_genotype() correctly labels WT vs Shank3."
        )
    if "condition" not in meta.columns:
        raise ValueError("Metadata is missing required column: 'condition'")

    # Filter to matched condition (HC5 default)
    meta = meta[meta["condition"] == args.condition].copy()

    wt_samples = meta.loc[meta["genotype"].astype(str).str.upper().eq("WT"), "sample"].tolist()
    sh_samples = meta.loc[meta["genotype"].astype(str).str.lower().str.contains("shank3"), "sample"].tolist()

    if len(wt_samples) == 0 or len(sh_samples) == 0:
        raise ValueError(
            f"Could not find both WT and Shank3 samples for condition {args.condition}.\n"
            "Open results/metadata.csv and inspect sample names; adjust infer_genotype() if needed."
        )

    # Optional gene filter
    if args.min_total_counts > 0:
        counts = counts.loc[counts.sum(axis=1) >= args.min_total_counts]

    logc = np.log2(counts + 1.0)

    log2fc = logc[sh_samples].mean(axis=1) - logc[wt_samples].mean(axis=1)

    pvals = []
    for gene_id in logc.index:
        p = ttest_ind(
            logc.loc[gene_id, sh_samples],
            logc.loc[gene_id, wt_samples],
            equal_var=False,
            nan_policy="omit",
        ).pvalue
        pvals.append(p)

    x_col = f"log2FC_Shank3_minus_WT_{args.condition}"

    res = pd.DataFrame(
        {
            x_col: log2fc.values,
            "pvalue": pvals,
        },
        index=logc.index,
    )

    res["pvalue"] = res["pvalue"].astype(float).clip(lower=1e-300, upper=1.0)
    res["fdr"] = multipletests(res["pvalue"].fillna(1.0), method="fdr_bh")[1]
    res = res.sort_values(["fdr", "pvalue"])

    out_csv = outdir / f"aim2_ranked_genes_{args.condition}.csv"
    res.to_csv(out_csv)
    print("Wrote:", out_csv)

    out_png = outdir / f"aim2_volcano_{args.condition}.png"
    volcano_plot(
        res,
        x_col=x_col,
        metric=args.volcano_metric,
        out_png=out_png,
        title=f"Aim 2: Genotype comparison within {args.condition}",
        xlabel=f"log2FC (Shank3 - WT) [{args.condition}]",
        ymax=args.ymax,
    )
    print("Wrote:", out_png)

    print("\nTop 10 genes by FDR:")
    print(res.head(10))


if __name__ == "__main__":
    main()

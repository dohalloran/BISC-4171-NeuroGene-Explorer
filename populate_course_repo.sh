#!/usr/bin/env bash
set -euo pipefail

# BISC 4171 - NeuroGene Explorer
# Populate the *current* repository directory with scaffold + Aim 1/Aim 2 scripts.
#
# Run this from the repo root (the folder that contains the .git directory):
#   bash populate_course_repo.sh
#
# Then commit + push:
#   git add .
#   git commit -m "Add course scaffold + Aim 1/Aim 2 scripts"
#   git push -u origin main

if [[ ! -d ".git" ]]; then
  echo "ERROR: I don't see a .git folder. Please cd into your cloned repo root and run again."
  exit 1
fi

echo "Populating repo in: $(pwd)"

# Folders
mkdir -p src docs assignments data results

# .gitignore
cat > .gitignore <<'EOF'
# Do not commit data or derived results
data/
results/

# Python
__pycache__/
*.py[cod]
*.egg-info/
dist/
build/

# Virtual environments
.venv/
venv/
.env

# OS/editor
.DS_Store
Thumbs.db
.vscode/
.idea/

# Jupyter
.ipynb_checkpoints/
EOF

# environment.yml
cat > environment.yml <<'EOF'
name: bisc4171-neurogene
channels:
  - conda-forge
dependencies:
  - python=3.11
  - pip
  - numpy
  - pandas
  - matplotlib
  - scipy
  - statsmodels
EOF

# LICENSE
cat > LICENSE <<'EOF'
MIT License

Copyright (c) 2026

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF

# README.md (Aim 1 + Aim 2 only; no AI log requirement)
cat > README.md <<'EOF'
# BISC 4171 — NeuroGene Explorer (Asynchronous)

A lightweight, semester-long undergraduate research project introducing **Python + bioinformatics** through a **neural gene-expression** dataset.

## Course structure (1 credit)
- **No scheduled meeting time** (asynchronous / independent research)
- Major deliverables:
  - **Midterm Progress Report** (Aim 1)
  - **Final Report** (Aim 1 + Aim 2)

## Dataset
Primary dataset: **NCBI GEO GSE113754** (mouse prefrontal cortex RNA-seq; sleep deprivation vs homecage; WT and Shank3 mutant).

## Aims (required)
1. **Aim 1 (WT sleep deprivation):** WT SD5 vs WT HC5  
2. **Aim 2 (genotype comparison):** Shank3 vs WT within a matched condition (minimum: HC5)

## Repository layout
- `src/` — analysis scripts (Python)
- `docs/` — PDFs for students (syllabus, project packet, rubric, etc.)
- `assignments/` — prompts for midterm/final reports (optional)
- `data/` — **not tracked** (raw downloads live here)
- `results/` — **not tracked** (generated outputs live here)

## Setup (students)
Recommended: Miniconda / Mambaforge

```bash
conda env create -f environment.yml
conda activate bisc4171-neurogene
```

## Running the pipeline (example)
```bash
# 1) Download processed counts (no FASTQ processing)
python src/00_download_gse113754.py --outdir data --series-matrix

# 2) Build metadata + QC figure
python src/01_qc_and_metadata.py --counts data/<COUNTS_FILE>.txt.gz --outdir results

# 3) Aim 1: WT SD5 vs WT HC5
python src/02_rank_genes_aim1_wt_sleepdep.py --counts data/<COUNTS_FILE>.txt.gz --metadata results/metadata.csv --outdir results

# 4) Aim 2: Shank3 vs WT within HC5 (or SD5 if desired)
python src/03_rank_genes_aim2_genotype.py --counts data/<COUNTS_WITH_GENOTYPE>.txt.gz --metadata results/metadata.csv --condition HC5 --outdir results
```

## Notes on AI tools
Students may use AI tools to help **debug code**, but there is **no AI log requirement**. Students remain responsible for understanding and explaining their final code and results.

## License
Code: MIT (see `LICENSE`).
EOF

# docs placeholder
cat > docs/README.md <<'EOF'
# Course documents

Put course PDFs here (syllabus, project packet, genotype module, rubric, etc.).
EOF

# assignments placeholder
cat > assignments/README.md <<'EOF'
# Assignments (optional)

Add prompts here for:
- Midterm Progress Report (Aim 1)
- Final Report (Aim 1 + Aim 2)
EOF

# ---------------------------
# Scripts (Aim 1 + Aim 2 only)
# ---------------------------

cat > src/00_download_gse113754.py <<'EOF'
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
EOF
chmod +x src/00_download_gse113754.py

cat > src/01_qc_and_metadata.py <<'EOF'
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
EOF
chmod +x src/01_qc_and_metadata.py

cat > src/02_rank_genes_aim1_wt_sleepdep.py <<'EOF'
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
EOF
chmod +x src/02_rank_genes_aim1_wt_sleepdep.py

cat > src/03_rank_genes_aim2_genotype.py <<'EOF'
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
EOF
chmod +x src/03_rank_genes_aim2_genotype.py

echo
echo "✅ Repo files written."
echo "Next:"
echo "  git add ."
echo "  git commit -m "Add course scaffold + Aim 1/Aim 2 scripts""
echo "  git push -u origin main"

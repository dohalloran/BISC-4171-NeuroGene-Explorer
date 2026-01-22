# BISC 4171 — NeuroGene Explorer 

An undergraduate research project introducing **Python + bioinformatics** through a **neural gene-expression** dataset.

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

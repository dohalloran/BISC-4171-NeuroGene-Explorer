# BISC 4171 — NeuroGene Explorer

An undergraduate research project introducing **Python + bioinformatics** through a **neural gene-expression** dataset.

> Course repo: https://github.com/dohalloran/BISC-4171-NeuroGene-Explorer  
> (Students can **download the ZIP** rather than using Git.)

---

## Course structure (1 credit)
- **No scheduled meeting time** (asynchronous / independent research)
- Major deliverables:
  - **Midterm Progress Report** (Aim 1)
  - **Final Report** (Aim 1 + Aim 2)
  - **Lightning talk** (one-slide PDF + 2-3minute recording)

---

## Dataset
Primary dataset: **NCBI GEO GSE113754**  
Mouse prefrontal cortex RNA-seq; sleep deprivation vs homecage; WT and Shank3 mutant.

---

## Aims (required)
1. **Aim 1 (WT sleep deprivation):** WT SD5 vs WT HC5  
2. **Aim 2 (genotype comparison):** Shank3 vs WT within a matched condition (**minimum: HC5**)

---

## Repository layout
- `src/` — analysis scripts (Python)
- `docs/` — course PDFs (syllabus, project packet, rubric, etc.)
- `assignments/` — prompts for midterm/final reports
- `data/` — **not tracked** (raw downloads live here)
- `results/` — **not tracked** (generated outputs live here)

---

## Quick start (students) — ZIP download (no Git)

### Step 1 — Download the repo ZIP
1. Go to the repo page: https://github.com/dohalloran/BISC-4171-NeuroGene-Explorer  
2. Click **Code** → **Download ZIP**
3. Unzip the folder somewhere you can find easily (Desktop is fine).

### Step 2 — Install Miniconda or Mambaforge
Recommended: Miniconda (or Mambaforge). Either works.

### Step 3 — Open a terminal in the unzipped folder

#### Mac (Terminal)
- Finder → open the unzipped folder
- Right-click inside the folder → **New Terminal at Folder** (or open Terminal and `cd` into the folder)

#### Windows (Anaconda Prompt)
- Start menu → open **Anaconda Prompt**
- `cd` into the unzipped folder (example):
  ```bat
  cd %USERPROFILE%\Desktop\BISC-4171-NeuroGene-Explorer
  ```

### Step 4 — Create and activate the environment
```bash
conda env create -f environment.yml
conda activate bisc4171-neurogene
```

### Step 5 — Create output folders
```bash
mkdir -p data results
```

Windows PowerShell (if `mkdir -p` fails):
```powershell
mkdir data
mkdir results
```

---

## Running the pipeline

### 1) Download processed counts (no FASTQ processing)
```bash
python src/00_download_gse113754.py --outdir data --series-matrix
```

### 2) Build metadata + QC figure
Replace `<COUNTS_FILE>.txt.gz` with the file created in `data/`.
```bash
python src/01_qc_and_metadata.py --counts data/<COUNTS_FILE>.txt.gz --outdir results
```

### 3) Aim 1: WT SD5 vs WT HC5
```bash
python src/02_rank_genes_aim1_wt_sleepdep.py \
  --counts data/GSE113754_RNASeq_S3_WT_SD5_HC5_Counts.txt.gz \
  --metadata results/metadata.csv \
  --outdir results \
  --volcano_metric pvalue \
  --ymax 8
```

### 4) Aim 2: Shank3 vs WT within HC5
If your downloaded counts include genotype labels, run:
```bash
python src/03_rank_genes_aim2_genotype.py \
  --counts data/<COUNTS_WITH_GENOTYPE>.txt.gz \
  --metadata results/metadata.csv \
  --condition HC5 \
  --outdir results
```

---

## What to submit (students)
- Midterm report (PDF)
- Final report (PDF)
- Lightning talk (one-slide PDF)
- Code + reproducibility package (ZIP) including:
  - `src/` scripts
  - `results/` outputs
  - a short README or text file with the exact commands you ran  
  (Do **not** submit large raw downloads unless asked.)

---

## Notes on AI tools
Students may use AI tools to help **debug code**.  
Students remain responsible for understanding and explaining their final code and results.

---

## License
Code: MIT (see `LICENSE`).

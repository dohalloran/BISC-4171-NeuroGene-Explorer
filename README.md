# BISC 4171 — NeuroGene Explorer

An undergraduate research project introducing **Python + bioinformatics** through a **neural gene-expression** dataset.

---

## Course structure (1 credit)
- **No scheduled meeting time** (asynchronous / independent research)
- Major deliverables:
  - **Midterm Progress Report** (Aim 1)
  - **Final Report** (Aim 1 + Aim 2)

---

## Dataset
Primary dataset: **NCBI GEO GSE113754** (mouse prefrontal cortex RNA-seq; sleep deprivation vs homecage; WT and Shank3 mutant).

---

## Aims (required)
1. **Aim 1 (WT sleep deprivation):** WT SD5 vs WT HC5  
2. **Aim 2 (genotype comparison):** Shank3 vs WT within a matched condition (minimum: HC5)

---

## Repository layout
- `src/` — analysis scripts (Python)
- `docs/` — PDFs for students (syllabus, project packet, rubric, etc.)
- `assignments/` — prompts for midterm/final reports 
- `data/` — **not tracked** (raw downloads live here)
- `results/` — **not tracked** (generated outputs live here)

---

# Quick Commands (copy/paste)

## Step 1 — Download + unzip the repo (no Git required)
From GitHub: **Code → Download ZIP**, unzip it.

You will get a folder named something like:
- `BISC-4171-NeuroGene-Explorer-main`

---

## Step 2 — Open a terminal and go into the project folder

### Mac (Terminal)
Example if you unzipped to Desktop:
```bash
cd ~/Desktop/BISC-4171-NeuroGene-Explorer-main
````

### Windows (Anaconda Prompt)

Example if you unzipped to Desktop:

```bash
cd Desktop\BISC-4171-NeuroGene-Explorer-main
```

**Confirm you’re in the right place:**

### Mac

```bash
ls
```

### Windows (Anaconda Prompt)

```bash
dir
```

You should see folders like `src`, `docs`, and a file called `environment.yml`.

---

## Step 3 — Create the environment (one time)

```bash
conda env create -f environment.yml
```

## Step 4 — Activate the environment (every time you work)

```bash
conda activate bisc4171-neurogene
```

Optional check:

```bash
python --version
```

---

## Step 5 — Download the processed counts (fast)

```bash
python src/00_download_gse113754.py --outdir data --series-matrix
```

## Step 6 — List files in `data/` and copy the counts filename

### Mac

```bash
ls data
```

### Windows

```bash
dir data
```

---

## Step 7 — QC + metadata (replace the filename!)

Replace `PASTE_COUNTS_FILENAME_HERE.txt.gz` with the exact filename you saw in Step 6.

```bash
python src/01_qc_and_metadata.py --counts data/PASTE_COUNTS_FILENAME_HERE.txt.gz --outdir results
```

---

## Step 8 — Aim 1 (required): WT sleep deprivation

```bash
python src/02_rank_genes_aim1_wt_sleepdep.py --counts data/PASTE_COUNTS_FILENAME_HERE.txt.gz --metadata results/metadata.csv --outdir results
```

---

## Step 9 — Aim 2 (required): Genotype comparison (Shank3 vs WT within HC5)

```bash
python src/03_rank_genes_aim2_genotype.py --counts data/PASTE_COUNTS_FILENAME_HERE.txt.gz --metadata results/metadata.csv --condition HC5 --outdir results
```

---

## Step 10 — Check your outputs

### Mac

```bash
ls results
```

### Windows

```bash
dir results
```

You should see:

* `metadata.csv`, `qc_library_size.png`
* `aim1_ranked_genes.csv`, `aim1_volcano.png`
* `aim2_ranked_genes_HC5.csv`, `aim2_volcano_HC5.png`

---

# Step-by-step: Student Quickstart (more detail)

## Step 0 — Install prerequisites (one time)

You need:

* **Miniconda** or **Mambaforge** (to manage Python packages)
* A terminal:

  * **Mac:** Terminal
  * **Windows:** Anaconda Prompt (recommended)

If you already have conda installed, move on to Step 1.

---

## Step 1 — Download the project as a ZIP

1. Go to the GitHub repository page (your instructor will provide the link).
2. Click the green **Code** button.
3. Click **Download ZIP**.
4. Unzip the folder somewhere you can find it (Desktop is fine).

---

## Step 2 — Navigate into the project folder

### Mac (Terminal)

Tip: type `cd ` (with a trailing space) and drag the folder into the Terminal window.

Example:

```bash
cd /Users/yourname/Desktop/BISC-4171-NeuroGene-Explorer-main
```

### Windows (Anaconda Prompt)

Example:

```bash
cd Desktop\BISC-4171-NeuroGene-Explorer-main
```

Confirm:

* Mac: `ls`
* Windows: `dir`

---

## Step 3 — Create and activate the course Python environment

Create (one time):

```bash
conda env create -f environment.yml
```

Activate (every session):

```bash
conda activate bisc4171-neurogene
```

---

## Step 4 — Download the dataset (processed counts)

We use **processed counts** (not raw FASTQ), so downloads are quick and beginner-friendly:

```bash
python src/00_download_gse113754.py --outdir data --series-matrix
```

This creates files inside `data/`.

---

## Step 5 — Build metadata + run a basic QC check

This step creates a `metadata.csv` file that the later steps use.

```bash
python src/01_qc_and_metadata.py --counts data/<COUNTS_FILE>.txt.gz --outdir results
```

### Important: Replace `<COUNTS_FILE>.txt.gz`

List the downloaded files:

* Mac: `ls data`
* Windows: `dir data`

Then copy/paste the exact filename into the command.

Outputs in `results/`:

* `metadata.csv`
* `qc_library_size.png`

If `metadata.csv` contains `UNKNOWN` values for genotype or condition, you may need to adjust the sample-name parsing rules inside:

* `src/01_qc_and_metadata.py`

(Instructions are commented inside that file.)

---

## Step 6 — Aim 1 (required): WT sleep deprivation (SD5 vs HC5)

Run:

```bash
python src/02_rank_genes_aim1_wt_sleepdep.py --counts data/<COUNTS_FILE>.txt.gz --metadata results/metadata.csv --outdir results
```

Outputs:

* `results/aim1_ranked_genes.csv`
* `results/aim1_volcano.png`

---

## Step 7 — Aim 2 (required): Genotype comparison (Shank3 vs WT)

Required minimum: run within **HC5**.

```bash
python src/03_rank_genes_aim2_genotype.py --counts data/<COUNTS_FILE>.txt.gz --metadata results/metadata.csv --condition HC5 --outdir results
```

Outputs:

* `results/aim2_ranked_genes_HC5.csv`
* `results/aim2_volcano_HC5.png`

---

## Step 8 — What to include in your reports

### Midterm Progress Report (Aim 1)

Include:

* `aim1_ranked_genes.csv` (or top 20–50 genes summarized)
* `aim1_volcano.png`
* A short interpretation: What changed with sleep deprivation? Which genes stand out and why?

### Final Report (Aim 1 + Aim 2)

Include:

* Aim 1 outputs (above)
* Aim 2 outputs:

  * `aim2_ranked_genes_HC5.csv`
  * `aim2_volcano_HC5.png`
* A comparison paragraph:

  * How do the “sleep deprivation” results (Aim 1) relate to the “genotype” results (Aim 2)?
  * Do any genes/pathways overlap?

---

## Notes on AI tools

Students may use AI tools to help **debug code**. Students remain responsible for understanding and explaining their final code and results.

---

## License

Code: MIT (see `LICENSE`).


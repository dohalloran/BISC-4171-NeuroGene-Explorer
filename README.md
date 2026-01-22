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
- `docs/` — PDFs for students
- `assignments/` — prompts for midterm/final reports 
- `data/` — raw downloads live here
- `results/` — generated outputs live here

---

# Step-by-step: Student Quickstart (ZIP download)

## Step 0 — Install prerequisites (one time)
You need:
- **Miniconda** or **Mambaforge** (to manage Python packages)
- A terminal:
  - **Mac:** Terminal
  - **Windows:** Anaconda Prompt (recommended) or PowerShell

If you already have conda installed, move on to Step 1.

---

## Step 1 — Download the project as a ZIP (no Git required)
1. Go to the GitHub repository page (your instructor will provide the link).
2. Click the green **Code** button.
3. Click **Download ZIP**.
4. Unzip the folder somewhere you can find it (Desktop is fine).

You should now have a folder named something like:
- `BISC-4171-NeuroGene-Explorer-main`

---

## Step 2 — Open a terminal in the project folder

### Mac
1. Open **Terminal**
2. Type `cd ` (with a trailing space), then drag the unzipped folder into the Terminal window and press Enter.

Example:
```bash
cd /Users/yourname/Desktop/BISC-4171-NeuroGene-Explorer-main
````

### Windows

1. Open **Anaconda Prompt**
2. Navigate to the folder, for example:

```bash
cd Desktop\BISC-4171-NeuroGene-Explorer-main
```

To confirm you’re in the right place, run:

```bash
ls
```

(or on Windows PowerShell you can use `dir`)

You should see folders like `src`, `docs`, and the file `environment.yml`.

---

## Step 3 — Create and activate the course Python environment

From the project folder:

```bash
conda env create -f environment.yml
conda activate bisc4171-neurogene
```

Verify Python is working:

```bash
python --version
```

---

## Step 4 — Create the folders for outputs (if needed)

These usually already exist, but it’s safe to run:

```bash
mkdir -p data results
```

---

## Step 5 — Download the dataset (processed counts)

We use **processed counts** (not raw FASTQ), so downloads are quick and beginner-friendly:

```bash
python src/00_download_gse113754.py --outdir data --series-matrix
```

This should create files inside `data/`.

---

## Step 6 — Build metadata + run a basic QC check

This step creates a `metadata.csv` file that the later steps use.

```bash
python src/01_qc_and_metadata.py --counts data/<COUNTS_FILE>.txt.gz --outdir results
```

### Important: Replace `<COUNTS_FILE>.txt.gz`

After Step 5, list the downloaded files:

```bash
ls data
```

Then copy/paste the exact filename into Step 6.

Outputs you should get in `results/`:

* `metadata.csv`
* `qc_library_size.png`

If `metadata.csv` contains `UNKNOWN` values for genotype or condition, you may need to adjust the sample-name parsing rules inside:

* `src/01_qc_and_metadata.py`

(Instructions are commented inside that file.)

---

## Step 7 — Aim 1 (required): WT sleep deprivation (SD5 vs HC5)

Run:

```bash
python src/02_rank_genes_aim1_wt_sleepdep.py \
  --counts data/<COUNTS_FILE>.txt.gz \
  --metadata results/metadata.csv \
  --outdir results
```

Outputs you should get:

* `results/aim1_ranked_genes.csv`
* `results/aim1_volcano.png`

---

## Step 8 — Aim 2 (required): Genotype comparison (Shank3 vs WT)

By default, this runs within **HC5** (required minimum). You can also run SD5 if you want.

```bash
python src/03_rank_genes_aim2_genotype.py \
  --counts data/<COUNTS_WITH_GENOTYPE>.txt.gz \
  --metadata results/metadata.csv \
  --condition HC5 \
  --outdir results
```

Outputs you should get:

* `results/aim2_ranked_genes_HC5.csv`
* `results/aim2_volcano_HC5.png`

---

## Step 9 — What to include in your reports

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

Students may use AI tools to help **debug code**, but there is **no AI log requirement**. Students remain responsible for understanding and explaining their final code and results.

---

## Optional (advanced): using Git instead of ZIP

If you want to learn Git later, you can clone the project instead of downloading a ZIP:

```bash
git clone https://github.com/dohalloran/BISC-4171-NeuroGene-Explorer.git
cd BISC-4171-NeuroGene-Explorer
```

---

## License

Code: MIT (see `LICENSE`).


#!/usr/bin/env bash
set -euo pipefail

# BISC 4171 - NeuroGene Explorer
# Repo scaffolding + optional GitHub creation/push (requires GitHub CLI "gh")
#
# Usage:
#   bash init_course_repo.sh
#
# Notes:
# - If you have the GitHub CLI installed and authenticated ("gh auth login"),
#   this script can create the repo and push automatically.
# - If not, it will still scaffold the repo locally and tell you what to do next.

REPO_DEFAULT="bisc4171-neurogene-explorer"

echo "=== BISC 4171 Course Repo Scaffolder ==="
read -rp "Repo name [$REPO_DEFAULT]: " REPO_NAME
REPO_NAME="${REPO_NAME:-$REPO_DEFAULT}"

read -rp "Visibility (public/private) [public]: " VIS
VIS="${VIS:-public}"
if [[ "$VIS" != "public" && "$VIS" != "private" ]]; then
  echo "Visibility must be 'public' or 'private'."
  exit 1
fi

if [[ -e "$REPO_NAME" ]]; then
  echo "Folder '$REPO_NAME' already exists. Refusing to overwrite."
  exit 1
fi

mkdir -p "$REPO_NAME"
cd "$REPO_NAME"

# --- folders ---
mkdir -p src docs templates assignments examples data results .github/ISSUE_TEMPLATE

# --- basic files ---
cat > .gitignore <<'EOF'
# Data + outputs (do not commit)
data/
results/

# Python
__pycache__/
*.pyc
*.pyo
*.pyd
*.egg-info/
dist/
build/

# Environments
.venv/
venv/
.env
.envrc

# OS/editor
.DS_Store
Thumbs.db
.vscode/
.idea/

# Jupyter
.ipynb_checkpoints/
EOF

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
  - pip:
      - requests
EOF

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

cat > docs/README.md <<'EOF'
# Course documents

Put your PDFs here (syllabus, project packet, genotype module, rubric, etc.).

Suggested naming:
- `BISC4171_Syllabus_Sp2026.pdf`
- `BISC4171_NeuroGene_Explorer_Student_Project_Packet.pdf`
- `BISC4171_NeuroGene_Explorer_Genotype_Comparison_Module.pdf`
- `BISC4171_NeuroGene_Explorer_Instructor_Rubric.pdf`
EOF

cat > templates/AI_Debugging_Log.md <<'EOF'
# AI Debugging Log (Template)

If you use AI tools to debug code, include a short entry per interaction.

| Date | Script / location | Error or goal | AI tool | Prompt (or summary) | AI response (or summary) | What you changed | What you learned |
|---|---|---|---|---|---|---|---|
| YYYY-MM-DD | 01_qc_and_metadata.py: infer_genotype() | Example: genotype labels not detected | ChatGPT | ... | ... | ... | ... |
EOF

cat > assignments/README.md <<'EOF'
# Assignments

This folder can hold:
- Midterm Progress Report prompt + rubric
- Final Report prompt + rubric
- Lightning talk instructions
EOF

cat > .github/ISSUE_TEMPLATE/bug_report.md <<'EOF'
---
name: Bug report
about: Report a problem running the pipeline
title: "[BUG] "
labels: bug
assignees: ""
---

## What happened?
Describe the error and what you expected.

## Command you ran
Paste the exact command.

## Error output
Paste the full error message/traceback.

## Your environment
- OS:
- Python version:
- Conda env name:
EOF

# --- README.md (if not provided separately) ---
if [[ ! -f README.md ]]; then
  cat > README.md <<'EOF'
# BISC 4171 — NeuroGene Explorer (Asynchronous)

A lightweight, semester-long undergraduate research project introducing **Python + bioinformatics** through a **neural gene-expression** dataset.

## Course structure (1 credit)
- **No scheduled meeting time** (asynchronous / independent research)
- A few high-impact deliverables:
  - **Midterm Progress Report** (Aim 1)
  - **Final Report** (Aim 1 + Aim 2 genotype required)
  - **Lightning talk** (one slide PDF; optional short recording)
  - **Code + reproducibility package**

## Dataset
Primary dataset: **NCBI GEO GSE113754** (mouse prefrontal cortex RNA-seq; sleep deprivation vs homecage; WT and Shank3 mutant).

## Aims (required)
1. **Aim 1 (WT sleep deprivation):** WT SD5 vs WT HC5
2. **Aim 2 (genotype comparison required):** Shank3 vs WT within a matched condition (minimum: HC5)

## Repository layout
- `src/` — analysis scripts (Python)
- `docs/` — PDFs for students (syllabus, packets, rubric)
- `templates/` — report/log templates (e.g., AI Debugging Log)
- `data/` — **not tracked** (raw downloads live here)
- `results/` — **not tracked** (generated outputs live here)

## Setup (students)
Recommended: Miniconda / Mambaforge
```bash
conda env create -f environment.yml
conda activate bisc4171-neurogene


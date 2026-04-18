# Step 2 — Data Preparation

**LB2 Project · Group 7 · Signal Peptide Prediction**

Raw sequence datasets from public databases almost always contain redundancy — many proteins share high sequence identity simply because they were independently deposited from related organisms or splice variants. Training on a redundant dataset inflates performance metrics and risks data leakage between train and test. This step removes that redundancy and produces the clean, non-overlapping splits that all downstream models are trained and evaluated on.

---

## Contents

| File                        | Description |
|-----------------------------|-------------|
| `step2_data_preparation.ipynb` | Full pipeline: MMseqs2 clustering, metadata filtering, labeling, 80/20 split, and 5-fold CV assignment |
| `filtered_positive.tsv`     | Positive cluster representatives with metadata (1,093 entries) |
| `filtered_negative.tsv`     | Negative cluster representatives with metadata (8,934 entries) |
| `training_with_folds.tsv`   | Training set (80%) with `label` and `fold` columns (8,021 entries) |
| `benchmarking_set.tsv`      | Held-out benchmark set (20%) with `label` column — **never used during training** (2,006 entries) |

---

## Pipeline Overview

1. **Redundancy removal (MMseqs2)** — Positive and negative sets are clustered independently using `easy-cluster`. One representative per cluster is retained. We use MMseqs2 rather than CD-HIT because it scales better to large datasets and allows more flexible coverage-mode settings.
2. **TSV filtering** — The metadata files from Step 1 are filtered to keep only the cluster representatives, so sequence and metadata files stay in sync.
3. **Labeling** — A binary `label` column is added (`1` = SP-positive, `0` = SP-negative).
4. **Stratified 80/20 split** — The split is performed separately within each class, preserving the natural class ratio in both train and benchmark sets. The benchmark set is sealed at this point and never touched again until final evaluation.
5. **5-fold cross-validation** — Fold assignments are created only for the training set, again stratified by label. These folds are reused across Steps 4, 5, and 6 to ensure all models are compared under identical data conditions.

**Reproducibility:** `RANDOM_SEED = 42` is used throughout. Do not change this value when re-running.

---

## MMseqs2 Clustering Parameters

| Parameter      | Value | Meaning |
|----------------|-------|---------|
| `--min-seq-id` | 0.3   | Minimum 30% sequence identity |
| `-c`           | 0.4   | Minimum 40% coverage |
| `--cov-mode`   | 0     | Coverage computed over both query and target |
| `--cluster-mode` | 1   | Connected-component clustering |

A 30% identity threshold is a standard choice for removing redundancy in protein classification tasks — stringent enough to prevent similar sequences from appearing in both train and test, while still leaving enough diversity for the model to generalise.

---

## Output Format

All TSV files share the same base columns, with extra fields depending on class:

| Column                | Description |
|-----------------------|-------------|
| `Accession`           | UniProt accession |
| `Organism`            | Source organism |
| `Kingdom`             | Metazoa / Viridiplantae / Fungi / Other |
| `Sequence length`     | Full protein length |
| `SP cleavage`         | Cleavage site (positive only) |
| `N-term transmembrane`| N-terminal TM helix flag (negative only) |
| `label`               | `1` = SP-positive, `0` = SP-negative |
| `fold`                | CV fold (0–4) — only in `training_with_folds.tsv` |

---

## Dataset Statistics

| Step                          | Positive | Negative | Total  |
|-------------------------------|----------|----------|--------|
| Raw (from Step 1)             | 2,932    | 20,615   | 23,547 |
| After MMseqs2 clustering      | **1,093**| **8,934**| 10,027 |
| Training set (80%)            | 874      | 7,147    | **8,021** |
| Benchmarking set (20%)        | 219      | 1,787    | **2,006** |

Class imbalance ≈ 8.2 negatives per positive. This imbalance is real and biologically expected — signal peptides are present in a minority of the proteome. It is handled explicitly in Steps 5 and 6 rather than resampled away here, to avoid distorting the evaluation conditions.

---

## How to Run

1. Upload the four files from the **Data Collection** folder:
   - `positive.fasta`, `negative.fasta`
   - `positive.tsv`, `negative.tsv`
2. Run `step2_data_preparation.ipynb` in Google Colab.
3. All output files will be generated automatically.

**Next step:** Proceed to the exploratory analysis notebook in `3- Data_analysis/`.

---

**Dependencies**
- Python 3.x
- `mmseqs2` (installed via `apt-get`)
- `biopython`
- `pandas`, `numpy`

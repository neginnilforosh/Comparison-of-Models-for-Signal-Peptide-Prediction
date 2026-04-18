# Step 5 — SVM Signal Peptide Classifier

**LB2 Project · Group 7 · Signal Peptide Prediction**

## Overview

The Von Heijne method in Step 4 works purely from the cleavage-site window — it has no knowledge of the broader N-terminal hydrophobic context that defines a signal peptide. This step asks whether hand-crafted biochemical features extracted from the N-terminal region can close that gap. An SVM with an RBF kernel is a natural fit here: it can capture non-linear relationships between features, handles moderate class imbalance reasonably well, and is well understood in terms of what it can and cannot learn. We also compare a full feature set (28 features) against a reduced set selected by Random Forest importance, to check whether simpler representations are competitive.

## Contents

| File | Description |
|------|-------------|
| `step5_svm.ipynb` | Main notebook — full pipeline from feature extraction to evaluation |
| `figures/svm_feature_importance.pdf/.png` | Top 20 features by Random Forest importance (training set) |
| `figures/svm_cv_pr_roc_all.pdf/.png` | 5-fold OOF PR and ROC curves (all features) |
| `figures/svm_benchmark_confusion.pdf/.png` | Benchmark confusion matrices (all vs selected features) |

## Method

1. **Data loading** — Training and benchmark splits are loaded and restricted to the filtered non-redundant accession set from Step 2.
2. **Feature extraction** — 28 biochemical features are computed from the N-terminal region of each sequence (see feature table below). All features were chosen to reflect known signal peptide biology: hydrophobicity, charge distribution, and secondary structure propensity of the N-terminus.
3. **Outer cross-validation loop** — For each fold:
   - Optionally select top features using **Random Forest importance** fit on the outer-training split only (no leakage from the validation fold).
   - **Z-score standardise** using a scaler fit on the outer-training split only.
   - Tune RBF kernel hyperparameters (**C**, **γ**) via **inner cross-validation**, optimising MCC.
   - Evaluate on the held-out fold to collect OOF predictions.
4. **Two configurations compared:** all 28 features vs. top 20 selected features. This tests whether the feature selection step is adding value or just adding complexity.
5. **Final model** — Hyperparameters re-selected on the full training set via 5-fold CV; model retrained and evaluated once on the blind benchmark.

## Features Extracted

| Feature group | Features | Region |
|---|---|---|
| AA composition | Frequency of each of 20 amino acids (`comp_*`) | First 20 aa |
| Hydrophobicity | Max and avg (Kyte-Doolittle, window=5) | First 40 aa |
| Charge | Max K/R abundance and its normalised position (window=3) | First 40 aa |
| α-helix propensity | Max and avg (Chou-Fasman, window=7) | First 40 aa |
| TM propensity | Max and avg (Kyte scale, window=7) | First 40 aa |

**Total: 28 features** (20 composition + 2 hydrophobicity + 2 charge + 2 α-helix + 2 TM)

### Top Features by RF Importance

The five most informative features were `max_tm_propensity` (~0.25), `max_hydrophobicity`, `avg_tm_propensity`, `avg_hydrophobicity`, and `comp_L` (~0.10). This is biologically intuitive — signal peptides are defined largely by their hydrophobic core — and suggests that a model built around hydrophobicity and TM propensity alone would recover most of the classification signal.

## Input Files (from previous steps)

| File | Description |
|------|-------------|
| `training_with_folds.tsv` | Training set with label and 5-fold assignment columns |
| `benchmarking_set.tsv` | Blind benchmark set with label column |
| `filtered_positive.tsv` | Clustered positive representatives (Step 2) |
| `filtered_negative.tsv` | Clustered negative representatives (Step 2) |
| `positive.fasta` / `negative.fasta` | Source protein sequences (filtered accessions only) |

> **Note:** All sequences are restricted to the filtered non-redundant set from Step 2, consistent with Steps 4 and 6.

## Results

### Hyperparameter Search

| Parameter | Value |
|-----------|-------|
| Kernel | RBF |
| C grid | [0.1, 1, 10] |
| γ grid | ['scale', 0.01, 0.001] |
| Selection metric | MCC (inner CV) |

MCC was chosen as the optimisation metric rather than accuracy because the class imbalance (~8:1) makes accuracy a misleading objective — a model predicting SP− for everything would score ~89% accuracy trivially.

### Performance

| Metric | CV (all 28) | Benchmark (all 28) | CV (sel 20) | Benchmark (sel 20) |
|--------|------------:|-------------------:|------------:|-------------------:|
| MCC | 0.852 | 0.854 | 0.853 | 0.846 |
| Precision | 0.867 | 0.872 | 0.867 | 0.866 |
| Recall | 0.870 | 0.868 | 0.872 | 0.858 |
| F1 | 0.868 | 0.870 | 0.869 | 0.862 |
| Accuracy | 0.971 | 0.972 | 0.971 | 0.970 |

Both configurations perform very similarly. The marginal edge of the full-feature model on the benchmark (MCC 0.854 vs 0.846) suggests that the 8 dropped features carry some signal, but not enough to justify the added complexity in most use cases. The near-perfect CV–benchmark agreement indicates the model is not overfitting.

### Benchmark Confusion Matrices

**All features (28)**

|  | Pred SP− | Pred SP+ |
|--|--------:|--------:|
| **True SP−** | 1759 (TN) | 28 (FP) |
| **True SP+** | 29 (FN) | 190 (TP) |

**Selected features (20)**

|  | Pred SP− | Pred SP+ |
|--|--------:|--------:|
| **True SP−** | 1758 (TN) | 29 (FP) |
| **True SP+** | 31 (FN) | 188 (TP) |

### Cross-Validation Curves (all features)

- **PR AUC:** 0.918 — OOF operating point: F1 = 0.868 (recall ≈ 0.87, precision ≈ 0.87)
- **ROC AUC:** 0.986

Compared to the Von Heijne baseline (F1 = 0.668, ROC-AUC = 0.954), the SVM delivers a substantial improvement simply by widening the feature scope beyond the cleavage-site window.

## Dependencies

```
biopython
pandas
numpy
scikit-learn
matplotlib
seaborn
```

Install with:
```bash
pip install biopython pandas numpy scikit-learn matplotlib seaborn
```

## Usage

Open and run `step5_svm.ipynb` cell by cell. Ensure the input `.tsv` and `.fasta` files from Steps 2 and 3 are available in the working directory. Output figures are saved to the `figures/` subdirectory.

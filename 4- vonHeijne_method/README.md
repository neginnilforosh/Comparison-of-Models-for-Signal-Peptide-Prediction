# Step 4 — Von Heijne Signal Peptide Classifier

**LB2 Project · Group 7 · Signal Peptide Prediction**

## Overview

Before reaching for machine learning, it makes sense to understand how well a classical, interpretable method performs on this dataset. The Von Heijne method is the canonical baseline for signal peptide detection — it uses a Position-Specific Weight Matrix (PSWM) built from known cleavage site windows to score any new sequence. Its simplicity is also its strength: every prediction can be directly traced back to position-specific amino acid preferences, and the PSWM itself is biologically interpretable.

## Contents

| File | Description |
|------|-------------|
| `step4_von_heijne.ipynb` | Main notebook — full pipeline from data loading to evaluation |
| `figures/vh_pswm_heatmap.pdf/.png` | PSWM log-odds weights heatmap (positions −13 to +2) |
| `figures/vh_cv_pr_roc.pdf/.png` | 5-fold cross-validation PR and ROC curves |
| `figures/vh_benchmark_confusion.pdf/.png` | Benchmark confusion matrix |

## Method

1. **Window extraction** — A 15-position window (positions **−13 to +2** relative to the cleavage site) is extracted from each positive training sequence. This window was chosen based on the sequence logo analysis in Step 3, where conserved positions were concentrated in this range.
2. **PSWM construction** — Log-odds scores are computed against SwissProt background amino acid frequencies, with pseudocount = 1.0 added to avoid log(0) for rare amino acids at any position.
3. **Scoring** — Each protein is scored by sliding the PSWM across N-terminal positions 15–100 and taking the maximum score. The upper bound of 100 is generous relative to the typical SP length (~22 aa) but avoids penalising proteins where the annotation does not perfectly match the retrieved sequence start.
4. **Threshold selection** — We use 5-fold cross-validation to select the classification threshold. For each fold, the threshold that maximises F1 on the out-of-fold (OOF) PR curve is recorded. The average OOF threshold (6.339) is then applied to the benchmark set — a clean separation between threshold selection and final evaluation.
5. **Final model** — Retrained on the full training set (873 positive sequences) and evaluated once on the blind benchmark.

## Input Files (from previous steps)

| File | Description |
|------|-------------|
| `filtered_positive.tsv` | Clustered positive representatives (Step 2) |
| `filtered_negative.tsv` | Clustered negative representatives (Step 2) |
| `training_with_folds.tsv` | Training set with label and 5-fold assignment columns |
| `benchmarking_set.tsv` | Blind benchmark set with label column |
| `positive.fasta` / `negative.fasta` | Source protein sequences (filtered accessions only) |

> **Note:** All sequences are restricted to the non-redundant set from Step 2. This is critical — evaluating on sequences that are highly similar to training sequences would artificially inflate performance.

## Results

### Model Parameters

| Parameter | Value |
|-----------|-------|
| Window size | 15 positions (−13 to +2) |
| Pseudocount | 1.0 |
| Background frequencies | SwissProt (ExPASy) |
| Cleavage site search range | positions 15–100 |
| Training positives | 873 |
| OOF threshold | 6.339 |

### Performance

| Metric | 5-fold CV | Benchmark |
|--------|----------:|----------:|
| Accuracy | 0.935 | 0.922 |
| Precision | 0.686 | 0.625 |
| Recall | 0.736 | 0.717 |
| F1 | 0.710 | 0.668 |
| MCC | 0.674 | 0.626 |

The CV–benchmark gap is modest and consistent, suggesting the threshold generalises reasonably well. The relatively low precision (0.625) reflects the difficulty of the problem for a pure positional scoring approach: the PSWM captures local cleavage-site patterns well (ROC-AUC 0.954) but cannot account for broader N-terminal context that distinguishes true signal peptides from superficially similar hydrophobic stretches.

### Benchmark Confusion Matrix

|  | Predicted SP− | Predicted SP+ |
|--|-------------:|-------------:|
| **Actual SP−** | 1693 (TN) | 94 (FP) |
| **Actual SP+** | 62 (FN) | 157 (TP) |

### Cross-Validation Curves

- **PR AUC:** 0.782
- **ROC AUC:** 0.954
- OOF operating point: F1 = 0.710 at recall ≈ 0.74, precision ≈ 0.69

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

Open and run `step4_von_heijne.ipynb` cell by cell. Ensure the input `.tsv` and `.fasta` files from Step 2 are available in the working directory. Output figures are saved to the `figures/` subdirectory.

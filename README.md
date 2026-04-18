# Signal Peptide Prediction — LB2 Project · Group 7

A reproducible machine learning pipeline for **eukaryotic signal peptide (SP) prediction**, built from scratch using UniProtKB data. The project implements and compares two core approaches — a classical rule-based method (Von Heijne) and an SVM classifier — with a deep learning extension (CNN + BiLSTM on ESM-2 embeddings) developed as an optional advancement. All methods are evaluated under identical cross-validation and blind benchmark conditions.

> **Course:** Laboratory of Bioinformatics 2 (LB2) · MSc Bioinformatics · University of Bologna  
> **Group 7**

---

## What This Project Does

Signal peptides are short N-terminal sequences that direct proteins into the secretory pathway. Predicting their presence is a fundamental problem in proteomics — it matters for annotating newly sequenced proteomes, designing recombinant proteins for secretion, and identifying candidate drug targets among secreted proteins.

This pipeline:
1. Collects a high-confidence labeled dataset from UniProtKB (Swiss-Prot), using experimental evidence filters to avoid training on computationally predicted annotations
2. Removes sequence redundancy using MMseqs2 clustering, so the train/benchmark split is biologically meaningful
3. Performs exploratory analysis to validate biological coherence and motivate modelling choices
4. Trains and evaluates three classifiers of increasing complexity:
   - **Von Heijne** — position-specific weight matrix (PSWM) baseline
   - **SVM** — biochemical features + Random Forest selection
   - **CNN + LSTM** — ESM-2 protein language model embeddings

---

## Repository Structure

```
.
├── 1- Data_Collection/        # Step 1 — UniProtKB retrieval and filtering
├── 2- Data_Preparation/       # Step 2 — MMseqs2 clustering, splits, CV folds
├── 3- Data_analysis/          # Step 3 — Exploratory analysis and visualizations
├── 4- vonHeijne_method/       # Step 4 — PSWM baseline classifier
├── 5- SVM-SPSelection/        # Step 5 — SVM with feature selection
└── 6- Deep_learning/          # Step 6 — CNN-LSTM on ESM-2 embeddings
```

Each folder contains its own `README.md` with detailed documentation, a Jupyter notebook, and all outputs (figures, TSVs, saved models).

---

## Folder Guide

### [`1- Data_Collection/`](./1-%20Data_Collection/README.md)
**Notebook:** `step1_data_collection.ipynb`

Retrieves eukaryotic proteins from UniProtKB via REST API and applies strict experimental evidence filters to produce two labeled sets. A deliberate choice was made here to require experimentally confirmed signal peptides (not just predicted ones) for the positive set, and experimentally confirmed subcellular localisation to non-secretory compartments for the negative set — this avoids any circularity with computational annotation.

| Set | Criteria | Final count |
|-----|----------|-------------|
| Positive (SP+) | Swiss-Prot, PE1, experimentally confirmed SP, cleavage site known, SP ≥ 14 aa | 2,932 |
| Negative (SP−) | Swiss-Prot, PE1, no SP annotation, localized to cytosol/nucleus/mitochondrion/plastid/peroxisome/membrane | 20,615 |

**Key outputs:** `positive.fasta`, `negative.fasta`, `positive.tsv`, `negative.tsv`

---

### [`2- Data_Preparation/`](./2-%20Data_Preparation/README.md)
**Notebook:** `step2_data_preparation.ipynb`

Removes redundant sequences and constructs train/benchmark splits with cross-validation fold assignments. The same splits are reused across all modelling steps (4, 5, 6) so that model comparisons reflect genuine performance differences, not lucky data partitions.

| Step | Tool/method | Detail |
|------|-------------|--------|
| Redundancy removal | MMseqs2 `easy-cluster` | 30% identity, 40% coverage, connected-component mode |
| Train/benchmark split | Stratified 80/20 | Per-class shuffle, `RANDOM_SEED=42` |
| CV fold assignment | Stratified 5-fold | Per-class, pre-assigned |

| Set | Positive | Negative | Total |
|-----|----------|----------|-------|
| After clustering | 1,093 | 8,934 | 10,027 |
| Training (80%) | ~874 | ~7,147 | 8,021 |
| Benchmark (20%) | ~219 | ~1,787 | 2,006 |

**Key outputs:** `training_with_folds.tsv`, `benchmarking_set.tsv`, `filtered_positive.tsv`, `filtered_negative.tsv`

---

### [`3- Data_analysis/`](./3-%20Data_analysis/README.md)
**Notebook:** `step3_data_analysis.ipynb`

Exploratory analysis performed on the full dataset before any model is trained. The goal is to confirm the dataset is biologically coherent and to inform downstream design choices — particularly, what N-terminal window lengths and feature types are likely to be informative.

| Figure | What it shows |
|--------|---------------|
| `01_protein_length` | SP+ proteins are shorter (median ~ 300 aa) than SP− (median ~ 450 aa) |
| `02_sp_length` | SP lengths tightly distributed around mean=22.9 aa, range 14–65 aa |
| `03_aa_composition` | SP regions enriched in Leu/Ala and depleted in charged residues vs SwissProt background |
| `04_kingdom_distribution` | Both sets dominated by Metazoa; different kingdom proportions noted as potential bias |
| `05_cleavage_site_logo` | Strong Ala conservation at −1 and −3 confirms the von Heijne −1/−3 rule |

All figures are exported as both PNG and PDF inside `3- Data_analysis/`.

---

### [`4- vonHeijne_method/`](./4-%20vonHeijne_method/README.md)
**Notebook:** `step4_von_heijne.ipynb`

Classical PSWM-based classifier serving as an interpretable rule-based baseline. Because every prediction can be traced back to position-specific amino acid frequencies, this method also acts as a sanity check on the dataset — if the PSWM heatmap did not show the expected hydrophobic core pattern, it would indicate a problem with data quality.

- **Window:** 15 positions (−13 to +2 relative to cleavage site)
- **Scoring:** log-odds vs SwissProt background, pseudocount=1.0
- **Inference:** sliding window over N-terminal positions 15–100, maximum score

| Metric | 5-Fold CV |
|--------|-----------|
| ROC-AUC | 0.954 |
| PR-AUC | 0.782 |
| F1 | 0.708 |

Benchmark confusion matrix: TP=157, FP=94, TN=1693, FN=62

---

### [`5- SVM-SPSelection/`](./5-%20SVM-SPSelection/README.md)
**Notebook:** `step5_svm.ipynb`

SVM classifier on 28 biochemical features from each protein's N-terminal region, with feature selection via Random Forest importance. The top features — hydrophobicity and TM propensity — confirm that the model is picking up on genuine signal peptide biology rather than dataset artefacts.

- **Features:** AA composition (first 20 aa), hydrophobicity, TM propensity, alpha-helix propensity, charge features (first 40 aa)
- **Selection:** top 20 features by mean RF importance across CV folds
- **SVM:** RBF kernel, grid search over C and γ, optimising MCC

Top features: `max_tm_propensity`, `max_hydrophobicity`, `avg_tm_propensity`, `avg_hydrophobicity`, `comp_L`

| Metric | 5-Fold CV (all 28) | Benchmark (all 28) | Benchmark (top 20) |
|--------|--------------------|--------------------|---------------------|
| ROC-AUC | 0.986 | — | — |
| PR-AUC | 0.918 | — | — |
| F1 | 0.868 | — | — |
| TP | — | 190 | 188 |
| FP | — | 28 | 29 |

---

### [`6- Deep_learning/`](./6-%20Deep_learning/README.md)
**Notebook:** `step6_deep_learning.ipynb`

CNN-LSTM model trained on ESM-2 protein language model embeddings — the best-performing method in the pipeline. The jump in performance over the SVM is largely attributable to the richer input representation: ESM-2 embeddings encode evolutionary and structural context that explicit biochemical features do not capture.

**Architecture:**
- Input: ESM-2 (`esm2_t12_35M_UR50D`, dim=480) embeddings of N-terminal 150 aa
- 3 × Conv1d blocks (64/128/128 filters, kernels 3/5/3) + BatchNorm + MaxPool
- 2-layer bidirectional LSTM (hidden=128)
- Classifier head with Dropout(0.3)

**Training:** Adam lr=1e-3, BCEWithLogitsLoss with class weights (pos_weight≈8.2), ReduceLROnPlateau, early stopping on val MCC (patience=5), gradient clipping (max norm=1.0)

| Metric | 5-Fold CV (mean ± std) | Benchmark |
|--------|------------------------|-----------|
| Precision | 0.976 ± 0.014 | 0.951 |
| Recall | 0.974 ± 0.009 | 0.973 |
| F1 | 0.975 ± 0.006 | 0.962 |
| MCC | 0.972 ± 0.007 | 0.957 |
| PR-AUC | 0.980 ± 0.009 | 0.984 |
| ROC-AUC | 0.996 ± 0.001 | 0.998 |

Benchmark confusion matrix: TP=213, FP=11, TN=1776, FN=6

**Saved model:** `deep_learning/cnn_signal_peptide_model.pt`

---

## Full Model Comparison (Blind Benchmark, n=2,006)

| Model | Precision | Recall | F1 | MCC | PR-AUC | ROC-AUC |
|-------|-----------|--------|----|-----|--------|---------|
| Von Heijne (PSWM) | 0.62 | 0.72 | 0.67 | 0.63 | 0.782 | 0.954 |
| SVM (20 features) | 0.87 | 0.86 | 0.87 | 0.85 | 0.918 | 0.986 |
| **CNN-LSTM (ESM-2)** | **0.951** | **0.973** | **0.962** | **0.957** | **0.984** | **0.998** |

All three methods were trained on identical data splits and evaluated on the same held-out benchmark, so the comparison is direct and fair.

---

## How to Reproduce

All notebooks are designed to run in **Google Colab**. Follow the steps in order:

1. Run `1- Data_Collection/DataCollection.ipynb` → produces `positive.fasta`, `negative.fasta`, `positive.tsv`, `negative.tsv`
2. Upload outputs to `2- Data_Preparation/step2_data_preparation.ipynb` → produces training and benchmark TSVs
3. Upload TSVs + FASTA to subsequent notebooks (steps 3–6); required inputs are listed in each folder's `README.md`

Each notebook installs its own dependencies via `pip` and `apt-get`.

> **Reproducibility:** All random operations use `RANDOM_SEED = 42`. Do not change this value across any step.

---

## Dependencies (summary)

| Tool | Used in |
|------|---------|
| `requests` | data_collection |
| `mmseqs2` | data_preparation |
| `biopython`, `pandas`, `numpy` | all steps |
| `matplotlib`, `seaborn`, `logomaker` | data_analysis |
| `scikit-learn` | von_heijne, SVM, deep_learning |
| `torch`, `fair-esm` | deep_learning |

Full dependency lists are in each folder's `README.md`.

---

## Data Source

- **UniProtKB / Swiss-Prot REST API:** <https://www.uniprot.org/help/api>
- **Swiss-Prot release statistics (ExPASy):** <https://web.expasy.org/docs/relnotes/relstat.html>
- **ESM-2 repository (Meta AI):** <https://github.com/facebookresearch/esm>

## Contact

For questions about this project, feel free to reach out:

- **Mahan Balooei** — mahan.balooei@studio.unibo.it
- **Negin Nillforoosh** — negin.nilforosh@studio.unibo.it
- **Kimia Kanouni** — kimia.kanouni@studio.unibo.it

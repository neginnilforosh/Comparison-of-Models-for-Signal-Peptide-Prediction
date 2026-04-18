# Step 3 — Data Analysis

**LB2 Project · Group 7 · Signal Peptide Prediction**

Before training any model, it is worth asking: does the dataset actually look like what we think it should look like biologically? This step performs exploratory data analysis (EDA) on the full labeled dataset (training + benchmark combined) to verify biological coherence, flag potential biases, and identify which sequence features are likely to be informative for classification. The analysis informs modelling choices in later steps — particularly the choice of N-terminal window length and feature engineering strategy.

---

## Contents

| File | Description |
|------|-------------|
| `step3_data_analysis.ipynb` | Full EDA pipeline with inline commentary |
| `01_protein_length.png/.pdf` | Protein length distribution — SP+ vs SP− |
| `02_sp_length.png/.pdf` | Signal peptide length distribution (SP+ only) |
| `03_aa_composition.png/.pdf` | Amino acid composition vs SwissProt background |
| `04_kingdom_distribution.png/.pdf` | Taxonomic kingdom distribution — SP+ vs SP− |
| `05_cleavage_site_logo.png/.pdf` | Sequence logo around the cleavage site (positions −13 to +2) |

---

## Analyses

### 1 — Protein length distribution

SP+ proteins tend to be shorter than SP− proteins (median ~300 aa vs ~450 aa). This is not surprising — cytosolic and nuclear proteins are on average larger, while secreted proteins tend to be more compact. The difference is large enough to be visible in density histograms but not so extreme as to create a trivially separable length-based classifier, which would have been a red flag.

### 2 — Signal peptide length distribution

Across all SP+ sequences, cleavage positions (i.e., SP lengths) are **tightly** clustered between 14 and 40 aa, with mean 22.9 aa and median 22 aa. This narrow distribution is reassuring — it matches the canonical hydrophobic core length needed for efficient translocon insertion, and it directly motivates truncating protein sequences to the first 40–150 aa in the modelling steps that follow.

### 3 — Amino acid composition vs SwissProt background

We compare per-residue frequencies in two regions: the SP region of SP+ sequences and the N-terminal 30 aa of SP− sequences, both normalised against the SwissProt proteome-wide background (ExPASy). The SP region is clearly enriched in hydrophobic residues (Leu, Ala) and depleted in charged and polar residues — exactly what the tripartite n-region/h-region/c-region model of signal peptides predicts. Seeing this pattern emerge from the data without any model fitting is a good sanity check that the annotation quality is high.

### 4 — Kingdom distribution

Both sets are dominated by Metazoa (SP+: 79.3%, SP−: 52.6%), but the proportions differ between classes. **This is worth flagging:** a model that learns kingdom-specific biases could partially exploit taxonomic identity rather than sequence content. We do not filter by kingdom here, but it is something to keep in mind when interpreting benchmark results, especially if evaluating on non-Metazoan proteins downstream.

### 5 — Cleavage site sequence logo

We extract a 15-residue window (positions −13 to +2 relative to the cleavage site) from all SP+ sequences with valid FASTA entries, build a position frequency matrix, and convert it to an information-content logo using `logomaker`. The result shows strong conservation of Ala at position −1 and Ala/Gly at position −3 — the classic **von Heijne −1/−3 rule** for signal peptidase I recognition. This validates both the annotation quality and the biological plausibility of the dataset, and directly motivates the 15-position window used in Step 4.

---

## Methods

| Library | Usage |
|---------|-------|
| `pandas`, `numpy` | Data loading and manipulation |
| `matplotlib`, `seaborn` | Histograms, boxplots, bar charts |
| `logomaker` | Sequence logo from information-content matrix |
| `biopython` (`SeqIO`) | FASTA parsing for sequence extraction |

All figures are saved in both PNG (for screen) and PDF (for publication) formats.

---

## Input Files

Upload the following files to Colab before running the notebook:

- `training_with_folds.tsv` — from [`2- Data_Preparation/`](../2-%20Data_Prepration/step2_data_preparation.ipynb)
- `benchmarking_set.tsv` — from [`2- Data_Preparation/`](../2-%20Data_Prepration/step2_data_preparation.ipynb)
- `positive.fasta` — from [`1- Data_Collection/`](../1-%20Data_Collection/DataCollection.ipynb)
- `negative.fasta` — from [`1- Data_Collection/`](../1-%20Data_Collection/DataCollection.ipynb)

---

## Dependencies

```
pip install biopython pandas numpy matplotlib seaborn logomaker
```

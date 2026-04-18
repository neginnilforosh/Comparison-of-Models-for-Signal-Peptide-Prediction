# Step 1 — Data Collection from UniProt

**LB2 Project · Group 7 · Signal Peptide Prediction**

This folder covers the first step of the pipeline: building a labeled dataset of signal peptide (SP) positive and negative protein sequences. Everything is retrieved programmatically from the UniProt REST API, with filtering baked directly into the query rather than applied post-hoc — this keeps the dataset clean without requiring an extra manual curation pass.

---

## Contents

| File                        | Description |
|-----------------------------|-------------|
| `step1_data_collection-2.ipynb` | End-to-end pipeline: UniProt querying, JSON parsing, filtering, and export to TSV + FASTA |
| `positive.tsv`              | Metadata for SP-positive sequences (2,932 entries) |
| `negative.tsv`              | Metadata for SP-negative sequences (20,615 entries) |
| `positive.fasta`            | SP-positive sequences in FASTA format |
| `negative.fasta`            | SP-negative sequences in FASTA format |

---

## Pipeline Overview

1. **Query UniProt REST API** — Rather than downloading a bulk dump and filtering locally, we use structured search URLs that apply all filters at query time. Pagination is handled automatically, with a `Retry` adapter to recover gracefully from transient network failures.
2. **Parse JSON responses** — Each entry is parsed to extract the accession, organism, taxonomic kingdom (Metazoa / Viridiplantae / Fungi / Other), sequence length, and SP-specific annotation fields.
3. **Filter entries** — Per-entry filter functions enforce the quality criteria described below. Only entries passing all checks are retained.
4. **Export** — Valid entries are written simultaneously to TSV (metadata) and FASTA (sequences), so both files are always in sync.

---

## Filtering Criteria

The filtering strategy is asymmetric by design: the two classes have different biological definitions, so they require different quality checks.

### Positive set (SP-present)

We restricted the positive set to cases where the signal peptide has been directly confirmed by experiment — not just predicted. This is stricter than many published datasets, but it avoids training the model on annotations that are themselves computational predictions, which would introduce circularity.

- Reviewed (Swiss-Prot) entries only
- Eukaryota (`taxonomy_id:2759`)
- Full-length proteins (no fragments)
- Sequence length ≥ 40 residues
- Protein existence at protein level (PE1)
- Experimentally confirmed signal peptide (`ft_signal_exp:*`)
- SP cleavage site confirmed; mature SP length > 13 residues

### Negative set (SP-absent)

For the negatives, it is not enough to simply exclude sequences with an SP annotation — some proteins are unannotated by default without actively being confirmed SP-absent. We therefore required experimental localization evidence placing the protein in a compartment incompatible with the secretory pathway.

- Reviewed (Swiss-Prot) entries only
- Eukaryota (`taxonomy_id:2759`)
- Full-length proteins (no fragments)
- Sequence length ≥ 40 residues
- No signal peptide annotation of any kind (`NOT ft_signal:*`)
- Experimentally confirmed localization to cytosol, nucleus, mitochondrion, plastid, peroxisome, or cell membrane

---

## Output Format

**TSV columns:**

| File            | Columns |
|-----------------|---------|
| `positive.tsv`  | `Accession`, `Organism`, `Kingdom`, `Sequence length`, `SP cleavage` |
| `negative.tsv`  | `Accession`, `Organism`, `Kingdom`, `Sequence length`, `N-term transmembrane` |

**FASTA headers:** Standard UniProt accession format (e.g. `>O00300`).

---

## Dataset Statistics (as of latest run)

| Set                  | Total Retrieved | After Filtering |
|----------------------|-----------------|-----------------|
| Positive (SP+)       | 2,949           | **2,932**       |
| Negative (SP−)       | 20,615          | **20,615**      |

**Note:** UniProt is updated regularly. These numbers reflect the exact run in `step1_data_collection-2.ipynb` and may shift slightly on a fresh run.

---

## How to Run

Open `step1_data_collection-2.ipynb` in Google Colab (or any Jupyter environment with internet access). All four output files will be generated automatically.

**Next step:** Run `data_preparation/step2_data_preparation.ipynb` to perform MMseqs2 clustering and the 80/20 train/benchmark split.

---

**Dependencies**
- Python 3.x
- `requests`
- `pandas`

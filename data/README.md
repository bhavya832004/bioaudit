# Data

This directory holds small local files used by BioAudit. Large raw datasets are not committed.

## Bundled files

| File | What it is | License / terms |
| --- | --- | --- |
| `sample_proteins.fasta` | Five synthetic peptides for an offline smoke test | MIT project code; sequences are not biological |
| `toy_classification.csv` | Created by `python scripts/prepare_dataset.py --source toy` | Synthetic labels and sequences; not a benchmark |

Gitignored (download locally; do not commit):

| File | What it is | Terms |
| --- | --- | --- |
| `deeploc_swissprot_train_validation.csv` | Official DeepLoc 2.0 train/validation CSV (~17 MB) | Cite DeepLoc 2.0 and UniProt |
| `deeploc_single_label.csv` | Processed single-location subset used by `configs/experiment_deeploc.yaml` | Same; multi-label proteins dropped |
| `uniprot_human_reviewed.fasta` | Optional UniProt FASTA | [UniProt license](https://www.uniprot.org/help/license) |

## How to obtain data

```bash
# Offline supervised toy fixture (default experiment)
python scripts/prepare_dataset.py --source toy --out data/toy_classification.csv

# DeepLoc 2.0 labeled localization set (real proteins; gitignored)
python scripts/prepare_dataset.py --source deeploc --limit 3000 --out data/deeploc_single_label.csv
bioaudit experiment configs/experiment_deeploc.yaml --out reports/experiment_deeploc

# Optional bounded UniProt FASTA (no classification labels)
python scripts/prepare_dataset.py --source uniprot --limit 1000 --out data/uniprot_human_reviewed.fasta
```

Each command writes a sibling `*.manifest.json` with source URL or generator name, retrieval time, record count, SHA-256, and preprocessing notes.

## Dataset selection for the main experiment

The supervised real-data path is **DeepLoc 2.0** (Swiss-Prot eukaryotic proteins with experimental subcellular location; homology partitions 0–4). BioAudit keeps proteins with **exactly one** of the ten location labels so the current logistic baseline is single-class. Multi-label proteins are excluded and documented.

Default `configs/experiment.yaml` still uses the synthetic toy fixture so tests stay offline. Use `configs/experiment_deeploc.yaml` after downloading.

## Citation

- UniProt Consortium. UniProt: the Universal Protein Knowledgebase. https://www.uniprot.org/ — follow UniProt attribution and [license](https://www.uniprot.org/help/license) if you download their data.

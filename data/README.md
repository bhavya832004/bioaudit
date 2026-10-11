# Data

This directory holds small local files used by BioAudit. Large raw datasets are not committed.

## Bundled files

| File | What it is | License / terms |
| --- | --- | --- |
| `sample_proteins.fasta` | Five synthetic peptides for an offline smoke test | MIT project code; sequences are not biological |
| `toy_classification.csv` | Created by `python scripts/prepare_dataset.py --source toy` | Synthetic labels and sequences; not a benchmark |

`uniprot_human_reviewed.fasta` is gitignored. Obtain it only if you accept UniProt terms.

## How to obtain data

```bash
# Offline supervised toy fixture (default experiment)
python scripts/prepare_dataset.py --source toy --out data/toy_classification.csv

# Optional bounded UniProt download (not the classification benchmark)
python scripts/prepare_dataset.py --source uniprot --limit 1000 --out data/uniprot_human_reviewed.fasta
```

Each command writes a sibling `*.manifest.json` with source URL or generator name, retrieval time, record count, SHA-256, and preprocessing notes.

## Dataset selection for the main experiment

Reviewed UniProt human proteins are useful for ingestion demos but **do not** provide a ready classification label, so they are not the supervised benchmark.

A suitable public task would be protein subcellular localization (for example DeepLoc / UniProt subcellular location annotations) with a documented license. That download is **not** performed by default and was not started in Phase 3. Until a licensed public set is prepared, the runnable experiment uses the synthetic toy fixture so split and model code stay testable offline.

## Citation

- UniProt Consortium. UniProt: the Universal Protein Knowledgebase. https://www.uniprot.org/ — follow UniProt attribution and [license](https://www.uniprot.org/help/license) if you download their data.

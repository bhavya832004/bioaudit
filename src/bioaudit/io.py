from pathlib import Path

import pandas as pd
from Bio import SeqIO


SUPPORTED_FASTA_SUFFIXES = {".fa", ".fasta", ".faa", ".fas"}
SUPPORTED_CSV_SUFFIXES = {".csv"}


def normalize_sequence(value: object) -> str:
    """Uppercase a sequence and drop all whitespace. Empty/missing values become ''."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return "".join(str(value).split()).upper()


def _require_file(path: str | Path) -> Path:
    resolved = Path(path)
    if not resolved.is_file():
        raise FileNotFoundError(f"Input file not found: {resolved}")
    return resolved


def read_fasta(path: str | Path) -> list[dict]:
    """Read FASTA into normalized records with id, description, and sequence."""
    path = _require_file(path)
    records = []
    try:
        for record in SeqIO.parse(str(path), "fasta"):
            records.append({
                "id": str(record.id),
                "description": str(record.description),
                "sequence": normalize_sequence(record.seq),
            })
    except ValueError as exc:
        raise ValueError(f"Could not parse FASTA file {path}: {exc}") from exc
    if not records:
        raise ValueError(
            f"No FASTA records found in: {path}. "
            "Expected one or more records starting with '>'."
        )
    return records


def read_csv_sequences(
    path: str | Path,
    sequence_column: str,
    label_column: str | None = None,
) -> list[dict]:
    path = _require_file(path)
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f"CSV file is empty: {path}") from exc
    if df.empty:
        raise ValueError(f"No rows found in: {path}")
    if sequence_column not in df.columns:
        raise ValueError(
            f"Sequence column {sequence_column!r} not found. "
            f"Available columns: {list(df.columns)}"
        )
    if label_column and label_column not in df.columns:
        raise ValueError(
            f"Label column {label_column!r} not found. "
            f"Available columns: {list(df.columns)}"
        )
    records = []
    for i, row in df.iterrows():
        seq = normalize_sequence(row[sequence_column])
        label = None
        if label_column:
            raw_label = row[label_column]
            if pd.isna(raw_label) or str(raw_label).strip() == "":
                label = None
            else:
                label = str(raw_label)
        if "id" in df.columns:
            record_id = row["id"]
        elif "ACC" in df.columns:
            record_id = row["ACC"]
        else:
            record_id = i
        if pd.isna(record_id):
            record_id = i
        record = {
            "id": str(record_id),
            "description": "",
            "sequence": seq,
            "label": label,
        }
        if "partition" in df.columns:
            raw_part = row["partition"]
            record["partition"] = None if pd.isna(raw_part) else str(raw_part)
        records.append(record)
    return records


def read_records(
    path: str | Path,
    sequence_column: str = "sequence",
    label_column: str | None = None,
) -> list[dict]:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in SUPPORTED_FASTA_SUFFIXES:
        return read_fasta(path)
    if suffix in SUPPORTED_CSV_SUFFIXES:
        return read_csv_sequences(path, sequence_column, label_column)
    raise ValueError(
        f"Unsupported input type {suffix!r} for {path}. "
        "Supported inputs: FASTA (.fa/.fasta/.faa/.fas) or CSV (.csv)"
    )

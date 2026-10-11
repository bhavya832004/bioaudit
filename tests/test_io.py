from pathlib import Path

import pytest

from bioaudit.io import normalize_sequence, read_csv_sequences, read_fasta, read_records


ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "data" / "sample_proteins.fasta"


def test_read_sample_fasta():
    records = read_fasta(SAMPLE)
    assert len(records) == 5
    assert records[0]["id"] == "P001"
    assert records[0]["sequence"].startswith("MKWVTF")
    assert all(r["sequence"] == r["sequence"].upper() for r in records)


def test_normalize_sequence_whitespace_and_lowercase():
    assert normalize_sequence(" mk wv\ttf\n") == "MKWVTF"
    assert normalize_sequence(None) == ""


def test_missing_file_error(tmp_path):
    missing = tmp_path / "absent.fasta"
    with pytest.raises(FileNotFoundError, match="Input file not found"):
        read_records(missing)


def test_empty_fasta_error(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="No FASTA records found"):
        read_fasta(path)


def test_malformed_fasta_without_records(tmp_path):
    path = tmp_path / "not_fasta.fasta"
    path.write_text("this is not fasta\nMKWVTF\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Could not parse FASTA file"):
        read_fasta(path)


def test_unsupported_extension(tmp_path):
    path = tmp_path / "proteins.txt"
    path.write_text(">A\nACDE\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported input type"):
        read_records(path)


def test_empty_csv_error(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="CSV file is empty"):
        read_csv_sequences(path, "sequence")


def test_csv_header_only_error(tmp_path):
    path = tmp_path / "header.csv"
    path.write_text("id,sequence\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No rows found"):
        read_csv_sequences(path, "sequence")


def test_csv_missing_sequence_column(tmp_path):
    path = tmp_path / "cols.csv"
    path.write_text("id,seq\nA,ACDE\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Sequence column 'sequence' not found"):
        read_csv_sequences(path, "sequence")


def test_csv_missing_label_column(tmp_path):
    path = tmp_path / "nolabel.csv"
    path.write_text("id,sequence\nA,ACDE\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Label column 'label' not found"):
        read_csv_sequences(path, "sequence", "label")


def test_csv_missing_values_and_generated_ids(tmp_path):
    path = tmp_path / "proteins.csv"
    path.write_text("sequence,label\nACDE,pos\n,neg\nGAVLI, \n", encoding="utf-8")
    records = read_csv_sequences(path, "sequence", "label")
    assert records[0]["id"] == "0"
    assert records[0]["sequence"] == "ACDE"
    assert records[1]["sequence"] == ""
    assert records[1]["label"] == "neg"
    assert len(records) == 3

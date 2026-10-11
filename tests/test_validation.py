from bioaudit.audit import audit_records
from bioaudit.io import read_csv_sequences


def _codes(result):
    return [item["code"] for item in result["warnings"]]


def test_lowercase_and_whitespace_via_csv(tmp_path):
    path = tmp_path / "mixed.csv"
    path.write_text("id,sequence\nA, mk wv tf \nB,acde\n", encoding="utf-8")
    records = read_csv_sequences(path, "sequence")
    result = audit_records(records)
    assert result["errors"] == []
    assert result["metrics"]["record_count"] == 2
    assert "unsupported_symbol" not in _codes(result)
    assert records[0]["sequence"] == "MKWVTF"


def test_gaps_and_stop_are_not_called_invalid():
    result = audit_records([
        {"id": "gapped", "sequence": "ACDE-FGH"},
        {"id": "stopped", "sequence": "ACDE*"},
    ])
    codes = _codes(result)
    assert "gap_symbol" in codes
    assert "stop_symbol" in codes
    assert "unsupported_symbol" not in codes
    assert result["invalid_symbol_records"] == []
    gap = next(w for w in result["warnings"] if w["code"] == "gap_symbol")
    assert gap["count"] == 1
    assert "gapped" in gap["ids"]


def test_ambiguous_and_extended_residues_warn_but_are_accepted():
    result = audit_records([
        {"id": "bzx", "sequence": "ACDEXBZ"},
        {"id": "jou", "sequence": "ACDEJOU"},
    ])
    assert "extended_residue" in _codes(result)
    assert "unsupported_symbol" not in _codes(result)
    assert result["record_count"] == 2


def test_duplicate_ids_and_empty_sequences():
    result = audit_records([
        {"id": "same", "sequence": "ACDE"},
        {"id": "same", "sequence": ""},
        {"id": "other", "sequence": "GAVL"},
    ])
    codes = _codes(result)
    assert "duplicate_record_ids" in codes
    assert "empty_sequence" in codes
    dup = next(w for w in result["warnings"] if w["code"] == "duplicate_record_ids")
    assert "same" in dup["ids"]
    assert result["errors"] == []


def test_empty_dataset_is_fatal_error():
    result = audit_records([])
    assert result["errors"]
    assert result["errors"][0]["code"] == "empty_dataset"
    assert result["metrics"]["records_read"] == 0
    assert result["metrics"]["filtering"] == "none"


def test_length_thresholds():
    result = audit_records(
        [
            {"id": "short", "sequence": "AC"},
            {"id": "ok", "sequence": "ACDEFG"},
            {"id": "long", "sequence": "ACDEFGHIKLMN"},
        ],
        min_length=4,
        max_length=8,
    )
    codes = _codes(result)
    assert "below_min_length" in codes
    assert "above_max_length" in codes


def test_csv_missing_labels(tmp_path):
    path = tmp_path / "labels.csv"
    path.write_text("id,sequence,label\nA,ACDE,pos\nB,GAVL,\n", encoding="utf-8")
    records = read_csv_sequences(path, "sequence", "label")
    result = audit_records(records)
    assert "missing_label" in _codes(result)
    assert result["label_counts"] == {"pos": 1}
    assert result["class_imbalance_ratio_max_to_min"] == 1.0


def test_findings_include_counts_and_schema_keys():
    result = audit_records([
        {"id": "a", "sequence": "ACDE"},
        {"id": "b", "sequence": "ACDE"},
    ])
    assert "errors" in result and "warnings" in result and "metrics" in result
    warning = next(w for w in result["warnings"] if w["code"] == "exact_duplicate_sequence")
    assert warning["count"] == 2
    assert "a" in warning["ids"] and "b" in warning["ids"]
    assert result["metrics"]["records_read"] == result["metrics"]["records_audited"]

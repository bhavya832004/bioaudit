from bioaudit.audit import audit_records
from bioaudit.similarity import jaccard_similarity, approximate_cross_split_hits


def test_duplicate_detection():
    records = [
        {"id": "a", "sequence": "MKWVTF"},
        {"id": "b", "sequence": "MKWVTF"},
        {"id": "c", "sequence": "GAVLIP"},
    ]
    result = audit_records(records)
    assert result["record_count"] == 3
    assert result["exact_duplicate_excess_records"] == 1
    assert result["unique_nonempty_sequences"] == 2
    assert result["exact_duplicate_groups"][0]["ids"] == ["a", "b"]


def test_empty_records_and_empty_sequences():
    result = audit_records([])
    assert result["record_count"] == 0
    assert result["sequence_length"]["min"] is None
    assert result["empty_sequences"] == []

    result = audit_records([
        {"id": "e", "sequence": ""},
        {"id": "ok", "sequence": "ACDEFG"},
    ])
    assert result["empty_sequences"] == ["e"]
    assert result["sequence_length"]["count"] == 1


def test_invalid_symbols():
    result = audit_records([{"id": "bad", "sequence": "MKW123"}])
    assert result["invalid_symbol_records"][0]["id"] == "bad"
    assert result["invalid_symbol_records"][0]["invalid_symbols"] == ["1", "2", "3"]
    assert any(w["code"] == "unsupported_symbol" for w in result["warnings"])


def test_label_imbalance():
    records = [
        {"id": "a", "sequence": "ACDE", "label": "pos"},
        {"id": "b", "sequence": "GAVL", "label": "pos"},
        {"id": "c", "sequence": "MNQP", "label": "neg"},
    ]
    result = audit_records(records)
    assert result["label_counts"] == {"pos": 2, "neg": 1}
    assert result["class_imbalance_ratio_max_to_min"] == 2.0


def test_similarity_identical_sequences():
    assert jaccard_similarity("MKWVTFIS", "MKWVTFIS") == 1.0


def test_cross_split_hit():
    train = [{"id": "train1", "sequence": "MKWVTFISLLFLFSSAYS"}]
    test = [{"id": "test1", "sequence": "MKWVTFISLLFLFSSAYS"}]
    hits = approximate_cross_split_hits(train, test, threshold=0.7, k=3)
    assert len(hits) == 1
    assert hits[0]["test_id"] == "test1"
    assert hits[0]["metric"] == "kmer_jaccard"

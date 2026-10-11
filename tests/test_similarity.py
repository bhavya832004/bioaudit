from pathlib import Path

import pytest

from bioaudit.similarity import (
    approximate_cross_split_hits,
    jaccard_similarity,
    kmer_split_audit,
    kmers,
    mmseqs_available,
    mmseqs_split_audit,
    parse_mmseqs_hits,
)

FIXTURE = Path(__file__).parent / "fixtures" / "mmseqs_hits.tsv"


def test_jaccard_identical_and_empty():
    assert jaccard_similarity("MKWVTFIS", "MKWVTFIS") == 1.0
    assert jaccard_similarity("", "") == 1.0
    assert jaccard_similarity("ACDE", "") == 0.0


def test_jaccard_dissimilar_sequences():
    score = jaccard_similarity("AAAAAA", "WWWWWW", k=3)
    assert score == 0.0


def test_short_sequence_kmers():
    assert kmers("AC", k=3) == {"AC"}
    assert kmers("", k=3) == set()


def test_invalid_k_rejected():
    with pytest.raises(ValueError, match="k-mer size"):
        kmers("ACDE", k=0)


def test_cross_split_no_train_records():
    hits = approximate_cross_split_hits([], [{"id": "t", "sequence": "ACDE"}], threshold=0.1)
    assert hits == []


def test_cross_split_threshold_validation():
    with pytest.raises(ValueError, match="threshold"):
        approximate_cross_split_hits(
            [{"id": "a", "sequence": "ACDE"}],
            [{"id": "b", "sequence": "ACDE"}],
            threshold=1.5,
        )


def test_parse_mmseqs_hits_valid_empty_and_malformed():
    text = "\n".join([
        "# comment",
        "Q1\tT1\t0.95\t100\t1e-20\t200",
        "truncated\tonly\tthree",
        "Q2\tT2\tnot_a_float\t10\t0.1\t5",
        "",
        "Q3\tT3\t0.40\t50\t0.01\t12",
    ])
    hits = parse_mmseqs_hits(text)
    assert len(hits) == 2
    assert hits[0]["test_id"] == "Q1"
    assert hits[0]["sequence_identity"] == 0.95
    assert hits[1]["test_id"] == "Q3"
    assert parse_mmseqs_hits("") == []


def test_parse_mmseqs_hits_with_coverage_fixture():
    hits = parse_mmseqs_hits(FIXTURE.read_text(encoding="utf-8"))
    assert len(hits) == 3
    assert hits[0]["query_coverage"] == 0.90
    assert hits[0]["target_coverage"] == 0.85
    assert hits[0]["sequence_identity"] == 0.95


def test_kmer_audit_best_matches_and_sensitivity():
    train = [
        {"id": "train_same", "sequence": "MKWVTFISLLFLFSSAYS"},
        {"id": "train_other", "sequence": "GAVLIPFYWSTCMNQDEKRH"},
    ]
    test = [
        {"id": "identical", "sequence": "MKWVTFISLLFLFSSAYS"},
        {"id": "short", "sequence": "MK"},
        {"id": "empty", "sequence": ""},
        {"id": "ambiguous", "sequence": "MKWVXXXXXX"},
    ]
    result = kmer_split_audit(train, test, threshold=0.7, k=3, extra_thresholds=[0.3, 0.7])
    assert result["metric"] == "kmer_jaccard"
    assert "percent identity" not in result["method"].lower() or "NOT percent identity" in result["method"]
    by_id = {row["test_id"]: row for row in result["best_matches"]}
    assert by_id["identical"]["above_threshold"] is True
    assert by_id["identical"]["score"] == 1.0
    assert by_id["empty"]["above_threshold"] is False
    assert len(result["threshold_sensitivity"]) == 2


def test_kmer_refuses_large_pairwise_scans():
    train = [{"id": f"t{i}", "sequence": "ACDE"} for i in range(20)]
    test = [{"id": f"q{i}", "sequence": "ACDE"} for i in range(20)]
    with pytest.raises(ValueError, match="max_pairs"):
        kmer_split_audit(train, test, max_pairs=10, allow_quadratic=False)
    result = kmer_split_audit(train, test, max_pairs=10, allow_quadratic=True, threshold=0.9)
    assert result["hits"]


@pytest.mark.skipif(not mmseqs_available(), reason="MMseqs2 executable not installed")
def test_mmseqs_integration_when_installed():
    train = [{"id": "train1", "sequence": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"}]
    test = [{"id": "test1", "sequence": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"}]
    result = mmseqs_split_audit(train, test, threshold=0.3)
    assert result["metric"] == "mmseqs_identity"
    assert result["mmseqs_version"]
    assert result["hits"]

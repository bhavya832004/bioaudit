"""Approximate k-mer Jaccard similarity. This is not percent sequence identity."""

from __future__ import annotations

from statistics import mean, median

DEFAULT_MAX_PAIRS = 250_000
KMER_LIMITATIONS = (
    "kmer_jaccard depends on k, sequence length, and residue composition. "
    "It is not an alignment, not percent sequence identity, and not a homology proof. "
    "The implementation is brute-force (quadratic in train×test size)."
)


def kmers(sequence: str, k: int = 3) -> set[str]:
    if k < 1:
        raise ValueError("k-mer size k must be >= 1")
    if len(sequence) < k:
        return {sequence} if sequence else set()
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def jaccard_similarity(a: str, b: str, k: int = 3) -> float:
    """Approximate similarity based on sets of k-mers; not percent identity."""
    ka, kb = kmers(a, k), kmers(b, k)
    if not ka and not kb:
        return 1.0
    if not ka or not kb:
        return 0.0
    return len(ka & kb) / len(ka | kb)


def _require_pair_budget(n_train: int, n_test: int, max_pairs: int, allow_quadratic: bool) -> int:
    n_pairs = n_train * n_test
    if n_pairs > max_pairs and not allow_quadratic:
        raise ValueError(
            f"k-mer brute-force would compare {n_pairs} train/test pairs, "
            f"exceeding max_pairs={max_pairs}. Reduce the input, raise --max-pairs, "
            "pass --allow-quadratic, or use --backend mmseqs."
        )
    return n_pairs


def kmer_split_audit(
    train_records: list[dict],
    test_records: list[dict],
    threshold: float = 0.7,
    k: int = 3,
    max_pairs: int = DEFAULT_MAX_PAIRS,
    allow_quadratic: bool = False,
    extra_thresholds: list[float] | None = None,
) -> dict:
    """Best k-mer Jaccard match for every test record, plus above-threshold hits."""
    if not 0 <= threshold <= 1:
        raise ValueError("k-mer Jaccard threshold must be in [0, 1]")
    n_pairs = _require_pair_budget(len(train_records), len(test_records), max_pairs, allow_quadratic)

    best_matches: list[dict] = []
    for test in test_records:
        best = None
        for train in train_records:
            score = jaccard_similarity(test["sequence"], train["sequence"], k=k)
            if best is None or score > best["score"]:
                best = {
                    "test_id": test["id"],
                    "train_id": train["id"],
                    "metric": "kmer_jaccard",
                    "score": round(score, 5),
                    "kmer_jaccard": round(score, 5),
                    "k": k,
                    "threshold": threshold,
                    "above_threshold": score >= threshold,
                }
        if best is None:
            best = {
                "test_id": test["id"],
                "train_id": None,
                "metric": "kmer_jaccard",
                "score": None,
                "kmer_jaccard": None,
                "k": k,
                "threshold": threshold,
                "above_threshold": False,
            }
        best_matches.append(best)

    hits = [m for m in best_matches if m["above_threshold"]]
    scores = [m["score"] for m in best_matches if m["score"] is not None]
    thresholds = sorted({threshold, *(extra_thresholds or [])})
    sensitivity = []
    for t in thresholds:
        if not 0 <= t <= 1:
            continue
        n_hit = sum(1 for m in best_matches if m["score"] is not None and m["score"] >= t)
        sensitivity.append({
            "threshold": t,
            "metric": "kmer_jaccard",
            "test_records_with_hits": n_hit,
            "test_fraction_with_hits": round(n_hit / len(test_records), 5) if test_records else 0.0,
        })

    warnings = [KMER_LIMITATIONS]
    if extra_thresholds:
        warnings.append("Threshold sensitivity reuses the same pairwise k-mer scores; it is not a new alignment.")

    return {
        "method": f"k-mer Jaccard (k={k}; approximate similarity, NOT percent identity)",
        "metric": "kmer_jaccard",
        "k": k,
        "threshold": threshold,
        "pair_comparisons": n_pairs,
        "max_pairs": max_pairs,
        "best_matches": best_matches,
        "hits": hits,
        "score_distribution": {
            "count": len(scores),
            "min": min(scores) if scores else None,
            "median": median(scores) if scores else None,
            "mean": round(mean(scores), 5) if scores else None,
            "max": max(scores) if scores else None,
        },
        "threshold_sensitivity": sensitivity,
        "warnings": warnings,
        "interpretation": (
            "A k-mer Jaccard hit is a risk indicator of train/test similarity, "
            "not proof of leakage or homology."
        ),
    }


def approximate_cross_split_hits(
    train_records: list[dict],
    test_records: list[dict],
    threshold: float = 0.7,
    k: int = 3,
    max_pairs: int = DEFAULT_MAX_PAIRS,
    allow_quadratic: bool = False,
) -> list[dict]:
    """Brute-force k-mer neighbour scan. Returns above-threshold best matches only."""
    result = kmer_split_audit(
        train_records,
        test_records,
        threshold=threshold,
        k=k,
        max_pairs=max_pairs,
        allow_quadratic=allow_quadratic,
    )
    return result["hits"]

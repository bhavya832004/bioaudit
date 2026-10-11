"""Protein sequence validation for BioAudit.

Alphabet policy (protein FASTA/CSV, this release):
- Canonical residues (no symbol warning): ACDEFGHIKLMNPQRSTVWY
- Extended residues (kept; warning ``extended_residue``): B, J, O, U, X, Z
  - B = Asx (D or N), Z = Glx (E or Q), J = Xle (I or L), X = unknown
  - U = selenocysteine, O = pyrrolysine
- Stop ``*`` (kept; warning ``stop_symbol``). Not treated as an "invalid" residue.
- Gap ``-`` (kept; warning ``gap_symbol``). Not treated as an "invalid" residue.
- Any other character (kept; warning ``unsupported_symbol``).

Records are never deleted, relabelled, or rebalanced. Filtering must be opt-in
and is not performed by ``audit_records``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from statistics import mean, median
from typing import Any

from .io import normalize_sequence

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
EXTENDED_AA = set("BJOUXZ")
STOP_SYMBOL = "*"
GAP_SYMBOL = "-"
ALLOWED_SYMBOLS = CANONICAL_AA | EXTENDED_AA | {STOP_SYMBOL, GAP_SYMBOL}

# Backward-compatible name: symbols that are not "unsupported".
VALID_AA = ALLOWED_SYMBOLS

MAX_IDS_IN_FINDING = 100


def _ids_payload(ids: list[str]) -> dict[str, Any]:
    unique = list(dict.fromkeys(ids))
    payload: dict[str, Any] = {"count": len(ids)}
    payload["ids"] = unique[:MAX_IDS_IN_FINDING]
    if len(unique) > MAX_IDS_IN_FINDING:
        payload["ids_truncated"] = True
        payload["id_count"] = len(unique)
    else:
        payload["ids_truncated"] = False
        payload["id_count"] = len(unique)
    return payload


def _finding(code: str, message: str, ids: list[str], **extra: Any) -> dict[str, Any]:
    item = {"code": code, "message": message, **_ids_payload(ids)}
    item.update(extra)
    return item


def audit_records(
    records: list[dict],
    *,
    min_length: int | None = None,
    max_length: int | None = None,
) -> dict:
    """Validate records and return errors, warnings, and metrics.

    Sequences are normalized (uppercase, whitespace removed) before checks.
    ``*`` and ``-`` are retained.
    """
    if min_length is not None and min_length < 0:
        raise ValueError("--min-length must be >= 0")
    if max_length is not None and max_length < 0:
        raise ValueError("--max-length must be >= 0")
    if min_length is not None and max_length is not None and min_length > max_length:
        raise ValueError("min-length cannot be greater than max-length")

    errors: list[dict] = []
    warnings: list[dict] = []

    if not records:
        errors.append(_finding(
            "empty_dataset",
            "No records were provided. This is a fatal input error for an audit.",
            [],
        ))

    normalized: list[tuple[str, str, Any]] = []
    for i, record in enumerate(records):
        rid = str(record.get("id", i))
        seq = normalize_sequence(record.get("sequence", ""))
        normalized.append((rid, seq, record.get("label")))

    ids = [rid for rid, _, _ in normalized]
    sequences = [seq for _, seq, _ in normalized]

    id_map: dict[str, list[int]] = defaultdict(list)
    for idx, rid in enumerate(ids):
        id_map[rid].append(idx)
    duplicate_id_groups = [
        {"id": rid, "count": len(idxs), "positions": idxs}
        for rid, idxs in id_map.items() if len(idxs) > 1
    ]
    if duplicate_id_groups:
        warnings.append(_finding(
            "duplicate_record_ids",
            "Record IDs are not unique. BioAudit does not assume FASTA/CSV IDs are globally unique.",
            [g["id"] for g in duplicate_id_groups],
            groups=duplicate_id_groups,
        ))

    empty_ids = [rid for rid, seq, _ in normalized if not seq]
    if empty_ids:
        warnings.append(_finding(
            "empty_sequence",
            "One or more records have an empty sequence after whitespace normalization.",
            empty_ids,
        ))

    exact_map: dict[str, list[str]] = defaultdict(list)
    for rid, seq, _ in normalized:
        if seq:
            exact_map[seq].append(rid)
    duplicate_groups = [
        {"sequence_length": len(seq), "ids": members, "count": len(members)}
        for seq, members in exact_map.items() if len(members) > 1
    ]
    duplicate_count = sum(g["count"] - 1 for g in duplicate_groups)
    if duplicate_groups:
        dup_ids = [rid for g in duplicate_groups for rid in g["ids"]]
        warnings.append(_finding(
            "exact_duplicate_sequence",
            "Exact duplicate sequences were found. Records were not removed.",
            dup_ids,
            groups=duplicate_groups,
        ))

    extended_ids: list[str] = []
    stop_ids: list[str] = []
    gap_ids: list[str] = []
    unsupported: list[dict] = []
    for rid, seq, _ in normalized:
        symbols = set(seq)
        if symbols & EXTENDED_AA:
            extended_ids.append(rid)
        if STOP_SYMBOL in symbols:
            stop_ids.append(rid)
        if GAP_SYMBOL in symbols:
            gap_ids.append(rid)
        bad = sorted(symbols - ALLOWED_SYMBOLS)
        if bad:
            unsupported.append({
                "id": rid,
                "invalid_symbols": bad,
                "unsupported_symbols": bad,
            })

    if extended_ids:
        warnings.append(_finding(
            "extended_residue",
            "Sequences contain extended amino-acid codes (B, J, O, U, X, and/or Z). "
            "These are accepted under the protein alphabet policy but are not canonical residues.",
            extended_ids,
        ))
    if stop_ids:
        warnings.append(_finding(
            "stop_symbol",
            "Sequences contain '*' (translation stop). The symbol is kept and is not labelled invalid.",
            stop_ids,
        ))
    if gap_ids:
        warnings.append(_finding(
            "gap_symbol",
            "Sequences contain '-' (gap). The symbol is kept and is not labelled invalid.",
            gap_ids,
        ))
    if unsupported:
        warnings.append(_finding(
            "unsupported_symbol",
            "Sequences contain symbols outside the documented protein alphabet "
            f"(canonical {''.join(sorted(CANONICAL_AA))}, extended {''.join(sorted(EXTENDED_AA))}, "
            "stop '*', gap '-'). Records were not removed.",
            [row["id"] for row in unsupported],
            records=unsupported,
        ))

    lengths = [len(seq) for seq in sequences if seq]
    short_ids: list[str] = []
    long_ids: list[str] = []
    if min_length is not None:
        short_ids = [rid for rid, seq, _ in normalized if seq and len(seq) < min_length]
        if short_ids:
            warnings.append(_finding(
                "below_min_length",
                f"Non-empty sequences shorter than min_length={min_length}.",
                short_ids,
            ))
    if max_length is not None:
        long_ids = [rid for rid, seq, _ in normalized if seq and len(seq) > max_length]
        if long_ids:
            warnings.append(_finding(
                "above_max_length",
                f"Sequences longer than max_length={max_length}.",
                long_ids,
            ))

    labels = [label for _, _, label in normalized if label is not None]
    label_counts = dict(Counter(labels))
    missing_label_ids = [
        rid for rid, _, label in normalized if label is None
    ]
    # Only warn about missing labels when at least one record has a label.
    if labels and missing_label_ids:
        warnings.append(_finding(
            "missing_label",
            "Some records are missing a label while others have labels. Missing labels were not imputed.",
            missing_label_ids,
        ))

    imbalance_ratio = None
    imbalance_note = (
        "Largest class count divided by smallest class count among labels that appear "
        "at least once. 1.0 means equal class sizes. Missing labels are excluded. "
        "This is a descriptive ratio, not a biological reliability score."
    )
    if label_counts:
        positive_counts = [v for v in label_counts.values() if v > 0]
        imbalance_ratio = round(max(positive_counts) / min(positive_counts), 4)

    metrics = {
        "records_read": len(records),
        "records_audited": len(records),
        "filtering": "none",
        "alphabet": {
            "name": "protein",
            "canonical": "".join(sorted(CANONICAL_AA)),
            "extended": "".join(sorted(EXTENDED_AA)),
            "stop": STOP_SYMBOL,
            "gap": GAP_SYMBOL,
        },
        "record_count": len(records),
        "unique_nonempty_sequences": len(exact_map),
        "empty_sequence_count": len(empty_ids),
        "duplicate_record_id_groups": len(duplicate_id_groups),
        "sequence_length": {
            "count": len(lengths),
            "min": min(lengths) if lengths else None,
            "median": median(lengths) if lengths else None,
            "mean": round(mean(lengths), 2) if lengths else None,
            "max": max(lengths) if lengths else None,
            "min_length_threshold": min_length,
            "max_length_threshold": max_length,
        },
        "sequence_lengths": lengths,
        "exact_duplicate_groups": duplicate_groups,
        "exact_duplicate_excess_records": duplicate_count,
        "label_counts": label_counts,
        "class_imbalance_ratio_max_to_min": imbalance_ratio,
        "class_imbalance_note": imbalance_note if label_counts else None,
    }

    # Backward-compatible keys used by existing tests and HTML cards.
    result = {
        "errors": errors,
        "warnings": warnings,
        "metrics": metrics,
        "record_count": metrics["record_count"],
        "unique_nonempty_sequences": metrics["unique_nonempty_sequences"],
        "empty_sequences": empty_ids,
        "invalid_symbol_records": unsupported,
        "sequence_length": metrics["sequence_length"],
        "sequence_lengths": lengths,
        "exact_duplicate_groups": duplicate_groups,
        "exact_duplicate_excess_records": duplicate_count,
        "label_counts": label_counts,
        "class_imbalance_ratio_max_to_min": imbalance_ratio,
    }
    return result

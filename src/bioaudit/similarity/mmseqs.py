"""MMseqs2 backend. Sequence identity is not a universal homology criterion."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

MMSEQS_LIMITATIONS = (
    "MMseqs2 sequence identity depends on alignment parameters, coverage mode, "
    "and the selected threshold. Identity alone is not a universal homology criterion "
    "and does not prove train/test leakage."
)


def mmseqs_available() -> bool:
    return shutil.which("mmseqs") is not None


def mmseqs_version() -> str | None:
    executable = shutil.which("mmseqs")
    if not executable:
        return None
    proc = subprocess.run(
        [executable, "version"],
        capture_output=True,
        text=True,
        check=False,
    )
    text = (proc.stdout or proc.stderr or "").strip()
    return text.splitlines()[0] if text else "unknown"


def parse_mmseqs_hits(text: str) -> list[dict]:
    """Parse MMseqs2 tabular output.

    Supported layouts (no header):
    - 6 columns: query, target, fident, alnlen, evalue, bits
    - 8 columns: query, target, fident, alnlen, qcov, tcov, evalue, bits
    Malformed rows are skipped.
    """
    hits = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        try:
            row = {
                "test_id": parts[0],
                "train_id": parts[1],
                "metric": "mmseqs_identity",
                "score": float(parts[2]),
                "sequence_identity": float(parts[2]),
                "alignment_length": int(parts[3]),
            }
            if len(parts) >= 8:
                row["query_coverage"] = float(parts[4])
                row["target_coverage"] = float(parts[5])
                row["evalue"] = float(parts[6])
                row["bits"] = float(parts[7])
            else:
                row["evalue"] = float(parts[4])
                row["bits"] = float(parts[5])
            hits.append(row)
        except ValueError:
            continue
    return hits


def _select_top_k(hits: list[dict], top_k: int | None) -> list[dict]:
    if top_k is None:
        return hits
    if top_k < 1:
        raise ValueError("--top-k must be >= 1")
    grouped: dict[str, list[dict]] = {}
    for hit in hits:
        grouped.setdefault(hit["test_id"], []).append(hit)
    selected = []
    for test_id, group in grouped.items():
        group.sort(key=lambda h: (-h["sequence_identity"], -h.get("bits", 0)))
        selected.extend(group[:top_k])
    return selected


def mmseqs_cross_split_hits(
    train_records: list[dict],
    test_records: list[dict],
    threshold: float = 0.3,
    top_k: int | None = None,
) -> list[dict]:
    """Use MMseqs2 easy-search; threshold is minimum sequence identity fraction."""
    result = mmseqs_split_audit(train_records, test_records, threshold=threshold, top_k=top_k)
    return result["hits"]


def mmseqs_split_audit(
    train_records: list[dict],
    test_records: list[dict],
    threshold: float = 0.3,
    top_k: int | None = None,
    coverage: float = 0.8,
    cov_mode: int = 0,
) -> dict:
    if not (0 < threshold <= 1):
        raise ValueError("MMseqs2 identity threshold must be in (0, 1]")
    if not (0 < coverage <= 1):
        raise ValueError("MMseqs2 coverage (-c) must be in (0, 1]")
    if not mmseqs_available():
        raise RuntimeError(
            "MMseqs2 executable 'mmseqs' was not found on PATH. "
            "Install MMseqs2 or use --backend kmer for small datasets."
        )
    version = mmseqs_version()
    with tempfile.TemporaryDirectory(prefix="bioaudit_") as tmp:
        tmp_path = Path(tmp)
        train_fa = tmp_path / "train.fasta"
        test_fa = tmp_path / "test.fasta"
        out_tsv = tmp_path / "hits.tsv"
        tmp_db = tmp_path / "tmp"
        SeqIO.write(
            [SeqRecord(Seq(r["sequence"]), id=str(r["id"]), description="") for r in train_records],
            train_fa,
            "fasta",
        )
        SeqIO.write(
            [SeqRecord(Seq(r["sequence"]), id=str(r["id"]), description="") for r in test_records],
            test_fa,
            "fasta",
        )
        cmd = [
            "mmseqs", "easy-search", str(test_fa), str(train_fa), str(out_tsv), str(tmp_db),
            "--min-seq-id", str(threshold),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode),
            "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(f"MMseqs2 failed:\n{proc.stderr}")
        raw = out_tsv.read_text(encoding="utf-8") if out_tsv.exists() else ""
        hits = _select_top_k(parse_mmseqs_hits(raw), top_k)

    best_by_test: dict[str, dict] = {}
    for hit in hits:
        current = best_by_test.get(hit["test_id"])
        if current is None or hit["sequence_identity"] > current["sequence_identity"]:
            best_by_test[hit["test_id"]] = hit

    best_matches = []
    for test in test_records:
        hit = best_by_test.get(str(test["id"]))
        if hit:
            best_matches.append({
                **hit,
                "threshold": threshold,
                "above_threshold": True,
            })
        else:
            best_matches.append({
                "test_id": test["id"],
                "train_id": None,
                "metric": "mmseqs_identity",
                "score": None,
                "sequence_identity": None,
                "threshold": threshold,
                "above_threshold": False,
            })

    return {
        "method": "MMseqs2 sequence identity",
        "metric": "mmseqs_identity",
        "threshold": threshold,
        "coverage": coverage,
        "cov_mode": cov_mode,
        "top_k": top_k,
        "mmseqs_version": version,
        "command": cmd,
        "hits": hits,
        "best_matches": best_matches,
        "warnings": [MMSEQS_LIMITATIONS],
        "interpretation": (
            "An MMseqs2 identity hit is a risk indicator of related sequences across "
            "splits, not proof of leakage."
        ),
    }

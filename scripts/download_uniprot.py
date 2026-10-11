"""Download a small reviewed human-protein FASTA subset from UniProtKB."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bioaudit.dataset import download_uniprot_fasta


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--out", default="data/uniprot_human_reviewed.fasta")
    parser.add_argument("--timeout", type=int, default=90)
    args = parser.parse_args(argv)
    try:
        path, manifest = download_uniprot_fasta(args.limit, args.out, timeout=args.timeout)
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Saved sequences to {path}")
    print(f"Manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

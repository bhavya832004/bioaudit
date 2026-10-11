"""Prepare BioAudit datasets with provenance manifests.

Default source is the local synthetic toy classification fixture (offline).
UniProt download is bounded and optional; it is not the supervised benchmark.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bioaudit.dataset import (
    download_deeploc_localization,
    download_uniprot_fasta,
    write_toy_classification_csv,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Prepare BioAudit datasets.")
    parser.add_argument("--source", choices=["toy", "uniprot", "deeploc"], default="toy")
    parser.add_argument("--limit", type=int, default=1000, help="Record cap for UniProt/DeepLoc downloads.")
    parser.add_argument("--timeout", type=int, default=90)
    parser.add_argument(
        "--out",
        default=None,
        help="Output path. Defaults depend on --source.",
    )
    args = parser.parse_args(argv)
    try:
        if args.source == "toy":
            out = Path(args.out or "data/toy_classification.csv")
            path, manifest = write_toy_classification_csv(out)
            print(f"Wrote {path}")
            print(f"Manifest: {manifest}")
            print("This toy file is not a scientific benchmark.")
        elif args.source == "deeploc":
            out = Path(args.out or "data/deeploc_single_label.csv")
            path, manifest = download_deeploc_localization(
                out, limit=args.limit, timeout=args.timeout
            )
            print(f"Wrote {path}")
            print(f"Manifest: {manifest}")
            print("Cite DeepLoc 2.0 and UniProt. Multi-label proteins were excluded.")
        else:
            out = Path(args.out or "data/uniprot_human_reviewed.fasta")
            path, manifest = download_uniprot_fasta(args.limit, out, timeout=args.timeout)
            print(f"Wrote {path}")
            print(f"Manifest: {manifest}")
            print("UniProt sequences are not used as the default supervised benchmark.")
    except (FileNotFoundError, ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

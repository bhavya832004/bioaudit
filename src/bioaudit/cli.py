import argparse
import json
import logging
import sys

from . import __version__
from .io import read_records
from .audit import audit_records
from .similarity import kmer_split_audit, mmseqs_split_audit
from .report import write_report
from .experiment import load_experiment_config, run_experiment

log = logging.getLogger("bioaudit")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioaudit",
        description="Audit biological sequence datasets.",
    )
    parser.add_argument("--version", action="version", version=f"bioaudit {__version__}")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_stdout",
        help="Print machine-readable output paths as JSON after a successful command.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    audit = sub.add_parser("audit", help="Audit one FASTA or CSV dataset.")
    audit.add_argument("input")
    audit.add_argument("--sequence-column", default="sequence")
    audit.add_argument("--label-column")
    audit.add_argument("--min-length", type=int, default=None, help="Warn if non-empty sequences are shorter than this.")
    audit.add_argument("--max-length", type=int, default=None, help="Warn if sequences are longer than this.")
    audit.add_argument("--out", default="reports/audit")
    audit.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=True)

    split = sub.add_parser("split-audit", help="Check train/test sequence similarity.")
    split.add_argument("--train", required=True)
    split.add_argument("--test", required=True)
    split.add_argument("--sequence-column", default="sequence")
    split.add_argument("--threshold", type=float, default=0.7)
    split.add_argument(
        "--thresholds",
        default=None,
        help="Optional extra k-mer thresholds as a comma-separated list, e.g. 0.3,0.5,0.7",
    )
    split.add_argument("--backend", choices=["kmer", "mmseqs"], default="kmer")
    split.add_argument("--k", type=int, default=3)
    split.add_argument("--max-pairs", type=int, default=250000)
    split.add_argument(
        "--allow-quadratic",
        action="store_true",
        help="Allow k-mer brute-force above --max-pairs.",
    )
    split.add_argument("--top-k", type=int, default=None, help="Keep at most this many MMseqs2 hits per test sequence.")
    split.add_argument("--out", default="reports/split_audit")
    split.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=True)

    exp = sub.add_parser("experiment", help="Run random vs homology-aware evaluation from a YAML config.")
    exp.add_argument("config")
    exp.add_argument("--out", default=None)
    exp.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=True)
    return parser


def _parse_extra_thresholds(raw: str | None) -> list[float]:
    if not raw:
        return []
    values = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        values.append(float(part))
    return values


def _exact_cross_split_ids(train: list[dict], test: list[dict]) -> list[str]:
    train_seqs = {r["sequence"] for r in train if r.get("sequence")}
    return [r["id"] for r in test if r.get("sequence") and r["sequence"] in train_seqs]


def _emit(paths, args: argparse.Namespace) -> None:
    print(f"JSON: {paths[0]}")
    print(f"HTML: {paths[1]}")
    if args.json_stdout:
        print(json.dumps({"json": str(paths[0]), "html": str(paths[1])}))


def _run(args: argparse.Namespace) -> None:
    overwrite = getattr(args, "overwrite", True)
    if args.command == "audit":
        records = read_records(args.input, args.sequence_column, args.label_column)
        report = {
            "tool": "BioAudit",
            "audit_type": "dataset_quality",
            "input": str(args.input),
            "summary": audit_records(
                records,
                min_length=args.min_length,
                max_length=args.max_length,
            ),
        }
        paths = write_report(report, args.out, overwrite=overwrite)
        print(f"Audited {len(records)} records")
        _emit(paths, args)
        return
    if args.command == "experiment":
        config = load_experiment_config(args.config)
        out = args.out or (config.get("output") or {}).get("prefix") or "reports/experiment"
        report = run_experiment(config, config_path=args.config)
        paths = write_report(report, out, overwrite=overwrite)
        print(f"Experiment complete. Mean macro-F1 difference (pp): {report.get('mean_macro_f1_difference_pp')}")
        _emit(paths, args)
        return

    train = read_records(args.train, args.sequence_column)
    test = read_records(args.test, args.sequence_column)
    extra = _parse_extra_thresholds(args.thresholds)
    if args.backend == "mmseqs":
        audit = mmseqs_split_audit(
            train,
            test,
            threshold=args.threshold,
            top_k=args.top_k,
        )
    else:
        audit = kmer_split_audit(
            train,
            test,
            threshold=args.threshold,
            k=args.k,
            max_pairs=args.max_pairs,
            allow_quadratic=args.allow_quadratic,
            extra_thresholds=extra,
        )
    hits = audit["hits"]
    hit_ids = {h["test_id"] for h in hits}
    exact_ids = _exact_cross_split_ids(train, test)
    report = {
        "tool": "BioAudit",
        "audit_type": "train_test_similarity",
        "train_file": str(args.train),
        "test_file": str(args.test),
        "method": audit["method"],
        "metric": audit["metric"],
        "threshold": args.threshold,
        "backend": args.backend,
        "parameters": {k: v for k, v in audit.items() if k not in {"hits", "best_matches"}},
        "summary": {
            "train_record_count": len(train),
            "test_record_count": len(test),
            "exact_duplicate_test_ids": exact_ids,
            "exact_duplicate_test_count": len(exact_ids),
            "cross_split_hit_count": len(hits),
            "test_records_with_hits": len(hit_ids),
            "test_fraction_with_hits": round(len(hit_ids) / len(test), 5) if test else 0.0,
            "score_distribution": audit.get("score_distribution"),
            "threshold_sensitivity": audit.get("threshold_sensitivity"),
        },
        "best_matches": audit.get("best_matches"),
        "hits": hits,
        "warnings": audit.get("warnings", []),
        "interpretation": audit.get("interpretation"),
    }
    paths = write_report(report, args.out, overwrite=overwrite)
    print(f"Found {len(hits)} cross-split hits using {audit['method']}")
    _emit(paths, args)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        logging.basicConfig(
            level=logging.DEBUG if args.verbose else logging.INFO,
            format="%(levelname)s %(name)s: %(message)s",
        )
        _run(args)
    except (FileNotFoundError, FileExistsError, ValueError, RuntimeError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

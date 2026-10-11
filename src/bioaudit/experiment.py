"""Random vs homology-aware evaluation using a lightweight k-mer + logistic baseline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline

from . import __version__
from .audit import audit_records
from .io import read_records
from .metrics import classification_metrics
from .provenance import environment_metadata, utc_now_iso
from .similarity.kmer import kmer_split_audit
from .splitting import cluster_records, homology_aware_split, random_stratified_split, split_summary


def load_experiment_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Experiment config must be a mapping: {path}")
    return payload


def _subset(records: list[dict], idxs: list[int]) -> list[dict]:
    return [records[i] for i in idxs]


def _pipeline(k: int, C: float, max_iter: int, seed: int) -> Pipeline:
    return Pipeline([
        ("vec", CountVectorizer(analyzer="char", ngram_range=(k, k))),
        (
            "clf",
            LogisticRegression(
                C=C,
                max_iter=max_iter,
                class_weight="balanced",
                random_state=seed,
            ),
        ),
    ])


def _evaluate_split(
    records: list[dict],
    parts: dict[str, list[int]],
    *,
    k: int,
    C: float,
    max_iter: int,
    seed: int,
    homology_threshold: float,
) -> dict[str, Any]:
    train = _subset(records, parts["train"])
    test = _subset(records, parts["test"])
    y_train = [r["label"] for r in train]
    y_test = [r["label"] for r in test]
    if len(set(y_train)) < 2:
        raise ValueError("Training split contains fewer than two classes; cannot fit logistic regression.")
    model = _pipeline(k, C, max_iter, seed)
    model.fit([r["sequence"] for r in train], y_train)
    y_pred = list(model.predict([r["sequence"] for r in test]))
    classes = list(model.named_steps["clf"].classes_)
    proba = model.predict_proba([r["sequence"] for r in test])
    metrics = classification_metrics(y_test, y_pred, y_proba=proba, labels=sorted(set(y_train) | set(y_test)))
    leakage = kmer_split_audit(train, test, threshold=homology_threshold, k=k, allow_quadratic=True)
    return {
        "split_sizes": split_summary(records, parts),
        "metrics": metrics,
        "leakage": {
            "metric": leakage["metric"],
            "threshold": homology_threshold,
            "test_records_with_hits": leakage and len(leakage["hits"]),
            "test_fraction_with_hits": round(len(leakage["hits"]) / len(test), 5) if test else 0.0,
            "score_distribution": leakage.get("score_distribution"),
        },
        "model_classes": classes,
        "val_count": len(parts.get("val") or []),
    }


def run_experiment(config: dict[str, Any], *, config_path: str | Path | None = None) -> dict[str, Any]:
    dataset = config["dataset"]
    split_cfg = config["split"]
    model_cfg = config["model"]
    seeds = list(config.get("seeds") or [split_cfg.get("seed", 0)])
    records = read_records(
        dataset["path"],
        sequence_column=dataset.get("sequence_column", "sequence"),
        label_column=dataset.get("label_column", "label"),
    )
    if any(r.get("label") is None for r in records):
        raise ValueError("Experiment records must all have labels.")
    quality = audit_records(records)
    k = int(model_cfg.get("k", 3))
    C = float(model_cfg.get("C", 1.0))
    max_iter = int(model_cfg.get("max_iter", 400))
    test_size = float(split_cfg.get("test_size", 0.3))
    val_size = float(split_cfg.get("val_size", 0.2))
    hom_threshold = float(split_cfg.get("homology_threshold", 0.7))
    hom_k = int(split_cfg.get("homology_k", k))
    clusters = cluster_records(records, threshold=hom_threshold, k=hom_k)

    per_seed = []
    for seed in seeds:
        random_parts = random_stratified_split(
            records, test_size=test_size, val_size=val_size, seed=int(seed)
        )
        homo_parts = homology_aware_split(
            records, clusters, test_size=test_size, val_size=val_size, seed=int(seed)
        )
        random_eval = _evaluate_split(
            records, random_parts, k=k, C=C, max_iter=max_iter, seed=int(seed), homology_threshold=hom_threshold
        )
        homo_eval = _evaluate_split(
            records, homo_parts, k=k, C=C, max_iter=max_iter, seed=int(seed), homology_threshold=hom_threshold
        )
        delta = None
        if random_eval["metrics"]["macro_f1"] is not None and homo_eval["metrics"]["macro_f1"] is not None:
            delta = round(
                (random_eval["metrics"]["macro_f1"] - homo_eval["metrics"]["macro_f1"]) * 100,
                4,
            )
        per_seed.append({
            "seed": int(seed),
            "random": random_eval,
            "homology_aware": {
                **homo_eval,
                "cluster_overlap": homo_parts.get("cluster_overlap"),
                "split_notes": homo_parts.get("notes"),
            },
            "macro_f1_difference_pp": delta,
        })

    deltas = [row["macro_f1_difference_pp"] for row in per_seed if row["macro_f1_difference_pp"] is not None]
    mean_delta = round(sum(deltas) / len(deltas), 4) if deltas else None
    return {
        "tool": "BioAudit",
        "audit_type": "experiment",
        "bioaudit_version": __version__,
        "generated_at_utc": utc_now_iso(),
        "config_path": str(config_path) if config_path else None,
        "config": config,
        "environment": environment_metadata(),
        "dataset_quality": {
            "record_count": quality["record_count"],
            "warnings": quality["warnings"],
            "errors": quality["errors"],
            "label_counts": quality["label_counts"],
        },
        "n_clusters": len(clusters),
        "cluster_sizes": sorted(len(c) for c in clusters),
        "seeds": per_seed,
        "mean_macro_f1_difference_pp": mean_delta,
        "interpretation": (
            "A higher random-split score than a homology-aware split suggests similarity-based "
            "shortcuts may inflate apparent performance. This does not prove prior published "
            "results are invalid. Toy fixtures must not be reported as biological benchmarks."
        ),
        "limitations": [
            "Baseline is character k-mer counts plus logistic regression, not a protein language model.",
            "Homology clusters use k-mer Jaccard connected components unless a later MMseqs2 clustering path is used.",
            "Validation partitions are held out and are not used for model selection in this baseline.",
        ],
    }

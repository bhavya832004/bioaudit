"""Classification metrics for BioAudit experiments."""

from __future__ import annotations

from typing import Any
import warnings

import numpy as np
from sklearn.metrics import (
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)


def classification_metrics(
    y_true: list[str],
    y_pred: list[str],
    *,
    y_proba: np.ndarray | None = None,
    labels: list[str] | None = None,
) -> dict[str, Any]:
    if labels is None:
        labels = sorted(set(y_true) | set(y_pred))
    report = classification_report(
        y_true,
        y_pred,
        labels=labels,
        output_dict=True,
        zero_division=0,
    )
    per_class = {}
    for label in labels:
        row = report.get(label) or {}
        per_class[label] = {
            "precision": row.get("precision"),
            "recall": row.get("recall"),
            "f1": row.get("f1-score"),
            "support": row.get("support"),
        }
    matrix = confusion_matrix(y_true, y_pred, labels=labels).tolist()
    auroc = None
    unique_true = sorted(set(y_true))
    if y_proba is not None and len(unique_true) == 2 and y_proba.ndim == 2 and y_proba.shape[1] >= 2:
        positive = unique_true[1]
        class_index = labels.index(positive) if positive in labels else 1
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                auroc = float(roc_auc_score(y_true, y_proba[:, class_index], labels=unique_true))
        except ValueError:
            auroc = None
    elif y_proba is not None and len(unique_true) > 2 and y_proba.shape[1] == len(labels):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                auroc = float(roc_auc_score(y_true, y_proba, multi_class="ovr", labels=labels, average="macro"))
        except ValueError:
            auroc = None
    return {
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", labels=labels, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "per_class": per_class,
        "confusion_matrix": {"labels": labels, "matrix": matrix},
        "auroc": auroc,
        "n_true": len(y_true),
        "missing_classes_in_pred": sorted(set(labels) - set(y_pred)),
        "missing_classes_in_true": sorted(set(labels) - set(y_true)),
    }

"""Train/validation/test splitting: random stratified vs homology-aware."""

from __future__ import annotations

import random
from collections import Counter, defaultdict
from typing import Any

from sklearn.model_selection import train_test_split

from .similarity.kmer import jaccard_similarity


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def cluster_records(
    records: list[dict],
    threshold: float = 0.7,
    k: int = 3,
) -> list[list[int]]:
    """Connected components on exact identity or k-mer Jaccard >= threshold."""
    n = len(records)
    uf = UnionFind(n)
    seqs = [r.get("sequence", "") for r in records]
    by_seq: dict[str, list[int]] = defaultdict(list)
    for i, seq in enumerate(seqs):
        if seq:
            by_seq[seq].append(i)
    for members in by_seq.values():
        for extra in members[1:]:
            uf.union(members[0], extra)
    for i in range(n):
        for j in range(i + 1, n):
            if not seqs[i] or not seqs[j]:
                continue
            if uf.find(i) == uf.find(j):
                continue
            if jaccard_similarity(seqs[i], seqs[j], k=k) >= threshold:
                uf.union(i, j)
    groups: dict[int, list[int]] = defaultdict(list)
    for i in range(n):
        groups[uf.find(i)].append(i)
    return [sorted(members) for members in groups.values()]


def random_stratified_split(
    records: list[dict],
    *,
    test_size: float,
    val_size: float,
    seed: int,
) -> dict[str, list[int]]:
    idx = list(range(len(records)))
    labels = [r.get("label") for r in records]
    stratify = labels if all(label is not None for label in labels) and len(set(labels)) > 1 else None
    try:
        train_val, test = train_test_split(
            idx, test_size=test_size, random_state=seed, stratify=stratify
        )
    except ValueError:
        train_val, test = train_test_split(idx, test_size=test_size, random_state=seed, stratify=None)
    remaining_labels = [labels[i] for i in train_val]
    stratify_val = (
        remaining_labels
        if remaining_labels and all(x is not None for x in remaining_labels) and len(set(remaining_labels)) > 1
        else None
    )
    relative_val = val_size / (1.0 - test_size) if 0 < test_size < 1 else val_size
    try:
        train, val = train_test_split(
            train_val, test_size=relative_val, random_state=seed, stratify=stratify_val
        )
    except ValueError:
        train, val = train_test_split(train_val, test_size=relative_val, random_state=seed, stratify=None)
    return {"train": sorted(train), "val": sorted(val), "test": sorted(test)}


def homology_aware_split(
    records: list[dict],
    clusters: list[list[int]],
    *,
    test_size: float,
    val_size: float,
    seed: int,
) -> dict[str, Any]:
    """Assign each cluster to exactly one partition."""
    rng = random.Random(seed)

    def majority_label(cluster: list[int]) -> str:
        counts = Counter(records[i].get("label") for i in cluster)
        return str(counts.most_common(1)[0][0])

    by_label: dict[str, list[list[int]]] = defaultdict(list)
    for cluster in clusters:
        by_label[majority_label(cluster)].append(cluster)
    for group in by_label.values():
        rng.shuffle(group)
    labels = list(by_label)
    rng.shuffle(labels)
    order: list[list[int]] = []
    while any(by_label[label] for label in labels):
        for label in labels:
            if by_label[label]:
                order.append(by_label[label].pop())
    n = len(records)
    target_test = int(round(n * test_size)) if n else 0
    target_val = int(round(n * val_size)) if n else 0
    parts: dict[str, list[int]] = {"train": [], "val": [], "test": []}
    notes: list[str] = []
    for cluster in order:
        if len(cluster) > target_test > 0:
            parts["train"].extend(cluster)
            notes.append(
                f"Cluster of size {len(cluster)} exceeds target test size {target_test}; assigned to train."
            )
            continue
        test_need = target_test - len(parts["test"])
        val_need = target_val - len(parts["val"])
        if test_need >= len(cluster) and test_need > 0:
            parts["test"].extend(cluster)
        elif val_need >= len(cluster) and val_need > 0:
            parts["val"].extend(cluster)
        else:
            parts["train"].extend(cluster)
    for key in parts:
        parts[key] = sorted(set(parts[key]))
    return {
        "train": parts["train"],
        "val": parts["val"],
        "test": parts["test"],
        "notes": notes,
        "cluster_overlap": cluster_partition_violations(clusters, parts),
    }


def cluster_partition_violations(clusters: list[list[int]], parts: dict[str, list[int]]) -> list[int]:
    membership = {}
    for name, idxs in parts.items():
        if name not in {"train", "val", "test"}:
            continue
        for i in idxs:
            membership[i] = name
    bad = []
    for cluster in clusters:
        places = {membership.get(i) for i in cluster}
        if len(places) > 1:
            bad.append(min(cluster))
    return bad


def split_summary(records: list[dict], parts: dict[str, list[int]]) -> dict[str, Any]:
    summary = {}
    for name in ("train", "val", "test"):
        idxs = parts.get(name, [])
        labels = [records[i].get("label") for i in idxs]
        summary[name] = {
            "count": len(idxs),
            "label_counts": dict(Counter(label for label in labels if label is not None)),
        }
    return summary

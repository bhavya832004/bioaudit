from bioaudit.dataset import generate_toy_classification_records
from bioaudit.splitting import (
    cluster_by_field,
    cluster_partition_violations,
    cluster_records,
    homology_aware_split,
    random_stratified_split,
)


def test_clusters_keep_near_identical_families_together():
    records = generate_toy_classification_records()
    clusters = cluster_records(records, threshold=0.7, k=3)
    family_indexes = {}
    for i, rec in enumerate(records):
        family_indexes.setdefault(rec["family"], []).append(i)
    for members in family_indexes.values():
        roots = {next(ci for ci, cluster in enumerate(clusters) if members[0] in cluster)}
        for idx in members[1:]:
            found = next(ci for ci, cluster in enumerate(clusters) if idx in cluster)
            roots.add(found)
        assert len(roots) == 1, members


def test_homology_split_does_not_split_clusters():
    records = generate_toy_classification_records()
    clusters = cluster_records(records, threshold=0.7, k=3)
    parts = homology_aware_split(records, clusters, test_size=0.3, val_size=0.2, seed=0)
    assert cluster_partition_violations(clusters, parts) == []
    assigned = set(parts["train"]) | set(parts["val"]) | set(parts["test"])
    assert assigned == set(range(len(records)))


def test_cluster_by_partition_field():
    records = [
        {"id": "a", "partition": "0"},
        {"id": "b", "partition": "0"},
        {"id": "c", "partition": "1"},
    ]
    clusters = cluster_by_field(records, "partition")
    sizes = sorted(len(c) for c in clusters)
    assert sizes == [1, 2]


def test_random_split_sizes():
    records = generate_toy_classification_records()
    parts = random_stratified_split(records, test_size=0.3, val_size=0.2, seed=0)
    n = len(records)
    assert abs(len(parts["test"]) - round(n * 0.3)) <= 1
    assert set(parts["train"]).isdisjoint(parts["test"])

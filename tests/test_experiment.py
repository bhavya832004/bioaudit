import json
from pathlib import Path

import yaml

from bioaudit.cli import main
from bioaudit.dataset import write_toy_classification_csv
from bioaudit.experiment import load_experiment_config, run_experiment


def _config(csv_path: Path, tmp_path: Path) -> dict:
    return {
        "dataset": {
            "path": str(csv_path),
            "sequence_column": "sequence",
            "label_column": "label",
            "source": "toy",
        },
        "split": {
            "test_size": 0.3,
            "val_size": 0.2,
            "homology_threshold": 0.7,
            "homology_k": 3,
        },
        "model": {"name": "kmer_logistic_regression", "k": 3, "C": 1.0, "max_iter": 400},
        "seeds": [0, 1],
        "output": {"prefix": str(tmp_path / "exp")},
    }


def test_experiment_runs_on_toy_fixture(tmp_path):
    csv_path, _ = write_toy_classification_csv(tmp_path / "toy.csv")
    result = run_experiment(_config(csv_path, tmp_path), config_path="memory")
    assert result["audit_type"] == "experiment"
    assert result["n_clusters"] >= 1
    assert len(result["seeds"]) == 2
    for row in result["seeds"]:
        assert row["homology_aware"]["cluster_overlap"] == []
        assert "macro_f1" in row["random"]["metrics"]
        assert "macro_f1" in row["homology_aware"]["metrics"]
    assert result["mean_macro_f1_difference_pp"] is None or isinstance(result["mean_macro_f1_difference_pp"], float)


def test_experiment_cli(tmp_path):
    csv_path, _ = write_toy_classification_csv(tmp_path / "toy.csv")
    cfg = _config(csv_path, tmp_path)
    cfg_path = tmp_path / "exp.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    loaded = load_experiment_config(cfg_path)
    assert loaded["dataset"]["source"] == "toy"
    code = main(["experiment", str(cfg_path), "--out", str(tmp_path / "out")])
    assert code == 0
    payload = json.loads((tmp_path / "out.json").read_text(encoding="utf-8"))
    assert payload["audit_type"] == "experiment"
    html = (tmp_path / "out.html").read_text(encoding="utf-8")
    assert "<svg" in html

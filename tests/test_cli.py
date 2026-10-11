import json
from pathlib import Path

import pytest

from bioaudit.cli import main


ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "data" / "sample_proteins.fasta"


def test_cli_audit_sample_writes_json_and_html(tmp_path):
    out = tmp_path / "sample"
    code = main(["audit", str(SAMPLE), "--out", str(out)])
    assert code == 0
    json_path = tmp_path / "sample.json"
    html_path = tmp_path / "sample.html"
    assert json_path.is_file()
    assert html_path.is_file()
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["audit_type"] == "dataset_quality"
    assert payload["summary"]["record_count"] == 5
    assert payload["summary"]["exact_duplicate_excess_records"] == 1


def test_cli_split_audit_kmer(tmp_path):
    out = tmp_path / "split"
    code = main([
        "split-audit",
        "--train", str(SAMPLE),
        "--test", str(SAMPLE),
        "--threshold", "0.7",
        "--out", str(out),
    ])
    assert code == 0
    payload = json.loads((tmp_path / "split.json").read_text(encoding="utf-8"))
    assert payload["audit_type"] == "train_test_similarity"
    assert "NOT percent identity" in payload["method"]
    assert payload["summary"]["cross_split_hit_count"] >= 1


def test_cli_missing_file_is_actionable(tmp_path, capsys):
    missing = tmp_path / "nope.fasta"
    code = main(["audit", str(missing), "--out", str(tmp_path / "out")])
    assert code == 1
    err = capsys.readouterr().err
    assert "error:" in err
    assert "not found" in err.lower()


def test_cli_malformed_input_is_actionable(tmp_path, capsys):
    bad = tmp_path / "bad.fasta"
    bad.write_text("not a fasta file", encoding="utf-8")
    code = main(["audit", str(bad), "--out", str(tmp_path / "out")])
    assert code == 1
    err = capsys.readouterr().err
    assert "error:" in err
    assert "FASTA" in err


def test_cli_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "bioaudit" in capsys.readouterr().out.lower()


def test_cli_json_stdout_and_no_overwrite(tmp_path, capsys):
    out = tmp_path / "sample"
    assert main(["--json", "audit", str(SAMPLE), "--out", str(out)]) == 0
    stdout = capsys.readouterr().out
    assert '"json"' in stdout
    code = main(["audit", str(SAMPLE), "--out", str(out), "--no-overwrite"])
    assert code == 1
    assert "overwrite" in capsys.readouterr().err.lower()

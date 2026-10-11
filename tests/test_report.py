import json

import pytest

from bioaudit.report import write_report


def test_write_report_json_and_html_escape(tmp_path):
    out = tmp_path / "nested" / "report"
    payload = {
        "tool": "BioAudit",
        "input": "<script>alert('xss')</script>",
        "summary": {
            "record_count": 1,
            "unique_nonempty_sequences": 1,
            "exact_duplicate_excess_records": 0,
        },
        "note": "<img src=x onerror=alert(1)>",
    }
    json_path, html_path = write_report(payload, out)
    assert json_path == tmp_path / "nested" / "report.json"
    assert html_path == tmp_path / "nested" / "report.html"
    assert json_path.is_file()
    assert html_path.is_file()

    loaded = json.loads(json_path.read_text(encoding="utf-8"))
    assert loaded["input"] == payload["input"]
    assert "generated_at_utc" in loaded

    html = html_path.read_text(encoding="utf-8")
    assert "<script>alert('xss')</script>" not in html
    assert "<img src=x onerror=alert(1)>" not in html
    assert "&lt;script&gt;" in html
    assert "Record Count" in html
    assert html.strip().startswith("<!doctype html>")


def test_write_report_refuses_overwrite(tmp_path):
    out = tmp_path / "report"
    write_report({"summary": {"record_count": 1}}, out)
    with pytest.raises(FileExistsError, match="overwrite"):
        write_report({"summary": {"record_count": 2}}, out, overwrite=False)

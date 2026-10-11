from __future__ import annotations

from io import StringIO
from pathlib import Path
import html
import json
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _findings_section(title: str, items: list) -> str:
    if not items:
        return ""
    rows = []
    for item in items:
        if not isinstance(item, dict):
            rows.append(f"<li>{html.escape(str(item))}</li>")
            continue
        code = html.escape(str(item.get("code", "")))
        message = html.escape(str(item.get("message", "")))
        count = html.escape(str(item.get("count", "")))
        ids = item.get("ids") or []
        shown = ", ".join(html.escape(str(i)) for i in ids[:20])
        extra = f" IDs: {shown}" if shown else ""
        rows.append(f"<li><strong>{code}</strong> (count={count}): {message}{extra}</li>")
    return f"<h2>{html.escape(title)}</h2><ul class=\"findings\">{''.join(rows)}</ul>"


def _svg_from_axes(fig) -> str:
    buf = StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    text = buf.getvalue()
    start = text.find("<svg")
    return text[start:] if start >= 0 else text


def _bar_svg(labels: list[str], values: list[float], title: str) -> str:
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    ax.bar(range(len(values)), values, color="#1f6feb")
    ax.set_xticks(range(len(labels)), labels, rotation=20, ha="right")
    ax.set_title(title)
    ax.set_ylabel("value")
    fig.tight_layout()
    return _svg_from_axes(fig)


def _hist_svg(values: list[float], title: str) -> str:
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    ax.hist(values, bins=min(12, max(3, len(set(values)))), color="#1f6feb", edgecolor="white")
    ax.set_title(title)
    ax.set_xlabel("value")
    ax.set_ylabel("count")
    fig.tight_layout()
    return _svg_from_axes(fig)


def _charts(report: dict) -> str:
    blocks = []
    summary = report.get("summary") or {}
    lengths = summary.get("sequence_lengths") or (summary.get("metrics") or {}).get("sequence_lengths")
    if lengths:
        blocks.append(_hist_svg([float(x) for x in lengths], "Sequence length"))
    labels = summary.get("label_counts") or (summary.get("metrics") or {}).get("label_counts") or {}
    if labels:
        keys = [str(k) for k in labels]
        blocks.append(_bar_svg(keys, [float(labels[k]) for k in labels], "Label distribution"))
    matches = report.get("best_matches") or []
    scores = [m["score"] for m in matches if isinstance(m, dict) and m.get("score") is not None]
    if scores:
        blocks.append(_hist_svg(scores, "Best-match similarity scores"))
    sensitivity = (summary.get("threshold_sensitivity") or report.get("threshold_sensitivity") or [])
    if sensitivity:
        blocks.append(_bar_svg(
            [str(row.get("threshold")) for row in sensitivity],
            [float(row.get("test_fraction_with_hits") or 0) for row in sensitivity],
            "Threshold sensitivity (fraction of test with hits)",
        ))
    seeds = report.get("seeds") or []
    if seeds:
        random_f1 = [float(s["random"]["metrics"]["macro_f1"]) for s in seeds]
        homo_f1 = [float(s["homology_aware"]["metrics"]["macro_f1"]) for s in seeds]
        labels_s = [str(s["seed"]) for s in seeds]
        fig, ax = plt.subplots(figsize=(6.2, 3.2))
        x = range(len(labels_s))
        ax.bar([i - 0.18 for i in x], random_f1, width=0.36, label="random", color="#1f6feb")
        ax.bar([i + 0.18 for i in x], homo_f1, width=0.36, label="homology-aware", color="#9a3412")
        ax.set_xticks(list(x), labels_s)
        ax.set_title("Macro-F1 by split strategy")
        ax.legend()
        fig.tight_layout()
        blocks.append(_svg_from_axes(fig))
        matrix = seeds[0]["random"]["metrics"].get("confusion_matrix")
        if matrix and matrix.get("matrix"):
            fig, ax = plt.subplots(figsize=(4.8, 4.0))
            ax.imshow(matrix["matrix"], cmap="Blues")
            ax.set_xticks(range(len(matrix["labels"])), matrix["labels"], rotation=20, ha="right")
            ax.set_yticks(range(len(matrix["labels"])), matrix["labels"])
            ax.set_title("Confusion matrix (random split, seed 0)")
            fig.tight_layout()
            blocks.append(_svg_from_axes(fig))
    if not blocks:
        return ""
    return "<h2>Charts</h2>" + "".join(f'<div class="chart">{svg}</div>' for svg in blocks)


def write_report(
    report: dict,
    out_prefix: str | Path,
    *,
    overwrite: bool = True,
) -> tuple[Path, Path]:
    prefix = Path(out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_suffix(".json")
    html_path = prefix.with_suffix(".html")
    if not overwrite and (json_path.exists() or html_path.exists()):
        raise FileExistsError(
            f"Refusing to overwrite existing report {json_path} / {html_path}. Pass --overwrite."
        )
    report = dict(report)
    report["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    json_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

    summary = report.get("summary", report)
    cards = []
    for key in (
        "record_count", "unique_nonempty_sequences", "exact_duplicate_excess_records",
        "cross_split_hit_count", "exact_duplicate_test_count", "test_records_with_hits",
        "mean_macro_f1_difference_pp", "n_clusters",
    ):
        if key in summary or key in report:
            value = summary.get(key, report.get(key))
            cards.append(
                f'<div class="card"><div class="key">{html.escape(key.replace("_", " ").title())}'
                f'</div><div class="value">{html.escape(str(value))}</div></div>'
            )
    payload = html.escape(json.dumps(report, indent=2, default=str))
    warning_items = summary.get("warnings") or report.get("warnings") or []
    exec_bits = []
    errors = summary.get("errors") or report.get("dataset_quality", {}).get("errors") or []
    if errors:
        exec_bits.append(f"{len(errors)} error(s) in input validation.")
    if warning_items:
        exec_bits.append(f"{len(warning_items)} warning(s).")
    if report.get("mean_macro_f1_difference_pp") is not None:
        exec_bits.append(
            f"Mean random minus homology-aware macro-F1: {report['mean_macro_f1_difference_pp']} percentage points "
            "(toy/fixture results are not a biological benchmark)."
        )
    executive = ""
    if exec_bits:
        executive = "<h2>Executive summary</h2><ul>" + "".join(
            f"<li>{html.escape(bit)}</li>" for bit in exec_bits
        ) + "</ul>"
    findings = _findings_section("Errors", summary.get("errors") or []) + _findings_section(
        "Warnings", warning_items
    )
    charts = _charts(report)
    interpretation = report.get("interpretation") or ""
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>BioAudit report</title>
<style>
body{{font:16px/1.5 system-ui,sans-serif;max-width:1050px;margin:36px auto;padding:0 18px;color:#17202a}}
h1{{margin-bottom:4px}} .muted{{color:#667085}} .cards{{display:flex;gap:12px;flex-wrap:wrap;margin:24px 0}}
.card{{border:1px solid #d0d5dd;border-radius:10px;padding:16px;min-width:170px}}
.key{{font-size:13px;color:#667085}} .value{{font-size:25px;font-weight:700}}
pre{{background:#f6f8fa;padding:18px;border-radius:10px;overflow:auto;white-space:pre-wrap}}
ul.findings{{padding-left:20px}} .chart{{margin:16px 0;overflow:auto}} svg{{max-width:100%;height:auto}}
@media (max-width:640px){{body{{margin:16px auto}} .value{{font-size:20px}}}}
</style></head><body>
<h1>BioAudit report</h1><div class="muted">Generated {html.escape(str(report.get("generated_at_utc","")))}</div>
{executive}
<div class="cards">{''.join(cards)}</div>
{findings}
{charts}
<p class="muted">{html.escape(str(interpretation))}</p>
<p><a href="{html.escape(json_path.name)}">Download JSON results</a> (same folder as this HTML file).</p>
<h2>Full report data</h2><pre>{payload}</pre>
<p class="muted">Similarity estimates are method-dependent. Review methodology before drawing scientific conclusions. Unusual amino-acid codes are classified using the documented protein alphabet policy; they are not all labelled invalid. Charts are generated offline from this report's data.</p>
</body></html>"""
    html_path.write_text(document, encoding="utf-8")
    return json_path, html_path

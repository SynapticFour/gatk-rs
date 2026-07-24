#!/usr/bin/env python3
"""Publish nightly trio equivalence results: JSON summary, dashboard MD/HTML,
baseline comparison, and regression flag file for the workflow.
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
from datetime import datetime, timezone
from html import escape
from typing import Any


def parse_happy_summary(path: pathlib.Path) -> dict[str, dict[str, float]]:
    """Return Type → {precision, recall, f1, truth_tp, truth_fn, query_fp} for ALL/PASS rows."""
    out: dict[str, dict[str, float]] = {}
    if not path.is_file():
        return out
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            typ = (row.get("Type") or "").strip().upper()
            filt = (row.get("Filter") or "ALL").strip().upper()
            # Prefer FILTER=ALL aggregate rows when present
            if filt not in ("ALL", "", "PASS"):
                continue
            if typ not in ("SNP", "INDEL", "ROWS"):
                # Some builds use " contrived" — keep SNP/INDEL only
                if typ not in ("SNP", "INDEL"):
                    continue
            key = typ
            # Prefer ALL over PASS when both exist
            if key in out and filt == "PASS":
                continue

            def fget(*names: str, default: float = 0.0) -> float:
                for n in names:
                    v = row.get(n)
                    if v is None or v == "" or v == "nan":
                        continue
                    try:
                        return float(v)
                    except ValueError:
                        continue
                return default

            tp = fget("TRUTH.TP")
            fn = fget("TRUTH.FN")
            fp = fget("QUERY.FP")
            prec = fget("METRIC.Precision", default=-1.0)
            rec = fget("METRIC.Recall", default=-1.0)
            f1 = fget("METRIC.F1_Score", "METRIC.F1", default=-1.0)
            if prec < 0 and (tp + fp) > 0:
                prec = tp / (tp + fp)
            if rec < 0 and (tp + fn) > 0:
                rec = tp / (tp + fn)
            if f1 < 0 and prec >= 0 and rec >= 0 and (prec + rec) > 0:
                f1 = 2 * prec * rec / (prec + rec)
            if prec < 0:
                prec = 0.0
            if rec < 0:
                rec = 0.0
            if f1 < 0:
                f1 = 0.0
            out[key] = {
                "precision": prec,
                "recall": rec,
                "f1": f1,
                "truth_tp": tp,
                "truth_fn": fn,
                "query_fp": fp,
                "filter": filt or "ALL",
            }
    return out


def load_region_results(run_dir: pathlib.Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for region_dir in sorted(run_dir.glob("regions/*")):
        if not region_dir.is_dir():
            continue
        meta_path = region_dir / "region.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.is_file() else {}
        rust_sum = region_dir / "happy_rust" / "prefix.summary.csv"
        java_sum = region_dir / "happy_java" / "prefix.summary.csv"
        rust = parse_happy_summary(rust_sum)
        java = parse_happy_summary(java_sum)
        rows.append(
            {
                "region": meta.get("name", region_dir.name),
                "kind": meta.get("kind", ""),
                "span_bp": meta.get("span_bp"),
                "status": meta.get("status", "unknown"),
                "rust": rust,
                "java": java,
                "rust_summary_csv": str(rust_sum) if rust_sum.is_file() else None,
                "java_summary_csv": str(java_sum) if java_sum.is_file() else None,
            }
        )
    return rows


def f1_of(row: dict[str, Any], engine: str, typ: str) -> float | None:
    block = row.get(engine) or {}
    m = block.get(typ)
    if not m:
        return None
    return float(m["f1"])


def find_regressions(
    rows: list[dict[str, Any]],
    baseline: dict[str, Any] | None,
    threshold: float,
) -> list[dict[str, Any]]:
    if not baseline:
        return []
    base_regions = {
        r["region"]: r for r in baseline.get("regions", []) if isinstance(r, dict)
    }
    regs: list[dict[str, Any]] = []
    for row in rows:
        name = row["region"]
        prev = base_regions.get(name)
        if not prev:
            continue
        for typ in ("SNP", "INDEL"):
            cur = f1_of(row, "rust", typ)
            old = f1_of(prev, "rust", typ)
            if cur is None or old is None:
                continue
            drop = old - cur
            if drop > threshold:
                regs.append(
                    {
                        "region": name,
                        "type": typ,
                        "baseline_f1": old,
                        "current_f1": cur,
                        "drop": drop,
                        "threshold": threshold,
                    }
                )
    return regs


def render_markdown(
    summary: dict[str, Any], regressions: list[dict[str, Any]]
) -> str:
    lines = [
        "# Equivalence dashboard (nightly trio E2E)",
        "",
        "Full spine: **HaplotypeCaller (GVCF) → CombineGVCFs → GenotypeGVCFs → VariantFiltration**",
        "for GIAB Ashkenazi trio HG002/HG003/HG004, scored with **hap.py** vs HG002 truth.",
        "",
        f"- **Generated (UTC):** {summary.get('generated_utc')}",
        f"- **Commit:** `{summary.get('commit_sha', '—')}`",
        f"- **Run dir:** `{summary.get('run_dir', '—')}`",
        f"- **Baseline:** `{summary.get('baseline_path', '—')}`",
        f"- **Regression threshold (|ΔF1| drop):** {summary.get('regression_threshold')}",
        f"- **Regressions:** {len(regressions)}",
        "",
        "## Regions",
        "",
        "| Region | Kind | Status | Rust SNP F1 | Rust INDEL F1 | Java SNP F1 | Java INDEL F1 |",
        "|--------|------|--------|------------:|--------------:|------------:|--------------:|",
    ]
    for r in summary.get("regions", []):
        lines.append(
            "| {region} | {kind} | {status} | {rs} | {ri} | {js} | {ji} |".format(
                region=r.get("region"),
                kind=r.get("kind") or "—",
                status=r.get("status") or "—",
                rs=_fmt_f1(f1_of(r, "rust", "SNP")),
                ri=_fmt_f1(f1_of(r, "rust", "INDEL")),
                js=_fmt_f1(f1_of(r, "java", "SNP")),
                ji=_fmt_f1(f1_of(r, "java", "INDEL")),
            )
        )
    lines += ["", "## Regressions vs last green", ""]
    if not regressions:
        lines.append("_None above threshold._")
    else:
        lines += [
            "| Region | Type | Baseline F1 | Current F1 | Drop |",
            "|--------|------|------------:|-----------:|-----:|",
        ]
        for g in regressions:
            lines.append(
                f"| {g['region']} | {g['type']} | {g['baseline_f1']:.4f} | "
                f"{g['current_f1']:.4f} | {g['drop']:.4f} |"
            )
    lines += [
        "",
        "## Notes",
        "",
        "- BAMs are **region-sliced** with `samtools view -L` (no full WGS download).",
        "- Hard regions are capped slices of GIAB stratification BEDs (segdups / TR / alldifficult / MHC).",
        "- Soft gate: regressions open a GitHub issue (`equivalence-regression`); the workflow does not hard-fail.",
        "",
    ]
    return "\n".join(lines) + "\n"


def _fmt_f1(v: float | None) -> str:
    return "—" if v is None else f"{v:.4f}"


def render_html(summary: dict[str, Any], regressions: list[dict[str, Any]]) -> str:
    body = []
    for r in summary.get("regions", []):
        body.append(
            "<tr><td>{region}</td><td>{kind}</td><td>{status}</td>"
            "<td>{rs}</td><td>{ri}</td><td>{js}</td><td>{ji}</td></tr>".format(
                region=escape(str(r.get("region"))),
                kind=escape(str(r.get("kind") or "")),
                status=escape(str(r.get("status") or "")),
                rs=escape(_fmt_f1(f1_of(r, "rust", "SNP"))),
                ri=escape(_fmt_f1(f1_of(r, "rust", "INDEL"))),
                js=escape(_fmt_f1(f1_of(r, "java", "SNP"))),
                ji=escape(_fmt_f1(f1_of(r, "java", "INDEL"))),
            )
        )
    reg_html = "<p><em>None above threshold.</em></p>"
    if regressions:
        rows = "".join(
            "<tr><td>{region}</td><td>{typ}</td><td>{b:.4f}</td>"
            "<td>{c:.4f}</td><td>{d:.4f}</td></tr>".format(
                region=escape(g["region"]),
                typ=escape(g["type"]),
                b=g["baseline_f1"],
                c=g["current_f1"],
                d=g["drop"],
            )
            for g in regressions
        )
        reg_html = (
            "<table><thead><tr><th>Region</th><th>Type</th>"
            "<th>Baseline F1</th><th>Current F1</th><th>Drop</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>"
        )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>gatk-rs nightly trio equivalence</title>
<style>
body {{ font-family: ui-sans-serif, system-ui, sans-serif; margin: 2rem; color: #122; background: #f4f1ea; }}
h1 {{ font-size: 1.45rem; }}
table {{ border-collapse: collapse; background: #fff; margin: 1rem 0; }}
th, td {{ border: 1px solid #ccd; padding: 0.4rem 0.65rem; font-size: 0.9rem; }}
th {{ background: #e8eef5; }}
.meta {{ background: #fff; padding: 0.75rem 1rem; border-left: 4px solid #345; max-width: 56rem; }}
footer {{ margin-top: 2rem; color: #556; font-size: 0.85rem; }}
</style>
</head>
<body>
<h1>gatk-rs · Nightly trio equivalence</h1>
<div class="meta">
<p><strong>Commit:</strong> <code>{escape(str(summary.get('commit_sha', '—')))}</code><br/>
<strong>Generated (UTC):</strong> {escape(str(summary.get('generated_utc', '')))}<br/>
<strong>Regression threshold:</strong> {escape(str(summary.get('regression_threshold')))}</p>
</div>
<h2>Regions</h2>
<table>
<thead><tr>
<th>Region</th><th>Kind</th><th>Status</th>
<th>Rust SNP F1</th><th>Rust INDEL F1</th>
<th>Java SNP F1</th><th>Java INDEL F1</th>
</tr></thead>
<tbody>
{"".join(body)}
</tbody>
</table>
<h2>Regressions vs last green</h2>
{reg_html}
<footer>Independent community project — HC→Combine→Genotype→Filter joint E2E via hap.py.</footer>
</body>
</html>
"""


def issue_body(summary: dict[str, Any], regressions: list[dict[str, Any]]) -> str:
    lines = [
        "## Nightly equivalence regression",
        "",
        f"- **Commit:** `{summary.get('commit_sha')}`",
        f"- **Generated (UTC):** {summary.get('generated_utc')}",
        f"- **Threshold (F1 drop):** {summary.get('regression_threshold')}",
        f"- **Run:** `{summary.get('run_dir')}`",
        "",
        "| Region | Type | Baseline F1 | Current F1 | Drop |",
        "|--------|------|------------:|-----------:|-----:|",
    ]
    for g in regressions:
        lines.append(
            f"| {g['region']} | {g['type']} | {g['baseline_f1']:.4f} | "
            f"{g['current_f1']:.4f} | {g['drop']:.4f} |"
        )
    lines += [
        "",
        "Label: `equivalence-regression`. Soft gate — workflow does not hard-fail.",
        "",
        "Artifacts: see the `nightly-equivalence-*` workflow run.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", type=pathlib.Path, required=True)
    ap.add_argument("--baseline", type=pathlib.Path, required=True)
    ap.add_argument("--dashboard-md", type=pathlib.Path, required=True)
    ap.add_argument("--pages-dir", type=pathlib.Path, required=True)
    ap.add_argument("--commit-sha", default="unknown")
    ap.add_argument("--regression-threshold", type=float, default=0.02)
    args = ap.parse_args()

    rows = load_region_results(args.run_dir)
    baseline = None
    if args.baseline.is_file():
        baseline = json.loads(args.baseline.read_text(encoding="utf-8"))

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    summary: dict[str, Any] = {
        "generated_utc": generated,
        "commit_sha": args.commit_sha,
        "run_dir": str(args.run_dir),
        "baseline_path": str(args.baseline),
        "regression_threshold": args.regression_threshold,
        "regions": rows,
    }

    regressions = find_regressions(rows, baseline, args.regression_threshold)
    summary["regressions"] = regressions
    summary["is_green"] = len(regressions) == 0 and any(
        r.get("status") == "ok" for r in rows
    )

    # JSON artifact (hap.py-derived)
    summary_path = args.run_dir / "happy_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    md = render_markdown(summary, regressions)
    args.dashboard_md.parent.mkdir(parents=True, exist_ok=True)
    args.dashboard_md.write_text(md, encoding="utf-8")

    args.pages_dir.mkdir(parents=True, exist_ok=True)
    (args.pages_dir / "index.html").write_text(
        render_html(summary, regressions), encoding="utf-8"
    )
    (args.pages_dir / "index.md").write_text(md, encoding="utf-8")
    (args.pages_dir / "happy_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    # Flags for the workflow
    flag_dir = args.run_dir / "flags"
    flag_dir.mkdir(parents=True, exist_ok=True)
    (flag_dir / "regressions.json").write_text(
        json.dumps(regressions, indent=2) + "\n", encoding="utf-8"
    )
    (flag_dir / "issue_body.md").write_text(
        issue_body(summary, regressions), encoding="utf-8"
    )
    (flag_dir / "is_green").write_text(
        "1\n" if summary["is_green"] else "0\n", encoding="utf-8"
    )

    if summary["is_green"]:
        # Promote baseline for the workflow to commit
        green = {
            "generated_utc": generated,
            "commit_sha": args.commit_sha,
            "regression_threshold": args.regression_threshold,
            "regions": rows,
        }
        green_path = args.run_dir / "flags" / "new_baseline.json"
        green_path.write_text(json.dumps(green, indent=2) + "\n", encoding="utf-8")

    print(f"[nightly-publish] regions={len(rows)} regressions={len(regressions)}")
    print(f"[nightly-publish] wrote {summary_path}")
    print(f"[nightly-publish] wrote {args.dashboard_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

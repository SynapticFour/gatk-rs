#!/usr/bin/env python3
"""Build Markdown + HTML dashboard from a GIAB genomewide equivalence run directory."""
from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from html import escape


def load_rows(run_dir: pathlib.Path) -> list[dict]:
    path = run_dir / "samples.jsonl"
    if not path.is_file():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def fmt_sec(v) -> str:
    if v is None:
        return "—"
    try:
        return f"{float(v):.1f}s"
    except (TypeError, ValueError):
        return "—"


def fmt_rss(kb) -> str:
    if kb is None:
        return "—"
    try:
        return f"{float(kb) / (1024 * 1024):.2f} GiB"
    except (TypeError, ValueError):
        return "—"


def max_abs_delta(row: dict) -> float | None:
    eq = row.get("equiv_results") or {}
    deltas = eq.get("f1_deltas") or []
    if not deltas:
        return eq.get("max_abs_delta")
    return max(float(d.get("abs_delta", 0.0)) for d in deltas)


def render_markdown(run_dir: pathlib.Path, rows: list[dict], scope: str) -> str:
    lines = [
        "# GIAB genome-wide equivalence dashboard",
        "",
        f"- **Run dir:** `{run_dir}`",
        f"- **Generated (UTC):** {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"- **Scope:** {scope.strip() or '(see SCOPE.txt)'}",
        "",
        "## Samples",
        "",
        "| Sample | Mode | Gate | Max\\|ΔF1\\| | Java wall | Java RSS | Rust wall | Rust RSS |",
        "|--------|------|:----:|----------:|----------:|---------:|----------:|---------:|",
    ]
    for r in rows:
        lines.append(
            "| {sample} | {mode} | {gate} | {delta} | {jw} | {jr} | {rw} | {rr} |".format(
                sample=r.get("sample"),
                mode=r.get("mode"),
                gate="PASS" if r.get("gate_passed") else "FAIL",
                delta=("—" if max_abs_delta(r) is None else f"{max_abs_delta(r):.4f}"),
                jw=fmt_sec((r.get("java_perf") or {}).get("wall_sec")),
                jr=fmt_rss((r.get("java_perf") or {}).get("max_rss_kb")),
                rw=fmt_sec((r.get("rust_perf") or {}).get("wall_sec")),
                rr=fmt_rss((r.get("rust_perf") or {}).get("max_rss_kb")),
            )
        )
    lines += [
        "",
        "## Notes",
        "",
        "- Primary equivalence metric is **Rust−Java F1 delta** via `gatk-rs-equiv` (hap.py / RTG), not absolute F1.",
        "- Wall / RSS from `/usr/bin/time -v` (Linux) or `/usr/bin/time -l` (macOS).",
        "- “Genome-wide” follows `GIAB_MODE` — see `SCOPE.txt` in the run directory.",
        "",
    ]
    return "\n".join(lines) + "\n"


def render_html(run_dir: pathlib.Path, rows: list[dict], scope: str) -> str:
    body_rows = []
    for r in rows:
        gate = "PASS" if r.get("gate_passed") else "FAIL"
        klass = "pass" if r.get("gate_passed") else "fail"
        d = max_abs_delta(r)
        body_rows.append(
            "<tr class='{klass}'><td>{sample}</td><td>{mode}</td><td>{gate}</td>"
            "<td>{delta}</td><td>{jw}</td><td>{jr}</td><td>{rw}</td><td>{rr}</td></tr>".format(
                klass=klass,
                sample=escape(str(r.get("sample"))),
                mode=escape(str(r.get("mode"))),
                gate=gate,
                delta="—" if d is None else f"{d:.4f}",
                jw=escape(fmt_sec((r.get("java_perf") or {}).get("wall_sec"))),
                jr=escape(fmt_rss((r.get("java_perf") or {}).get("max_rss_kb"))),
                rw=escape(fmt_sec((r.get("rust_perf") or {}).get("wall_sec"))),
                rr=escape(fmt_rss((r.get("rust_perf") or {}).get("max_rss_kb"))),
            )
        )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>gatk-rs GIAB equivalence</title>
<style>
body {{ font-family: ui-sans-serif, system-ui, sans-serif; margin: 2rem; color: #122; background: #f7f6f2; }}
h1 {{ font-size: 1.4rem; }}
.scope {{ max-width: 52rem; padding: 0.75rem 1rem; background: #fff; border-left: 4px solid #345; }}
table {{ border-collapse: collapse; margin-top: 1.5rem; background: #fff; }}
th, td {{ border: 1px solid #ccd; padding: 0.45rem 0.7rem; font-size: 0.92rem; }}
th {{ background: #e8eef5; text-align: left; }}
tr.pass td:nth-child(3) {{ color: #0a6; font-weight: 600; }}
tr.fail td:nth-child(3) {{ color: #b00; font-weight: 600; }}
footer {{ margin-top: 2rem; color: #556; font-size: 0.85rem; }}
</style>
</head>
<body>
<h1>gatk-rs · GIAB equivalence dashboard</h1>
<p class="scope"><strong>Scope:</strong> {escape(scope.strip() or "(see SCOPE.txt)")}</p>
<p>Run: <code>{escape(str(run_dir))}</code><br/>
Generated (UTC): {datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}</p>
<table>
<thead><tr>
<th>Sample</th><th>Mode</th><th>Gate</th><th>Max|ΔF1|</th>
<th>Java wall</th><th>Java RSS</th><th>Rust wall</th><th>Rust RSS</th>
</tr></thead>
<tbody>
{"".join(body_rows)}
</tbody>
</table>
<footer>
Independent community project — not Broad-affiliated. Metric: Rust−Java F1 delta via hap.py/RTG (gatk-rs-equiv).
</footer>
</body>
</html>
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", required=True, type=pathlib.Path)
    ap.add_argument("--out-dir", required=True, type=pathlib.Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = load_rows(args.run_dir)
    scope = ""
    scope_path = args.run_dir / "SCOPE.txt"
    if scope_path.is_file():
        scope = scope_path.read_text(encoding="utf-8")
    md = render_markdown(args.run_dir, rows, scope)
    html = render_html(args.run_dir, rows, scope)
    (args.out_dir / "index.md").write_text(md, encoding="utf-8")
    (args.out_dir / "index.html").write_text(html, encoding="utf-8")
    (args.out_dir / "summary.json").write_text(
        json.dumps({"run_dir": str(args.run_dir), "scope": scope, "samples": rows}, indent=2)
        + "\n",
        encoding="utf-8",
    )
    print(f"[giab-dashboard] wrote {args.out_dir / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path


def load_tsv(path: Path) -> dict[str, float]:
    rows: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        case_id, value = line.split("\t")[:2]
        rows[case_id] = float(value)
    return rows


def load_classes(path: Path) -> dict[str, str]:
    classes: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        classes[cols[0]] = cols[5]
    return classes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--java", required=True)
    parser.add_argument("--rust", required=True)
    parser.add_argument("--json-out", required=True)
    parser.add_argument("--md-out", required=True)
    parser.add_argument("--warn-threshold", type=float, default=1.0)
    parser.add_argument("--fail-threshold", type=float, default=5.0)
    parser.add_argument("--gap-open", type=float)
    parser.add_argument("--gap-extend", type=float)
    parser.add_argument("--ins-emission", type=float)
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    java_path = Path(args.java)
    rust_path = Path(args.rust)
    json_out = Path(args.json_out)
    md_out = Path(args.md_out)

    classes = load_classes(matrix_path)
    java = load_tsv(java_path)
    rust = load_tsv(rust_path)

    rows = []
    missing = []
    for case_id, java_ll in java.items():
        if case_id not in rust:
            missing.append(case_id)
            continue
        rust_ll = rust[case_id]
        delta = rust_ll - java_ll
        rows.append(
            {
                "case_id": case_id,
                "class": classes.get(case_id, "unknown"),
                "java_log10_likelihood": java_ll,
                "rust_log10_likelihood": rust_ll,
                "delta": delta,
                "abs_delta": abs(delta),
            }
        )

    rows.sort(key=lambda r: r["abs_delta"], reverse=True)
    abs_values = [r["abs_delta"] for r in rows]
    max_abs = max(abs_values) if abs_values else 0.0
    mean_abs = (sum(abs_values) / len(abs_values)) if abs_values else 0.0
    warn_count = sum(1 for v in abs_values if v > args.warn_threshold)
    fail_count = sum(1 for v in abs_values if v > args.fail_threshold)
    median_delta = statistics.median([r["delta"] for r in rows]) if rows else 0.0
    centered_abs_values = [abs(r["delta"] - median_delta) for r in rows]
    centered_max_abs = max(centered_abs_values) if centered_abs_values else 0.0
    centered_mean_abs = (
        (sum(centered_abs_values) / len(centered_abs_values)) if centered_abs_values else 0.0
    )

    by_class: dict[str, dict[str, float | int]] = {}
    for row in rows:
        cls = row["class"]
        bucket = by_class.setdefault(cls, {"count": 0, "max_abs_delta": 0.0, "mean_abs_delta": 0.0, "sum_abs_delta": 0.0})
        bucket["count"] += 1
        bucket["sum_abs_delta"] += row["abs_delta"]
        bucket["max_abs_delta"] = max(bucket["max_abs_delta"], row["abs_delta"])
    for bucket in by_class.values():
        count = int(bucket["count"])
        bucket["mean_abs_delta"] = bucket["sum_abs_delta"] / count if count else 0.0
        del bucket["sum_abs_delta"]

    status = "pass"
    if fail_count > 0 or missing:
        status = "fail"
    elif warn_count > 0:
        status = "warn"

    payload = {
        "label": "phase6-live-pairhmm-drift",
        "status": status,
        "rust_params": {
            "gap_open_prob": args.gap_open,
            "gap_extend_prob": args.gap_extend,
            "insertion_emission_prob": args.ins_emission,
        },
        "thresholds": {"warn_abs_delta": args.warn_threshold, "fail_abs_delta": args.fail_threshold},
        "summary": {
            "cases_total": len(rows),
            "missing_cases": missing,
            "max_abs_delta": max_abs,
            "mean_abs_delta": mean_abs,
            "median_delta": median_delta,
            "centered_max_abs_delta": centered_max_abs,
            "centered_mean_abs_delta": centered_mean_abs,
            "warn_count": warn_count,
            "fail_count": fail_count,
        },
        "by_class": by_class,
        "rows": rows,
    }
    json_out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md_lines = [
        "# P6 Live PairHMM Drift",
        "",
        "## Config",
        f"- rust gap_open_prob: `{args.gap_open}`",
        f"- rust gap_extend_prob: `{args.gap_extend}`",
        f"- rust insertion_emission_prob: `{args.ins_emission}`",
        "",
        "## Summary",
        f"- status: `{status}`",
        f"- cases: `{len(rows)}`",
        f"- max abs delta: `{max_abs:.6f}`",
        f"- mean abs delta: `{mean_abs:.6f}`",
        f"- median delta: `{median_delta:.6f}`",
        f"- centered max abs delta: `{centered_max_abs:.6f}`",
        f"- centered mean abs delta: `{centered_mean_abs:.6f}`",
        f"- warn count (>{args.warn_threshold}): `{warn_count}`",
        f"- fail count (>{args.fail_threshold}): `{fail_count}`",
        "",
        "## Class Summary",
    ]
    for cls in sorted(by_class):
        b = by_class[cls]
        md_lines.append(
            f"- `{cls}`: count=`{int(b['count'])}` max_abs_delta=`{b['max_abs_delta']:.6f}` mean_abs_delta=`{b['mean_abs_delta']:.6f}`"
        )
    md_lines += ["", "## Top Drift Cases"]
    for row in rows[:8]:
        md_lines.append(
            f"- `{row['case_id']}` ({row['class']}): java=`{row['java_log10_likelihood']:.6f}` rust=`{row['rust_log10_likelihood']:.6f}` abs_delta=`{row['abs_delta']:.6f}`"
        )
    md_out.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return 1 if status == "fail" else 0


if __name__ == "__main__":
    raise SystemExit(main())

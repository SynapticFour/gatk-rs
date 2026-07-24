#!/usr/bin/env python3
"""Append a fair HC comparison run into docs/parity-site/data/perf_history.json."""
from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

ROOT = pathlib.Path(__file__).resolve().parents[2]


def load(path: pathlib.Path) -> dict[str, Any]:
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    return {
        "meta": {
            "title": "gatk-rs performance dashboard",
            "primary_java_baseline": "FASTEST_AVAILABLE",
            "notes": (
                "Speedups are median(Java FASTEST_AVAILABLE wall) / "
                "median(gatk-rs SIMD wall). Never vs LOGLESS_CACHING alone."
            ),
            "updated_utc": None,
        },
        "runs": [],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-json", required=True, type=pathlib.Path)
    ap.add_argument(
        "--site-dir",
        type=pathlib.Path,
        default=ROOT / "docs" / "parity-site",
    )
    ap.add_argument("--commit-sha", default="local")
    ap.add_argument("--workflow-run-url", default="")
    args = ap.parse_args()

    summary = json.loads(args.summary_json.read_text(encoding="utf-8"))
    site = args.site_dir
    hist_path = site / "data" / "perf_history.json"
    latest_path = site / "data" / "perf_latest.json"
    hist = load(hist_path)

    run = {
        "id": summary.get("generated_utc") or datetime.now(timezone.utc).isoformat(),
        "generated_utc": summary.get("generated_utc"),
        "commit_sha": args.commit_sha,
        "workflow_run_url": args.workflow_run_url,
        "workflow": "benchmark",
        "primary_java_baseline": "FASTEST_AVAILABLE",
        "host": summary.get("host") or {},
        "meta": summary.get("meta") or {},
        "speedups": [
            s
            for s in summary.get("speedups") or []
            if not s.get("secondary")
            and s.get("java_baseline_pair_hmm") == "FASTEST_AVAILABLE"
        ],
        "secondary_speedups": [
            s for s in summary.get("speedups") or [] if s.get("secondary")
        ],
        "cells": [
            {
                "region_size": c["region_size"],
                "config_id": c["config_id"],
                "pair_hmm": c.get("pair_hmm"),
                "wall_median_sec": (c.get("wall_sec") or {}).get("median"),
                "wall_stdev_sec": (c.get("wall_sec") or {}).get("stdev"),
                "user_median_sec": (c.get("user_sec") or {}).get("median"),
                "sys_median_sec": (c.get("sys_sec") or {}).get("median"),
                "peak_rss_kb_median": (c.get("peak_rss_kb") or {}).get("median"),
                "energy_joules_median": (c.get("energy_joules") or {}).get("median"),
                "n_ok": c.get("n_ok"),
                "interval": c.get("interval"),
            }
            for c in summary.get("cells") or []
        ],
    }

    # De-dupe by generated_utc + commit
    runs = [
        r
        for r in hist.get("runs") or []
        if not (
            r.get("generated_utc") == run["generated_utc"]
            and r.get("commit_sha") == run["commit_sha"]
        )
    ]
    runs.append(run)
    hist["runs"] = runs
    hist.setdefault("meta", {})
    hist["meta"]["updated_utc"] = datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    hist["meta"]["primary_java_baseline"] = "FASTEST_AVAILABLE"
    hist["meta"]["java_gatk_version"] = (summary.get("meta") or {}).get(
        "gatk_docker", "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    )

    hist_path.parent.mkdir(parents=True, exist_ok=True)
    hist_path.write_text(json.dumps(hist, indent=2) + "\n", encoding="utf-8")
    latest_path.write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {hist_path} ({len(runs)} runs)")
    print(f"Wrote {latest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

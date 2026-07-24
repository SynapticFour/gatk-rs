#!/usr/bin/env python3
"""Aggregate fair HC comparison repeats → median + stddev report (MD + JSON)."""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# Primary honest baseline for speedup charts (not LOGLESS_CACHING).
PRIMARY_JAVA_BASELINE = "java_fastest_available"
PRIMARY_RUST = "rust_simd"


def median(xs: list[float]) -> float | None:
    if not xs:
        return None
    return float(statistics.median(xs))


def stdev(xs: list[float]) -> float | None:
    if len(xs) < 2:
        return 0.0 if xs else None
    return float(statistics.stdev(xs))


def collect(raw_dir: Path) -> dict[str, Any]:
    trials: list[dict[str, Any]] = []
    for p in sorted(raw_dir.glob("**/metrics.json")):
        trials.append(json.loads(p.read_text(encoding="utf-8")))
    return {"trials": trials}


def summarize(trials: list[dict[str, Any]], host: dict[str, Any], meta: dict[str, Any]) -> dict[str, Any]:
    # group by region_size × config_id
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for t in trials:
        key = (t["region_size"], t["config_id"])
        groups.setdefault(key, []).append(t)

    cells: list[dict[str, Any]] = []
    for (region, cfg), rows in sorted(groups.items()):
        def col(name: str) -> list[float]:
            out = []
            for r in rows:
                v = r.get(name)
                if v is not None:
                    out.append(float(v))
            return out

        wall = col("wall_sec")
        user = col("user_sec")
        sys_t = col("sys_sec")
        rss = col("peak_rss_kb")
        energy = col("energy_joules")
        cell = {
            "region_size": region,
            "config_id": cfg,
            "label": rows[0].get("label", cfg),
            "engine": rows[0].get("engine"),
            "pair_hmm": rows[0].get("pair_hmm"),
            "n": len(rows),
            "n_ok": sum(1 for r in rows if r.get("ok")),
            "wall_sec": {"median": median(wall), "stdev": stdev(wall), "samples": wall},
            "user_sec": {"median": median(user), "stdev": stdev(user), "samples": user},
            "sys_sec": {"median": median(sys_t), "stdev": stdev(sys_t), "samples": sys_t},
            "peak_rss_kb": {"median": median(rss), "stdev": stdev(rss), "samples": rss},
            "energy_joules": {
                "median": median(energy),
                "stdev": stdev(energy),
                "samples": energy,
                "available": bool(energy),
            },
            "interval": rows[0].get("interval"),
        }
        cells.append(cell)

    # speedups vs primary Java baseline (per region)
    speedups: list[dict[str, Any]] = []
    by_region: dict[str, dict[str, dict[str, Any]]] = {}
    for c in cells:
        by_region.setdefault(c["region_size"], {})[c["config_id"]] = c

    for region, cfgs in sorted(by_region.items()):
        base = cfgs.get(PRIMARY_JAVA_BASELINE)
        rust = cfgs.get(PRIMARY_RUST)
        if not base or not rust:
            continue
        bw = (base.get("wall_sec") or {}).get("median")
        rw = (rust.get("wall_sec") or {}).get("median")
        if bw and rw and rw > 0:
            speedups.append(
                {
                    "region_size": region,
                    "metric": "wall_sec",
                    "numerator": PRIMARY_JAVA_BASELINE,
                    "denominator": PRIMARY_RUST,
                    "java_baseline_pair_hmm": "FASTEST_AVAILABLE",
                    "speedup": bw / rw,
                    "java_wall_median_sec": bw,
                    "rust_wall_median_sec": rw,
                    "note": (
                        "Speedup = median(Java FASTEST_AVAILABLE wall) / "
                        "median(gatk-rs SIMD wall). Not vs LOGLESS_CACHING."
                    ),
                }
            )
        # Also record vs LOGLESS_CACHING as secondary (clearly labeled cherry-pick risk)
        soft = cfgs.get("java_logless_caching")
        if soft and rust:
            sw = (soft.get("wall_sec") or {}).get("median")
            if sw and rw and rw > 0:
                speedups.append(
                    {
                        "region_size": region,
                        "metric": "wall_sec",
                        "numerator": "java_logless_caching",
                        "denominator": PRIMARY_RUST,
                        "java_baseline_pair_hmm": "LOGLESS_CACHING",
                        "speedup": sw / rw,
                        "java_wall_median_sec": sw,
                        "rust_wall_median_sec": rw,
                        "note": (
                            "SECONDARY reference only — Java software PairHMM path. "
                            "Do not use as the primary marketing speedup."
                        ),
                        "secondary": True,
                    }
                )

    return {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "primary_java_baseline": {
            "config_id": PRIMARY_JAVA_BASELINE,
            "pair_hmm": "FASTEST_AVAILABLE",
            "why": "Native AVX/GKL path when available; honest cloud comparison target.",
        },
        "primary_rust": {
            "config_id": PRIMARY_RUST,
            "pair_hmm": "AVX/SIMD",
        },
        "host": host,
        "meta": meta,
        "cells": cells,
        "speedups": speedups,
        "trials": trials,
    }


def fmt_ms(sec: float | None) -> str:
    if sec is None:
        return "—"
    if sec < 1:
        return f"{sec * 1000:.1f} ms"
    return f"{sec:.3f} s"


def write_md(summary: dict[str, Any], path: Path) -> None:
    host = summary.get("host") or {}
    meta = summary.get("meta") or {}
    lines = [
        "# Fair HaplotypeCaller comparison (dedicated host)",
        "",
        f"**Generated (UTC):** `{summary.get('generated_utc')}`  ",
        f"**Repeats per cell:** `{meta.get('repeats', '?')}` (report = median ± sample stdev)  ",
        f"**Threads:** `{meta.get('threads', 1)}`  ",
        f"**Commit:** `{meta.get('commit_sha', 'local')}`  ",
        "",
        "## Host",
        "",
        f"| Field | Value |",
        f"|-------|--------|",
        f"| CPU | {host.get('cpu_model', 'see HOST_SPECS')} |",
        f"| Logical CPUs | {host.get('logical_cpus', '?')} |",
        f"| RAM GiB | {host.get('mem_gib', '?')} |",
        f"| Kernel | {host.get('kernel', '?')} |",
        f"| AVX2 / AVX-512 | {host.get('simd', {})} |",
        f"| Governor / SMT | {host.get('governor', '?')} / {host.get('smt_control', '?')} |",
        "",
        "## Configurations",
        "",
        "| config_id | Engine | PairHMM | Role |",
        "|-----------|--------|---------|------|",
        "| `rust_logless_scalar` | gatk-rs | `LOGLESS_HMM` (scalar) | Rust scalar baseline |",
        "| `rust_simd` | gatk-rs | `AVX` / SIMD | Rust SIMD under test |",
        "| `java_fastest_available` | Java GATK 4.4 | **`FASTEST_AVAILABLE`** | **Primary fair baseline** |",
        "| `java_logless_caching` | Java GATK 4.4 | `LOGLESS_CACHING` | Secondary software reference |",
        "",
        "> Primary speedup uses **Java `FASTEST_AVAILABLE`** (native AVX when loaded), "
        "not `LOGLESS_CACHING`. Comparing only against the Java software fallback would "
        "be cherry-picking.",
        "",
        "## Results (median ± stdev)",
        "",
    ]

    # table per region
    regions = sorted({c["region_size"] for c in summary["cells"]})
    for region in regions:
        cells = [c for c in summary["cells"] if c["region_size"] == region]
        if not cells:
            continue
        interval = cells[0].get("interval", "?")
        lines += [
            f"### Region `{region}` (`{interval}`)",
            "",
            "| Config | n | Wall | User | Sys | Peak RSS | Energy (J) |",
            "|--------|---|------|------|-----|----------|------------|",
        ]
        for c in cells:
            def cell_s(block: dict[str, Any], unit: str = "s") -> str:
                med = block.get("median")
                sd = block.get("stdev")
                if med is None:
                    return "—"
                if unit == "s":
                    return f"{fmt_ms(med)} ± {fmt_ms(sd) if sd is not None else '—'}"
                if unit == "kb":
                    return f"{med:.0f} ± {(sd or 0):.0f} KiB"
                if unit == "j":
                    if not block.get("available"):
                        return "n/a"
                    return f"{med:.2f} ± {(sd or 0):.2f}"
                return str(med)

            lines.append(
                "| `{cfg}` | {n}/{nok} | {w} | {u} | {s} | {r} | {e} |".format(
                    cfg=c["config_id"],
                    n=c["n"],
                    nok=c["n_ok"],
                    w=cell_s(c["wall_sec"]),
                    u=cell_s(c["user_sec"]),
                    s=cell_s(c["sys_sec"]),
                    r=cell_s(c["peak_rss_kb"], "kb"),
                    e=cell_s(c["energy_joules"], "j"),
                )
            )
        lines.append("")

    lines += ["## Speedups (primary = vs Java FASTEST_AVAILABLE)", ""]
    prim = [s for s in summary["speedups"] if not s.get("secondary")]
    sec = [s for s in summary["speedups"] if s.get("secondary")]
    if prim:
        lines += [
            "| Region | Speedup (Java FASTEST / Rust SIMD) | Java wall | Rust wall |",
            "|--------|--------------------------------------|-----------|-----------|",
        ]
        for s in prim:
            lines.append(
                f"| {s['region_size']} | **{s['speedup']:.3f}×** | "
                f"{fmt_ms(s['java_wall_median_sec'])} | {fmt_ms(s['rust_wall_median_sec'])} |"
            )
        lines.append("")
    else:
        lines.append("_No primary speedup cells (missing baseline or Rust SIMD)._\n")

    if sec:
        lines += [
            "### Secondary (vs Java LOGLESS_CACHING — not for marketing)",
            "",
            "| Region | Speedup | Note |",
            "|--------|---------|------|",
        ]
        for s in sec:
            lines.append(
                f"| {s['region_size']} | {s['speedup']:.3f}× | software Java path |"
            )
        lines.append("")

    lines += [
        "## Repro",
        "",
        "```bash",
        "./scripts/perf/run_fair_hc_comparison.sh",
        "```",
        "",
        "Workflow: `.github/workflows/benchmark.yml` on label `gatk-rs-benchmark`.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True, type=Path)
    ap.add_argument("--host-json", type=Path, default=None)
    ap.add_argument("--meta-json", type=Path, default=None)
    ap.add_argument("--out-json", required=True, type=Path)
    ap.add_argument("--out-md", required=True, type=Path)
    args = ap.parse_args()

    host = {}
    if args.host_json and args.host_json.is_file():
        host = json.loads(args.host_json.read_text(encoding="utf-8"))
    meta = {}
    if args.meta_json and args.meta_json.is_file():
        meta = json.loads(args.meta_json.read_text(encoding="utf-8"))

    packed = collect(args.raw_dir)
    summary = summarize(packed["trials"], host, meta)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_md(summary, args.out_md)
    print(f"Wrote {args.out_json}")
    print(f"Wrote {args.out_md}")
    # Exit non-zero if any required cell missing ok trials
    bad = [c for c in summary["cells"] if c["n_ok"] < 1]
    if bad:
        print(f"ERROR: {len(bad)} cells with zero successful trials", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

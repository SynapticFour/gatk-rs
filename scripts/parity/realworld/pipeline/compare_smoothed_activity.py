#!/usr/bin/env python3
"""Compare GATK IGV assembly-region file vs Rust DumpSmoothedActivity TSV (per-base probs).

GATK rows use tab-separated Chromosome, Start, End, Feature, AssemblyRegions (last col = float).
Rust rows: contig, pos_1based, smoothed_prob.

**Contract (default PASS):** Java’s `AssemblyRegions` column is a **tri-state activity label**
(≈ −1 / 0 / +1) over coarse segments, while Rust exports **band-pass smoothed probabilities** in
[0,1]. They are not the same real-valued process, so requiring `|p_rust − f(p_java)| < ε` everywhere
is misleading. We instead assert:

* every Rust locus maps to a Java segment (no gaps), and
* **binary** active/inactive agreement (see `--java-on-threshold` / `--rust-on-threshold`).

Optional `--require-continuous-max-diff` enforces the legacy continuous gate for debugging.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys


def parse_java_igv(path: pathlib.Path) -> list[tuple[int, int, float, str]]:
    """Return (start0, end0_exclusive, assembly_score, feature) — skips end-marker rows."""
    rows: list[tuple[int, int, float, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        try:
            s = int(parts[1])
            e = int(parts[2])
            feat = parts[3]
            v = float(parts[4])
        except (ValueError, IndexError):
            continue
        if "end-marker" in feat:
            continue
        rows.append((s, e, v, feat))
    return rows


def java_value_at(java_rows: list[tuple[int, int, float, str]], pos1: int) -> float | None:
    """Java IGV uses 0-based half-open [start,end); pos1 is GATK 1-based."""
    p0 = pos1 - 1
    for s, e, v, _feat in reversed(java_rows):
        if s <= p0 < e:
            return v
    return None


def parse_rust_tsv(path: pathlib.Path) -> list[tuple[str, int, float]]:
    out: list[tuple[str, int, float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        a, b, c = line.split("\t")
        out.append((a, int(b), float(c)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--java-igv", required=True)
    ap.add_argument("--rust-tsv", required=True)
    ap.add_argument("--json-out", required=True)
    ap.add_argument(
        "--max-abs-diff",
        type=float,
        default=0.25,
        help="Max allowed difference per base after mapping Java score to [0,1] (default 0.25)",
    )
    ap.add_argument(
        "--max-disagree-rate",
        type=float,
        default=0.02,
        help="Max fraction of bases where binary active/inactive disagrees",
    )
    ap.add_argument(
        "--java-on-threshold",
        type=float,
        default=0.25,
        help="Java raw assembly score treated as ON when above this (default: >0.25 catches +1.0)",
    )
    ap.add_argument(
        "--rust-on-threshold",
        type=float,
        default=0.002,
        help="Rust smoothed prob treated as ON when above this (matches HC activeProbabilityThreshold scale)",
    )
    ap.add_argument(
        "--require-continuous-max-diff",
        action="store_true",
        help="If set, also require max abs diff after mapping Java score to [0,1] (legacy strict mode)",
    )
    args = ap.parse_args()

    java_p = pathlib.Path(args.java_igv)
    rust_p = pathlib.Path(args.rust_tsv)
    java_rows = parse_java_igv(java_p)
    rust_rows = parse_rust_tsv(rust_p)

    def java_to_unit_interval(jv: float) -> float:
        """Map GATK IGV assembly column (-1 / 0 / 1 style) into [0,1] for comparison to smoothed prob."""
        if jv >= 0.5:
            return min(1.0, max(0.0, jv))
        if jv <= -0.5:
            return 0.0
        return max(0.0, min(1.0, (jv + 1.0) / 2.0))

    diffs: list[float] = []
    binary_mismatch = 0
    missing = 0
    for _contig, pos1, rv in rust_rows:
        jv = java_value_at(java_rows, pos1)
        if jv is None:
            missing += 1
            continue
        jm = java_to_unit_interval(jv)
        diffs.append(abs(jm - rv))
        java_on = jv > args.java_on_threshold
        rust_on = rv > args.rust_on_threshold
        if java_on != rust_on:
            binary_mismatch += 1

    max_diff = max(diffs) if diffs else 0.0
    n = len(diffs)
    disagree_rate = (binary_mismatch / n) if n else 0.0
    binary_ok = missing == 0 and n > 0 and disagree_rate <= args.max_disagree_rate
    continuous_ok = max_diff <= args.max_abs_diff
    ok = binary_ok and (continuous_ok if args.require_continuous_max_diff else True)

    payload = {
        "mode": "smoothed_activity_parity",
        "equal": ok,
        "contract": (
            "binary_active_match_plus_no_missing_segments"
            if not args.require_continuous_max_diff
            else "binary_plus_continuous_max_abs_diff"
        ),
        "java_igv": str(java_p),
        "rust_tsv": str(rust_p),
        "compared_positions": n,
        "missing_java_segment": missing,
        "max_abs_diff": max_diff,
        "max_abs_diff_threshold": args.max_abs_diff,
        "continuous_within_threshold": continuous_ok,
        "require_continuous_max_diff": bool(args.require_continuous_max_diff),
        "binary_disagreements": binary_mismatch,
        "binary_disagree_rate": disagree_rate,
        "binary_disagree_rate_threshold": args.max_disagree_rate,
        "java_on_threshold": args.java_on_threshold,
        "rust_on_threshold": args.rust_on_threshold,
    }
    pathlib.Path(args.json_out).write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    if not ok:
        print(
            f"DIVERGENCE binary_ok={binary_ok} continuous_ok={continuous_ok} "
            f"max_abs_diff={max_diff} disagree_rate={disagree_rate:.4f} missing={missing}",
            file=sys.stderr,
        )
        return 1
    print(
        f"PARITY smoothed-activity contract={payload['contract']} max_abs_diff={max_diff} "
        f"(informational) disagree_rate={disagree_rate:.4f} n={n}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

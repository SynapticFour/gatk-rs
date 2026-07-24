#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path


def load_tsv(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        rows.append((parts[0], int(parts[1])))
    rows.sort(key=lambda x: (-x[1], x[0]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expected", required=True)
    ap.add_argument("--actual", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--json-out", required=True)
    args = ap.parse_args()

    expected = load_tsv(Path(args.expected))
    actual = load_tsv(Path(args.actual))

    exp_set = {(s, c) for s, c in expected}
    act_set = {(s, c) for s, c in actual}
    exact_equal = exp_set == act_set

    exp_count = len(expected)
    act_count = len(actual)
    abs_drift = abs(exp_count - act_count)
    drift_stats = {"median_abs_drift": float(abs_drift), "p95_abs_drift": float(abs_drift)}

    result = {
        "label": args.label,
        "expected_count": exp_count,
        "actual_count": act_count,
        "abs_drift": abs_drift,
        "exact_equal": exact_equal,
        "missing": sorted(list(exp_set - act_set)),
        "extra": sorted(list(act_set - exp_set)),
        "drift_stats": drift_stats,
    }
    out = Path(args.json_out)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if exact_equal else 1


if __name__ == "__main__":
    raise SystemExit(main())

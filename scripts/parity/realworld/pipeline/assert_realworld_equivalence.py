#!/usr/bin/env python3
"""Fail fast if `equivalence_report.json` under OUT_DIR violates strict real-world gates.

Use after `run_paired_realworld_pipeline.sh` + `realworld_equivalence_report.py <OUT_DIR>`.
Does not re-run Docker — reads JSON only (reproducible, fast).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "out_dir",
        type=pathlib.Path,
        help="Pipeline OUT_DIR (contains equivalence_report.json)",
    )
    args = ap.parse_args()
    od = args.out_dir.resolve()
    js = od / "equivalence_report.json"
    if not js.is_file():
        print(f"MISSING {js}", file=sys.stderr)
        return 2
    data = json.loads(js.read_text(encoding="utf-8"))
    steps = data.get("steps") or {}
    errs: list[str] = []

    def req(key: str, sub: str, want: str) -> None:
        st = steps.get(key) or {}
        got = st.get(sub)
        if got != want:
            errs.append(f"{key}.{sub}: expected {want!r}, got {got!r}")

    # 02 — operational only
    if (steps.get("02_validate") or {}).get("verdict") not in ("PASS_operational",):
        errs.append("02_validate: expected PASS_operational")

    req("03_count_reads", "verdict", "PARITY")
    req("04_count_bases", "verdict", "PARITY")

    s05 = steps.get("05_filter_reads") or {}
    if s05.get("parity_json"):
        if s05.get("verdict") != "PARITY":
            errs.append(f"05_filter_reads: expected PARITY when parity_json present, got {s05.get('verdict')!r}")

    s06 = steps.get("06_assembly_regions") or {}
    if s06.get("smoothed_parity_json"):
        if s06.get("smoothed_activity_verdict") != "PARITY":
            errs.append(
                f"06 smoothed: expected PARITY when smoothed JSON present, got {s06.get('smoothed_activity_verdict')!r}"
            )
        if s06.get("overall_verdict") not in ("PARITY_JAVA_IGV_AND_SMOOTHED",):
            errs.append(f"06 overall_verdict: expected PARITY_JAVA_IGV_AND_SMOOTHED, got {s06.get('overall_verdict')!r}")

    s07 = steps.get("07_haplotypecaller") or {}
    vs = s07.get("variant_set_verdict")
    if vs not in ("PARITY",):
        errs.append(f"07 variant_set_verdict: expected PARITY, got {vs!r}")

    if errs:
        print("assert_realworld_equivalence: FAILED", file=sys.stderr)
        for e in errs:
            print(f"  - {e}", file=sys.stderr)
        return 1
    print(f"assert_realworld_equivalence: OK ({od})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

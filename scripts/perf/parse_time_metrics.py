#!/usr/bin/env python3
"""Parse GNU/macOS time logs (+ optional perf energy) into JSON metrics."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


def parse_time_log(text: str) -> dict:
    wall = None
    user = None
    sys_t = None
    rss_kb = None

    m = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", text
    )
    if m:
        parts = m.group(1).split(":")
        if len(parts) == 3:
            wall = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        elif len(parts) == 2:
            wall = int(parts[0]) * 60 + float(parts[1])

    m = re.search(r"User time \(seconds\):\s*([0-9.]+)", text)
    if m:
        user = float(m.group(1))
    m = re.search(r"System time \(seconds\):\s*([0-9.]+)", text)
    if m:
        sys_t = float(m.group(1))
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    if m:
        rss_kb = int(m.group(1))

    # macOS /usr/bin/time -l
    if wall is None:
        m = re.search(r"(\d+\.\d+)\s+real\b", text) or re.search(
            r"^\s*real\s+(\d+\.\d+)", text, re.M
        )
        if m:
            wall = float(m.group(1))
    if user is None:
        m = re.search(r"(\d+\.\d+)\s+user\b", text)
        if m:
            user = float(m.group(1))
    if sys_t is None:
        m = re.search(r"(\d+\.\d+)\s+sys\b", text)
        if m:
            sys_t = float(m.group(1))
    if rss_kb is None:
        m = re.search(r"(\d+)\s+maximum resident set size\b", text)
        if m:
            rss_kb = int(m.group(1)) // 1024

    return {
        "wall_sec": wall,
        "user_sec": user,
        "sys_sec": sys_t,
        "peak_rss_kb": rss_kb,
    }


def parse_perf_energy(text: str) -> dict:
    """Best-effort Joules from `perf stat` RAPL counters."""
    joules = None
    # e.g. "12.34 Joules power/energy-pkg/"
    m = re.search(
        r"([0-9]+(?:\.[0-9]+)?)\s+Joules\s+power/energy-(?:pkg|cores)/",
        text,
    )
    if m:
        joules = float(m.group(1))
    else:
        # alternate formatting
        m = re.search(r"power/energy-pkg/\s*#\s*([0-9.]+)\s*Joules", text)
        if m:
            joules = float(m.group(1))
    watt_hours = (joules / 3600.0) if joules is not None else None
    return {"energy_joules": joules, "energy_watt_hours": watt_hours}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("time_log")
    ap.add_argument("--perf-log", default="")
    args = ap.parse_args()
    text = Path(args.time_log).read_text(encoding="utf-8", errors="replace")
    out = parse_time_log(text)
    out["energy_joules"] = None
    out["energy_watt_hours"] = None
    if args.perf_log and Path(args.perf_log).is_file():
        ptext = Path(args.perf_log).read_text(encoding="utf-8", errors="replace")
        out.update(parse_perf_energy(ptext))
    json.dump(out, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

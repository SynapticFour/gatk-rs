#!/usr/bin/env python3
"""Summarize HC_RSS_TRACE wall deltas + RT-first / cache phenotypes."""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("trace", type=Path)
    args = ap.parse_args()
    text = args.trace.read_text(encoding="utf-8", errors="replace")
    phase_ms: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    for line in text.splitlines():
        if "HC_RSS_TRACE phase=" not in line:
            continue
        m = re.search(r"phase=(\S+)", line)
        if not m:
            continue
        phase = m.group(1)
        counts[phase] += 1
        d = re.search(r"delta_ms=([0-9.]+)", line)
        if d:
            phase_ms[phase] += float(d.group(1))
    print(f"file={args.trace}")
    print(
        "rt_first_hit={hits} miss={miss} cache_hit={ch} cache_store={cs} "
        "skip_empty={sk} bushy_skip_dangling={bd} graph_build_begin={gb} "
        "spine_indel_skip={sis}".format(
            hits=counts.get("rt_first_configured_hit", 0),
            miss=counts.get("rt_first_configured_miss", 0),
            ch=counts.get("rt_extract_cache_hit", 0),
            cs=counts.get("rt_extract_cache_store", 0),
            sk=counts.get("rt_supplement_skip_empty_configured", 0)
            + counts.get("merge_rt_skip_empty_configured", 0),
            bd=counts.get("rt_build_skip_dangling_bushy", 0),
            gb=counts.get("rt_graph_build_begin", 0),
            sis=counts.get("parity_spine_indel_skip_cigar_complete", 0),
        )
    )
    top = sorted(phase_ms.items(), key=lambda x: -x[1])[:20]
    print("top phases by sum delta_ms:")
    for p, ms in top:
        print(f"  {p}: {ms / 1000.0:.2f}s (n={counts[p]})")
    print(f"sum_delta_s={sum(phase_ms.values()) / 1000.0:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Normalize GATK HaplotypeCaller --assembly-region-out IGV lines for stable diffs."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def normalize_text(text: str) -> str:
    out: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        if parts[0] == "Chromosome":
            continue
        chrom, start_s, end_s, feature, val_s = parts[:5]
        try:
            val_out = f"{float(val_s):.6f}"
        except ValueError:
            val_out = val_s
        out.append("\t".join([chrom, start_s, end_s, feature, val_out]))
    return "\n".join(out) + ("\n" if out else "")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("input", type=Path, nargs="?", help="IGV file (default: stdin)")
    p.add_argument("-o", "--output", type=Path, help="write here instead of stdout")
    args = p.parse_args()
    if args.input is None:
        text = sys.stdin.read()
    else:
        text = args.input.read_text(encoding="utf-8", errors="replace")
    normalized = normalize_text(text)
    if args.output:
        args.output.write_text(normalized, encoding="utf-8")
    else:
        sys.stdout.write(normalized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

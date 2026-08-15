#!/usr/bin/env python3
"""Concatenate non-overlapping VCF shards (header from first file, body from all).

Used by GIAB HC sharding so timed-out CI jobs can resume completed window shards
without re-running Java/Rust HaplotypeCaller on finished intervals.
Portable fallback when bcftools is missing or refuses plain `.vcf`
(`concat -a` on bcftools ≥1.19 requires bgzip — see `giab_concat_vcfs`).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def concat_vcfs(inputs: list[Path], output: Path) -> None:
    if not inputs:
        raise SystemExit("concat_vcfs: no input VCFs")
    for p in inputs:
        if not p.is_file():
            raise SystemExit(f"concat_vcfs: missing input {p}")

    header: list[str] = []
    bodies: list[str] = []
    chrom_line: str | None = None

    for idx, path in enumerate(inputs):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        file_header: list[str] = []
        file_body: list[str] = []
        for line in lines:
            if line.startswith("##"):
                file_header.append(line)
            elif line.startswith("#CHROM"):
                if chrom_line is None:
                    chrom_line = line
                elif line != chrom_line:
                    raise SystemExit(
                        f"concat_vcfs: #CHROM mismatch between {inputs[0]} and {path}"
                    )
            elif line.strip():
                file_body.append(line)
        if idx == 0:
            header = file_header
        bodies.extend(file_body)

    if chrom_line is None:
        raise SystemExit("concat_vcfs: no #CHROM line in inputs")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        for line in header:
            fh.write(line + "\n")
        fh.write(chrom_line + "\n")
        for line in bodies:
            fh.write(line + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", required=True, type=Path)
    ap.add_argument("inputs", nargs="+", type=Path)
    args = ap.parse_args()
    concat_vcfs(args.inputs, args.output)
    print(f"concat_vcfs: wrote {args.output} from {len(args.inputs)} shard(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

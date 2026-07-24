#!/usr/bin/env python3
"""Compare Java vs Rust HC VCFs from the real-world paired pipeline (no bcftools required)."""
from __future__ import annotations

import pathlib
import sys


def parse_variants(path: pathlib.Path) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        rows.append((chrom, pos, ref, alt))
    return rows


def main() -> int:
    if len(sys.argv) != 4:
        print(
            "usage: compare_realworld_hc_vcfs.py <java.vcf> <rust.vcf> <out.txt>",
            file=sys.stderr,
        )
        return 2
    java_p = pathlib.Path(sys.argv[1])
    rust_p = pathlib.Path(sys.argv[2])
    out_p = pathlib.Path(sys.argv[3])
    if not java_p.is_file() or not rust_p.is_file():
        print("missing input VCF", file=sys.stderr)
        return 1
    j = parse_variants(java_p)
    r = parse_variants(rust_p)
    js, rs = set(j), set(r)
    shared = js & rs
    j_only = sorted(js - rs)
    r_only = sorted(rs - js)
    lines = [
        "# Java vs Rust HaplotypeCaller (CHROM, POS, REF, ALT)",
        f"java_variants: {len(j)} unique ({java_p})",
        f"rust_variants: {len(r)} unique ({rust_p})",
        f"shared: {len(shared)}",
        f"java_only: {len(j_only)}",
        f"rust_only: {len(r_only)}",
        "",
    ]
    if j_only:
        lines.append("## Java-only")
        for t in j_only:
            lines.append("\t".join(t))
        lines.append("")
    if r_only:
        lines.append("## Rust-only")
        for t in r_only:
            lines.append("\t".join(t))
    out_p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out_p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

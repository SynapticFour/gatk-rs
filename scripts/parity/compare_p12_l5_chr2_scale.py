#!/usr/bin/env python3
"""L5.4 chr2 scale gate: variant set + rust-only on a 2k interval."""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

from compare_p12_l5_gvcf import (
    compare_l4_format_gvcf,
    compare_semantic_blocks,
    gvcf_blocks,
    index_variants,
    is_variant_row,
    normalize_alt,
    parse_format_sample,
    parse_vcf_rows,
    variant_key,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--java", type=pathlib.Path, required=True)
    parser.add_argument("--rust", type=pathlib.Path, required=True)
    parser.add_argument("--json-out", type=pathlib.Path)
    parser.add_argument("--md-out", type=pathlib.Path)
    args = parser.parse_args()

    j_rows = parse_vcf_rows(args.java)
    r_rows = parse_vcf_rows(args.rust)
    j_vars = index_variants(j_rows)
    r_vars = index_variants(r_rows)
    j_keys = set(j_vars)
    r_keys = set(r_vars)
    shared = j_keys & r_keys

    l4_mismatch = 0
    for key in sorted(shared, key=lambda k: (k[0], int(k[1]), k[2], k[3])):
        j_fmt = parse_format_sample(j_vars[key].get("format", ""), j_vars[key].get("sample", ""))
        r_fmt = parse_format_sample(r_vars[key].get("format", ""), r_vars[key].get("sample", ""))
        if compare_l4_format_gvcf(r_fmt, j_fmt):
            l4_mismatch += 1

    j_blocks = gvcf_blocks(j_rows)
    r_blocks = gvcf_blocks(r_rows)
    semantic_ok, semantic_stats = compare_semantic_blocks(j_blocks, r_blocks)

    rust_only = r_keys - j_keys
    java_only = j_keys - r_keys
    variant_gate = not rust_only and len(shared) > 0
    block_gate = semantic_ok
    status = "pass" if variant_gate and block_gate else ("variant_pass" if variant_gate else "fail")

    failures: list[str] = []
    if rust_only:
        failures.append(f"rust_only_variants={len(rust_only)}")
    if l4_mismatch:
        failures.append(f"l4_format_mismatch={l4_mismatch}/{len(shared)}")
    if not semantic_ok:
        failures.append(
            f"semantic_block_mismatch partitions={semantic_stats['per_base_gq_partition_mismatch']}"
        )

    payload = {
        "label": "p12-l5-chr2-scale",
        "status": status,
        "interval": "2:92300000-92302000",
        "variant_gate": variant_gate,
        "block_gate": block_gate,
        "block_gate_semantic": semantic_ok,
        "java_variant_count": len(j_keys),
        "rust_variant_count": len(r_keys),
        "shared_variant_count": len(shared),
        "java_only_variants": len(java_only),
        "rust_only_variants": len(rust_only),
        "l4_format_mismatch_count": l4_mismatch,
        "failures": failures,
        "java_gvcf": str(args.java),
        "rust_gvcf": str(args.rust),
    }
    payload.update(semantic_stats)

    md_lines = [
        "# P12 L5.4 chr2 scale",
        "",
        f"- status: **{status}**",
        f"- interval: `2:92300000-92302000`",
        f"- variants java/rust/shared: `{len(j_keys)}` / `{len(r_keys)}` / `{len(shared)}`",
        f"- rust-only: `{len(rust_only)}`",
        f"- L4 mismatches on shared: `{l4_mismatch}`",
        f"- semantic block gate: `{semantic_ok}`",
    ]
    if failures:
        md_lines.append("- failures:")
        for f in failures:
            md_lines.append(f"  - `{f}`")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if args.md_out:
        args.md_out.parent.mkdir(parents=True, exist_ok=True)
        args.md_out.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(
        f"[p12-l5-scale] status={status} variants java={len(j_keys)} rust={len(r_keys)} "
        f"shared={len(shared)} rust_only={len(rust_only)} semantic_block={semantic_ok}",
        flush=True,
    )
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

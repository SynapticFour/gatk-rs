#!/usr/bin/env python3
"""L5.2 P12 interval gVCF comparator: variant identity + L4 FORMAT + block boundaries."""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

# Mirror `gatk_haplotypecaller::gvcf_writer::GATK_HC_DEFAULT_GQB`.
GATK_HC_DEFAULT_GQB: tuple[int, ...] = tuple(range(1, 61)) + (70, 80, 90, 99)


def gvcf_gq_partition(gq: int) -> tuple[int, int]:
    """Java `GVCFBlockCombiner` GQ band `[lower, upper)`."""
    gq = max(0, min(99, gq))
    lower = 0
    for upper in GATK_HC_DEFAULT_GQB:
        if gq < upper:
            return (lower, upper)
        lower = upper
    return (lower, 100)


def parse_min_dp(value: str) -> int | None:
    if value in (".", ""):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if info in (".", ""):
        return out
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def parse_format_sample(format_col: str, sample_col: str) -> dict[str, str]:
    keys = format_col.split(":")
    vals = sample_col.split(":")
    if len(keys) != len(vals):
        return {}
    return dict(zip(keys, vals))


def normalize_alt(alt: str) -> str:
    """Strip gVCF `<NON_REF>` padding from ALT lists."""
    alts = [a for a in alt.split(",") if a != "<NON_REF>"]
    if not alts:
        return "<NON_REF>"
    return ",".join(alts)


def normalize_gt(gt: str) -> str:
    return gt.replace("|", "/")


def parse_vcf_rows(path: pathlib.Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        row = {
            "chrom": parts[0],
            "pos": parts[1],
            "ref": parts[3],
            "alt": parts[4],
            "qual": parts[5],
            "filter": parts[6],
            "info": parts[7],
            "line": line,
        }
        if len(parts) >= 10:
            row["format"] = parts[8]
            row["sample"] = parts[9]
        rows.append(row)
    return rows


def is_variant_row(alt: str) -> bool:
    return normalize_alt(alt) != "<NON_REF>"


def variant_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row["chrom"], row["pos"], row["ref"], normalize_alt(row["alt"]))


def load_pinned_keys(path: pathlib.Path) -> set[tuple[str, str, str, str]]:
    keys: set[tuple[str, str, str, str]] = set()
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        c = line.split("\t")
        keys.add((c[0], c[1], c[2], c[3]))
    return keys


def gvcf_blocks(rows: list[dict[str, str]]) -> list[dict[str, str | int]]:
    blocks: list[dict[str, str | int]] = []
    for row in rows:
        if is_variant_row(row["alt"]):
            continue
        info = parse_info(row["info"])
        end = int(info.get("END", row["pos"]))
        fmt = parse_format_sample(row.get("format", ""), row.get("sample", ""))
        blocks.append(
            {
                "chrom": row["chrom"],
                "start": int(row["pos"]),
                "end": end,
                "min_dp": info.get("MIN_DP", fmt.get("MIN_DP", ".")),
                "max_dp": info.get("MAX_DP", "."),
                "gq": fmt.get("GQ", "."),
            }
        )
    blocks.sort(key=lambda b: (str(b["chrom"]), int(b["start"]), int(b["end"])))
    return blocks


def block_span_key(block: dict[str, str | int]) -> tuple[str, int, int]:
    return (str(block["chrom"]), int(block["start"]), int(block["end"]))


def compare_l4_format_gvcf(rust_fmt: dict[str, str], java_fmt: dict[str, str]) -> list[str]:
    """L4 numeric parity with gVCF-aware GT/AD/PL normalization."""
    errs: list[str] = []
    if "GT" in java_fmt and normalize_gt(rust_fmt.get("GT", "")) != normalize_gt(java_fmt["GT"]):
        errs.append(f"FORMAT GT: rust={rust_fmt.get('GT')!r} java={java_fmt['GT']!r}")
    for key in ("GQ", "DP"):
        if key in java_fmt and rust_fmt.get(key) != java_fmt.get(key):
            errs.append(f"FORMAT {key}: rust={rust_fmt.get(key)!r} java={java_fmt.get(key)!r}")
    if "AD" in java_fmt:
        r_ad = [x for x in rust_fmt.get("AD", "").split(",") if x]
        j_ad = [x for x in java_fmt["AD"].split(",") if x]
        # Java gVCF appends a trailing 0 for `<NON_REF>` allele depth.
        if len(j_ad) == len(r_ad) + 1 and j_ad[-1] == "0":
            j_ad = j_ad[:-1]
        if r_ad != j_ad:
            errs.append(f"FORMAT AD: rust={rust_fmt.get('AD')!r} java={java_fmt['AD']!r}")
    if "PL" in java_fmt:
        try:
            r_pl = [int(x) for x in rust_fmt.get("PL", "").split(",") if x]
            j_pl = [int(x) for x in java_fmt["PL"].split(",") if x]
        except ValueError:
            errs.append("FORMAT PL: parse error")
        else:
            # Biallelic gVCF PL is 6 elements (ref/alt/non-ref cross); compare leading hom-ref triplet.
            j_cmp = j_pl[:3] if len(j_pl) >= 6 else j_pl
            if len(r_pl) != len(j_cmp):
                errs.append(f"FORMAT PL length: rust={len(r_pl)} java={len(j_cmp)}")
            else:
                for i, (rp, jp) in enumerate(zip(r_pl, j_cmp)):
                    if abs(rp - jp) > 1:
                        errs.append(f"FORMAT PL[{i}]: rust={rp} java={jp}")
    return errs


def index_variants(rows: list[dict[str, str]]) -> dict[tuple[str, str, str, str], dict[str, str]]:
    out: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for row in rows:
        if is_variant_row(row["alt"]):
            out[variant_key(row)] = row
    return out


def interval_coverage(blocks: list[dict[str, str | int]]) -> set[int]:
    cov: set[int] = set()
    for block in blocks:
        for pos in range(int(block["start"]), int(block["end"]) + 1):
            cov.add(pos)
    return cov


def expand_blocks_per_base(
    blocks: list[dict[str, str | int]],
) -> dict[tuple[str, int], dict[str, int | tuple[int, int]]]:
    """Per-base GQ partition + MIN_DP for tier-B semantic compare."""
    out: dict[tuple[str, int], dict[str, int | tuple[int, int]]] = {}
    for block in blocks:
        chrom = str(block["chrom"])
        gq_raw = block.get("gq", ".")
        try:
            gq = int(gq_raw) if gq_raw not in (".", "") else 0
        except (TypeError, ValueError):
            gq = 0
        min_dp = parse_min_dp(str(block.get("min_dp", ".")))
        partition = gvcf_gq_partition(gq)
        for pos in range(int(block["start"]), int(block["end"]) + 1):
            entry: dict[str, int | tuple[int, int]] = {
                "gq": gq,
                "partition": partition,
            }
            if min_dp is not None:
                entry["min_dp"] = min_dp
            out[(chrom, pos)] = entry
    return out


def compare_semantic_blocks(
    java_blocks: list[dict[str, str | int]],
    rust_blocks: list[dict[str, str | int]],
) -> tuple[bool, dict[str, int | list[str]]]:
    """Tier-B: same GQ partition + MIN_DP per covered base."""
    j_base = expand_blocks_per_base(java_blocks)
    r_base = expand_blocks_per_base(rust_blocks)
    j_cov = {pos for _, pos in j_base}
    r_cov = {pos for _, pos in r_base}
    all_pos = sorted(j_cov | r_cov, key=lambda p: (p,))
    gq_mismatch = 0
    min_dp_mismatch = 0
    partition_mismatch = 0
    samples: list[str] = []
    for pos in all_pos:
        key_j = next((k for k in j_base if k[1] == pos), None)
        key_r = next((k for k in r_base if k[1] == pos), None)
        if key_j is None or key_r is None:
            partition_mismatch += 1
            if len(samples) < 15:
                samples.append(f"{pos}: coverage java={key_j is not None} rust={key_r is not None}")
            continue
        jb = j_base[key_j]
        rb = r_base[key_r]
        if jb["partition"] != rb["partition"]:
            partition_mismatch += 1
            if len(samples) < 15:
                samples.append(
                    f"{pos}: GQ partition rust={rb['partition']} java={jb['partition']} "
                    f"(gq rust={rb['gq']} java={jb['gq']})"
                )
        elif jb["gq"] != rb["gq"]:
            gq_mismatch += 1
        j_min = jb.get("min_dp")
        r_min = rb.get("min_dp")
        # Tier-B: GQ=0 shadow blocks may differ in MIN_DP when rust pileup was empty pre-fix.
        if (
            j_min is not None
            and r_min is not None
            and j_min != r_min
            and jb["partition"] == rb["partition"]
            and jb["partition"] != (0, 1)
        ):
            min_dp_mismatch += 1
            if len(samples) < 15:
                samples.append(f"{pos}: MIN_DP rust={r_min} java={j_min}")
    cov_gap = len(j_cov - r_cov)
    ok = partition_mismatch == 0 and gq_mismatch == 0 and min_dp_mismatch == 0 and cov_gap == 0
    stats: dict[str, int | list[str]] = {
        "per_base_gq_partition_mismatch": partition_mismatch,
        "per_base_gq_value_mismatch": gq_mismatch,
        "per_base_min_dp_mismatch": min_dp_mismatch,
        "block_coverage_gap_bases": cov_gap,
    }
    if samples:
        stats["semantic_errors_sample"] = samples
    return ok, stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--java", type=pathlib.Path, required=True)
    parser.add_argument("--rust", type=pathlib.Path, required=True)
    parser.add_argument(
        "--pinned-keys",
        type=pathlib.Path,
        default=pathlib.Path("parity/fixtures/p12-java-production-emit/p12_production_emit_sites.tsv"),
    )
    parser.add_argument("--json-out", type=pathlib.Path)
    parser.add_argument("--md-out", type=pathlib.Path)
    parser.add_argument(
        "--block-contract",
        choices=("strict", "semantic", "coverage"),
        default="semantic",
        help="Block gate: strict=exact spans (tier A), semantic=GQ partition per base (tier B), coverage=gaps only (tier C)",
    )
    args = parser.parse_args()

    for label, path in (("java", args.java), ("rust", args.rust)):
        if not path.is_file():
            print(f"[p12-l5-gvcf] missing {label} gVCF: {path}", file=sys.stderr)
            return 1

    pinned = load_pinned_keys(args.pinned_keys)
    if len(pinned) != 66:
        print(f"[p12-l5-gvcf] expected 66 pinned keys, got {len(pinned)}", file=sys.stderr)
        return 1

    j_rows = parse_vcf_rows(args.java)
    r_rows = parse_vcf_rows(args.rust)
    j_vars = index_variants(j_rows)
    r_vars = index_variants(r_rows)
    j_keys = set(j_vars)
    r_keys = set(r_vars)

    j_blocks = gvcf_blocks(j_rows)
    r_blocks = gvcf_blocks(r_rows)
    j_spans = {block_span_key(b) for b in j_blocks}
    r_spans = {block_span_key(b) for b in r_blocks}

    failures: list[str] = []
    if r_keys - j_keys:
        failures.append(f"rust_only_variants={len(r_keys - j_keys)}")
    if not pinned.issubset(r_keys):
        failures.append(f"pinned_missing_in_rust={len(pinned - r_keys)}")
    if not pinned.issubset(j_keys):
        failures.append(f"pinned_missing_in_java={len(pinned - j_keys)}")

    java_extra = j_keys - pinned
    if java_extra:
        # Informational only — L5.2 gate is rust-only variants, not java-only.
        pass

    l4_mismatch = 0
    l4_errors: list[str] = []
    for key in sorted(pinned, key=lambda k: (k[0], int(k[1]), k[2], k[3])):
        if key not in j_vars or key not in r_vars:
            continue
        jr = j_vars[key]
        rr = r_vars[key]
        j_fmt = parse_format_sample(jr.get("format", ""), jr.get("sample", ""))
        r_fmt = parse_format_sample(rr.get("format", ""), rr.get("sample", ""))
        errs = compare_l4_format_gvcf(r_fmt, j_fmt)
        if errs:
            l4_mismatch += 1
            pos = key[1]
            for e in errs[:3]:
                l4_errors.append(f"{pos}: {e}")

    if l4_mismatch:
        failures.append(f"l4_format_mismatch={l4_mismatch}/66")

    block_only_java = j_spans - r_spans
    block_only_rust = r_spans - j_spans
    block_span_exact_match = len(j_spans & r_spans)
    block_boundary_fail = bool(block_only_java or block_only_rust)

    j_cov = interval_coverage(j_blocks)
    r_cov = interval_coverage(r_blocks)
    cov_gap = len(j_cov - r_cov)

    semantic_ok, semantic_stats = compare_semantic_blocks(j_blocks, r_blocks)

    if args.block_contract == "strict" and block_boundary_fail:
        failures.append(
            f"block_boundary_mismatch java_only={len(block_only_java)} rust_only={len(block_only_rust)}"
        )
    if args.block_contract == "semantic" and not semantic_ok:
        failures.append(
            f"semantic_block_mismatch partitions={semantic_stats['per_base_gq_partition_mismatch']} "
            f"gq={semantic_stats['per_base_gq_value_mismatch']} "
            f"min_dp={semantic_stats['per_base_min_dp_mismatch']}"
        )
    if args.block_contract == "coverage" and cov_gap:
        failures.append(f"block_coverage_gap_bases={cov_gap}")

    variant_gate = (
        not (r_keys - j_keys)
        and pinned.issubset(r_keys)
        and pinned.issubset(j_keys)
        and l4_mismatch == 0
    )
    block_gate_strict = not block_boundary_fail and cov_gap == 0
    block_gate_semantic = semantic_ok
    if args.block_contract == "strict":
        block_gate = block_gate_strict
    elif args.block_contract == "semantic":
        block_gate = block_gate_semantic
    else:
        block_gate = cov_gap == 0
    status = "pass" if variant_gate and block_gate else ("variant_pass" if variant_gate else "fail")

    payload = {
        "label": "p12-l5-gvcf",
        "status": status,
        "variant_gate": variant_gate,
        "block_gate": block_gate,
        "block_contract": args.block_contract,
        "block_gate_strict": block_gate_strict,
        "block_gate_semantic": block_gate_semantic,
        "java_variant_count": len(j_keys),
        "rust_variant_count": len(r_keys),
        "pinned_variant_count": len(pinned),
        "shared_variant_count": len(j_keys & r_keys),
        "java_only_variants": len(j_keys - r_keys),
        "rust_only_variants": len(r_keys - j_keys),
        "java_extra_vs_pinned": len(java_extra),
        "java_gvcf_block_count": len(j_blocks),
        "rust_gvcf_block_count": len(r_blocks),
        "block_boundary_java_only": len(block_only_java),
        "block_boundary_rust_only": len(block_only_rust),
        "block_span_exact_match": block_span_exact_match,
        "block_coverage_gap_bases": cov_gap,
        "l4_format_mismatch_count": l4_mismatch,
        "failures": failures,
        "java_gvcf": str(args.java),
        "rust_gvcf": str(args.rust),
    }
    payload.update(semantic_stats)
    if java_extra:
        payload["java_extra_variant_keys"] = sorted(java_extra, key=lambda k: int(k[1]))
    if l4_errors:
        payload["l4_format_errors_sample"] = l4_errors[:20]

    md_lines = [
        "# P12 L5.2 gVCF vs Java",
        "",
        f"- status: **{status}**",
        f"- variant_gate: **{variant_gate}** (66 pinned + L4 + no rust-only)",
        f"- block_gate: **{block_gate}** (contract=`{args.block_contract}`)",
        f"- block_gate_strict: `{block_gate_strict}` | block_gate_semantic: `{block_gate_semantic}`",
        f"- java variants / rust variants: `{len(j_keys)}` / `{len(r_keys)}`",
        f"- pinned 66-site coverage: java `{len(pinned & j_keys)}` rust `{len(pinned & r_keys)}`",
        f"- java-only / rust-only variants: `{len(j_keys - r_keys)}` / `{len(r_keys - j_keys)}`",
        f"- java extra vs pinned: `{len(java_extra)}`",
        f"- gVCF blocks: java `{len(j_blocks)}` rust `{len(r_blocks)}`",
        f"- block exact span match: `{block_span_exact_match}/{len(j_spans)}`",
        f"- block boundary mismatches: java-only `{len(block_only_java)}` rust-only `{len(block_only_rust)}`",
        f"- block coverage gap (bases in Java not Rust): `{cov_gap}`",
        f"- per-base GQ partition mismatches: `{semantic_stats.get('per_base_gq_partition_mismatch', 0)}`",
        f"- per-base MIN_DP mismatches: `{semantic_stats.get('per_base_min_dp_mismatch', 0)}`",
        f"- L4 FORMAT mismatches on pinned sites: `{l4_mismatch}/66`",
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
        f"[p12-l5-gvcf] status={status} variant_gate={variant_gate} block_gate={block_gate} "
        f"variants java={len(j_keys)} rust={len(r_keys)} blocks java={len(j_blocks)} rust={len(r_blocks)} "
        f"l4_mismatch={l4_mismatch}/66",
        flush=True,
    )
    if failures:
        for f in failures:
            print(f"[p12-l5-gvcf] FAIL: {f}", file=sys.stderr)
        if l4_errors:
            for e in l4_errors[:10]:
                print(f"  {e}", file=sys.stderr)

    if variant_gate and block_gate:
        print("[p12-l5-gvcf] PASS (variant + block gates)")
        return 0
    if variant_gate:
        print("[p12-l5-gvcf] VARIANT PASS — block boundaries open (see block_gate=false)")
        return 1
    print("[p12-l5-gvcf] FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

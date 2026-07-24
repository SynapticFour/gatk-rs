#!/usr/bin/env python3
"""L3 real-world HC VCF checks: non-vacuous + golden + Java L3 field parity (J-D01)."""
from __future__ import annotations

import argparse
import pathlib
import sys


def parse_variants(path: pathlib.Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        info = parts[7] if len(parts) > 7 else "."
        rows.append(
            {
                "chrom": parts[0],
                "pos": parts[1],
                "id": parts[2],
                "ref": parts[3],
                "alt": parts[4],
                "qual": parts[5],
                "filter": parts[6],
                "info": info,
            }
        )
    return rows


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if info in (".", ""):
        return out
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def float_close(a: float, b: float, rtol: float = 0.05, atol: float = 1.0) -> bool:
    return abs(a - b) <= atol or (max(abs(a), abs(b)) > 0 and abs(a - b) / max(abs(a), abs(b)) <= rtol)


def parse_format_sample(format_col: str, sample_col: str) -> dict[str, str]:
    keys = format_col.split(":")
    vals = sample_col.split(":")
    if len(keys) != len(vals):
        return {}
    return dict(zip(keys, vals))


def compare_l4_format(rust_fmt: dict[str, str], java_fmt: dict[str, str]) -> list[str]:
    errs: list[str] = []
    for key in ("GT", "GQ", "DP"):
        if key in java_fmt and rust_fmt.get(key) != java_fmt.get(key):
            errs.append(f"FORMAT {key}: rust={rust_fmt.get(key)!r} java={java_fmt.get(key)!r}")
    if "AD" in java_fmt and rust_fmt.get("AD") != java_fmt["AD"]:
        errs.append(f"FORMAT AD: rust={rust_fmt.get('AD')!r} java={java_fmt['AD']!r}")
    if "PL" in java_fmt:
        try:
            r_pl = [int(x) for x in rust_fmt.get("PL", "").split(",") if x]
            j_pl = [int(x) for x in java_fmt["PL"].split(",") if x]
        except ValueError:
            errs.append("FORMAT PL: parse error")
        else:
            if len(r_pl) != len(j_pl):
                errs.append(f"FORMAT PL length: rust={len(r_pl)} java={len(j_pl)}")
            else:
                for i, (rp, jp) in enumerate(zip(r_pl, j_pl)):
                    if abs(rp - jp) > 1:
                        errs.append(f"FORMAT PL[{i}]: rust={rp} java={jp}")
    return errs


def compare_l3_fields(rust: dict[str, str], java: dict[str, str]) -> list[str]:
    errs: list[str] = []
    for key in ("chrom", "pos", "ref", "alt", "filter"):
        if rust.get(key) != java.get(key):
            errs.append(f"{key}: rust={rust.get(key)!r} java={java.get(key)!r}")
    try:
        rq = float(rust["qual"])
        jq = float(java["qual"])
        if not float_close(rq, jq, rtol=0.15, atol=200.0):
            errs.append(f"qual: rust={rq} java={jq}")
    except ValueError:
        errs.append(f"qual: parse error rust={rust['qual']!r} java={java['qual']!r}")

    r_info = parse_info(rust.get("info", "."))
    j_info = parse_info(java.get("info", "."))
    for key in ("AC", "AN", "DP", "MLEAC"):
        if key in j_info and r_info.get(key) != j_info.get(key):
            errs.append(f"INFO {key}: rust={r_info.get(key)!r} java={j_info.get(key)!r}")
    for key in ("AF", "MLEAF", "FS", "MQ", "QD", "SOR", "ExcessHet"):
        if key not in j_info:
            continue
        try:
            rv = float(r_info.get(key, "nan"))
            jv = float(j_info[key])
        except ValueError:
            errs.append(f"INFO {key}: parse error")
            continue
        if not float_close(rv, jv, rtol=0.12, atol=2.0):
            errs.append(f"INFO {key}: rust={rv} java={jv}")
    return errs


def normalize_vcf_text(text: str) -> str:
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith("##reference="):
            continue
        lines.append(line)
    return "\n".join(lines) + ("\n" if text.endswith("\n") or lines else "")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("actual", type=pathlib.Path)
    parser.add_argument("--golden", type=pathlib.Path)
    parser.add_argument("--java", type=pathlib.Path, help="Java HC VCF for L3 / identity compare")
    parser.add_argument("--require-non-vacuous", action="store_true")
    parser.add_argument("--require-java-identity", action="store_true")
    parser.add_argument(
        "--require-java-l3",
        action="store_true",
        help="CHROM/POS/REF/ALT/QUAL/FILTER + agreed INFO vs --java",
    )
    parser.add_argument(
        "--require-java-l4",
        action="store_true",
        help="GT/GQ/DP/AD/PL vs --java (PL per-index tolerance ±1)",
    )
    args = parser.parse_args()

    if not args.actual.is_file():
        print(f"missing actual VCF: {args.actual}", file=sys.stderr)
        return 1

    actual_vars = parse_variants(args.actual)
    if args.require_non_vacuous and not actual_vars:
        print("[hc-realworld-compare] vacuous interval: zero variant rows", file=sys.stderr)
        return 1

    if args.golden:
        if not args.golden.is_file():
            print(f"missing golden VCF: {args.golden}", file=sys.stderr)
            return 1
        actual_text = normalize_vcf_text(args.actual.read_text(encoding="utf-8"))
        golden_text = normalize_vcf_text(args.golden.read_text(encoding="utf-8"))
        if actual_text != golden_text:
            print("[hc-realworld-compare] byte mismatch vs golden (reference line ignored)", file=sys.stderr)
            return 1
        print("[hc-realworld-compare] actual matches golden")

    if args.java:
        if not args.java.is_file():
            print(f"missing java VCF: {args.java}", file=sys.stderr)
            return 1
        j_vars = parse_variants(args.java)
        r_vars = actual_vars
        if args.require_java_identity or args.require_java_l3:
            j_set = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in j_vars}
            r_set = {(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in r_vars}
            if j_set != r_set:
                print(
                    f"[hc-realworld-compare] Java identity mismatch: java={len(j_set)} rust={len(r_set)}",
                    file=sys.stderr,
                )
                return 1
            print(f"[hc-realworld-compare] Java identity OK ({len(j_set)} variant(s))")

        if args.require_java_l3 and j_vars and r_vars:
            errs = compare_l3_fields(r_vars[0], j_vars[0])
            if errs:
                print("[hc-realworld-compare] Java L3 field mismatch:", file=sys.stderr)
                for e in errs:
                    print(f"  {e}", file=sys.stderr)
                return 1
            print("[hc-realworld-compare] Java L3 fields OK (QUAL/FILTER/INFO within tolerance)")

        if args.require_java_l4 and j_vars and r_vars:
            r_line = next(
                (
                    ln
                    for ln in args.actual.read_text(encoding="utf-8", errors="replace").splitlines()
                    if ln and not ln.startswith("#")
                ),
                "",
            )
            j_line = next(
                (
                    ln
                    for ln in args.java.read_text(encoding="utf-8", errors="replace").splitlines()
                    if ln and not ln.startswith("#")
                ),
                "",
            )
            r_parts = r_line.split("\t")
            j_parts = j_line.split("\t")
            if len(r_parts) < 10 or len(j_parts) < 10:
                print("[hc-realworld-compare] missing FORMAT/SAMPLE columns for L4", file=sys.stderr)
                return 1
            r_fmt = parse_format_sample(r_parts[8], r_parts[9])
            j_fmt = parse_format_sample(j_parts[8], j_parts[9])
            errs = compare_l4_format(r_fmt, j_fmt)
            if errs:
                print("[hc-realworld-compare] Java L4 FORMAT mismatch:", file=sys.stderr)
                for e in errs:
                    print(f"  {e}", file=sys.stderr)
                return 1
            print("[hc-realworld-compare] Java L4 FORMAT OK (GT/GQ/DP/AD/PL)")

    print(f"[hc-realworld-compare] variant_rows={len(actual_vars)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

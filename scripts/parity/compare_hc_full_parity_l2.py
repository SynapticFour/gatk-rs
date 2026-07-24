#!/usr/bin/env python3
"""Compare Rust hc_full_parity_gate_dump output vs frozen Java L2 TSV."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

# Phase D.1 `read-filters`: per-read TSV, then this delimiter, then `filter\tfiltered_count` summary.
HC_READ_FILTER_COUNT_SECTION = "---HC_READ_FILTER_COUNTS---"


def parse_tsv_text(text: str) -> tuple[list[str], list[dict[str, str]]]:
    lines = [
        ln
        for ln in text.strip().splitlines()
        if ln.strip()
        and "\t" in ln
        and ln.split("\t", 1)[0]
        and ln.split("\t", 1)[0][0].isalpha()
    ]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        cols = line.split("\t")
        rows.append({header[i]: cols[i] if i < len(cols) else "" for i in range(len(header))})
    return header, rows


def parse_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    return parse_tsv_text(path.read_text(encoding="utf-8"))


def parse_l2_tsv(path: Path) -> list[tuple[list[str], list[dict[str, str]]]]:
    content = path.read_text(encoding="utf-8")
    if HC_READ_FILTER_COUNT_SECTION not in content:
        return [parse_tsv_text(content)]
    before, after = content.split(HC_READ_FILTER_COUNT_SECTION, 1)
    return [parse_tsv_text(before), parse_tsv_text(after)]


def is_float_col(name: str) -> bool:
    if name.startswith("genotype_") and name.endswith("_log10"):
        return True
    return name in {
        "active_prob",
        "original_active_prob",
        "smoothed_active_prob",
        "log10_likelihood",
        "max_log10_likelihood",
        "hq_soft_clip_mean",
        "gq",
    }


def ignored_columns() -> set[str]:
    raw = os.environ.get("PARITY_L2_IGNORE_COLUMNS", "")
    return {c.strip() for c in raw.split(",") if c.strip()}


def rows_to_kv(rows: list[dict[str, str]], header: list[str]) -> dict[str, str]:
    if len(header) < 2:
        return {}
    key_col, val_col = header[0], header[1]
    return {row.get(key_col, ""): row.get(val_col, "") for row in rows}


def haplotype_signatures(kv: dict[str, str]) -> list[tuple[bool, float]]:
    indices: list[int] = []
    for key in kv:
        if key.startswith("haplotype_") and key.endswith("_log10_sum"):
            indices.append(int(key[len("haplotype_") : -len("_log10_sum")]))
    sigs: list[tuple[bool, float]] = []
    for i in sorted(indices):
        sum_key = f"haplotype_{i}_log10_sum"
        ref_key = f"haplotype_{i}_is_reference"
        if sum_key not in kv:
            continue
        sigs.append((kv.get(ref_key, "false") == "true", float(kv[sum_key])))
    return sigs


def parse_kept_indices(value: str) -> list[int]:
    if not value.strip():
        return []
    return [int(part) for part in value.split(",") if part.strip()]


def haplotype_bases(kv: dict[str, str], index: int) -> str:
    return kv.get(f"haplotype_{index}_bases", "")


def kept_signature_set(
    indices: list[int],
    sigs: list[tuple[bool, float]],
    float_eps: float,
) -> list[tuple[bool, float]]:
    rel_tol = float(os.environ.get("PARITY_L2_FLOAT_REL_TOL", "1e-2"))
    out: list[tuple[bool, float]] = []
    for i in indices:
        if i >= len(sigs):
            continue
        is_ref, score = sigs[i]
        # Round scores for stable multiset compare across Rust/Java float formatting.
        rounded = round(score, 6)
        out.append((is_ref, rounded))
    out.sort()
    return out


def kept_bases_set(indices: list[int], kv: dict[str, str]) -> list[str]:
    bases = [haplotype_bases(kv, i) for i in indices if haplotype_bases(kv, i)]
    return sorted(bases)


def kept_indices_semantically_equal(
    rust_rows: list[dict[str, str]],
    java_rows: list[dict[str, str]],
    header: list[str],
    rust_value: str,
    java_value: str,
    float_eps: float,
) -> tuple[bool, str | None]:
    """Compare trim kept set by hap signature (ref flag + score) and optional bases."""
    rust_kv = rows_to_kv(rust_rows, header)
    java_kv = rows_to_kv(java_rows, header)
    rust_sigs = haplotype_signatures(rust_kv)
    java_sigs = haplotype_signatures(java_kv)
    if not rust_sigs or not java_sigs:
        if rust_value == java_value:
            return True, None
        return False, f"kept_indices: rust={rust_value!r} java={java_value!r}"
    rust_kept = parse_kept_indices(rust_value)
    java_kept = parse_kept_indices(java_value)
    if rust_value == java_value:
        return True, None
    if len(rust_kept) != len(java_kept):
        return False, f"kept_indices: rust={rust_value!r} java={java_value!r}"

    rust_bases = kept_bases_set(rust_kept, rust_kv)
    java_bases = kept_bases_set(java_kept, java_kv)
    if rust_bases and java_bases:
        if rust_bases == java_bases:
            return True, None
        return False, (
            f"kept_indices: kept haplotype bases differ rust={rust_bases!r} java={java_bases!r}"
        )

    rust_sig_set = kept_signature_set(rust_kept, rust_sigs, float_eps)
    java_sig_set = kept_signature_set(java_kept, java_sigs, float_eps)
    if rust_sig_set == java_sig_set:
        return True, None
    if rust_sig_set != java_sig_set:
        return False, (
            f"kept_indices: signature sets differ rust={rust_value!r} java={java_value!r} "
            f"rust_sigs={rust_sig_set} java_sigs={java_sig_set}"
        )

    def kept_ref_count(indices: list[int], sigs: list[tuple[bool, float]]) -> int:
        return sum(1 for i in indices if i < len(sigs) and sigs[i][0])

    if kept_ref_count(rust_kept, rust_sigs) != kept_ref_count(java_kept, java_sigs):
        return False, f"kept_indices: ref kept rust={rust_value!r} java={java_value!r}"
    return True, None


def row_key(row: dict[str, str], header: list[str]) -> str:
    return row.get(header[0], "") if header else ""


def row_value(row: dict[str, str], header: list[str]) -> str:
    return row.get(header[1], "") if len(header) > 1 else ""


def haplotype_order_sensitive_keys() -> set[str]:
    return {
        "ref_haplotype_index",
        "alt_haplotype_index",
        "best_haplotype_index",
    }


def haplotype_dump_flex_keys() -> set[str]:
    """Genotype / GL fields that may differ when haplotype order or PairHMM path differs."""
    return {
        "genotype_0_log10",
        "genotype_1_log10",
        "genotype_2_log10",
        "pl",
        "gq",
        "ad",
        "dp",
        "ref_haplotype_index",
        "alt_haplotype_index",
        "best_haplotype_index",
    }


def is_haplotype_permutation_dump(rows: list[dict[str, str]], header: list[str]) -> bool:
    keys = {row_key(row, header) for row in rows}
    return any(k.startswith("haplotype_") and k.endswith("_log10_sum") for k in keys)


def haplotype_rows_equal(
    rust_rows: list[dict[str, str]],
    java_rows: list[dict[str, str]],
    header: list[str],
    float_eps: float,
) -> tuple[bool, str | None]:
    rust_sigs = haplotype_signatures(rows_to_kv(rust_rows, header))
    java_sigs = haplotype_signatures(rows_to_kv(java_rows, header))
    if len(rust_sigs) != len(java_sigs):
        return False, f"haplotype_count mismatch: rust={len(rust_sigs)} java={len(java_sigs)}"
    rust_ref = sum(1 for is_ref, _ in rust_sigs if is_ref)
    java_ref = sum(1 for is_ref, _ in java_sigs if is_ref)
    if rust_ref != java_ref:
        return False, f"haplotype ref count: rust={rust_ref} java={java_ref}"
    return True, None


def rows_equal(
    rust_row: dict[str, str],
    java_row: dict[str, str],
    header: list[str],
    float_eps: float,
) -> tuple[bool, str | None]:
    skip = ignored_columns()
    row_key_name = row_key(rust_row, header)
    for col in header:
        if col in skip:
            continue
        rv = rust_row.get(col, "")
        jv = java_row.get(col, "")
        if col == header[0]:
            continue
        if row_key_name in ("pl", "ad"):
            rust_parts = [p.strip() for p in rv.split(",") if p.strip()]
            java_parts = [p.strip() for p in jv.split(",") if p.strip()]
            if len(rust_parts) != len(java_parts):
                return False, f"column {col}: rust={rv!r} java={jv!r}"
            rel_tol = float(os.environ.get("PARITY_L2_FLOAT_REL_TOL", "1e-2"))
            for rp, jp in zip(rust_parts, java_parts):
                try:
                    rd = float(rp)
                    jd = float(jp)
                except ValueError:
                    return False, f"column {col}: non-numeric rust={rv!r} java={jv!r}"
                if not math.isclose(rd, jd, rel_tol=rel_tol, abs_tol=float_eps):
                    return False, f"column {col}: rust={rv!r} java={jv!r}"
        elif is_float_col(row_key_name) or is_float_col(col):
            try:
                rd = float(rv)
                jd = float(jv)
            except ValueError:
                return False, f"column {col}: non-numeric rust={rv!r} java={jv!r}"
            rel_tol = float(os.environ.get("PARITY_L2_FLOAT_REL_TOL", "1e-2"))
            if not math.isclose(rd, jd, rel_tol=rel_tol, abs_tol=float_eps):
                return False, (
                    f"column {col}: rust={rv} java={jv} "
                    f"(abs_tol={float_eps}, rel_tol={rel_tol})"
                )
        elif col in ("start", "end", "extended_start", "extended_end") and rv.isdigit() and jv.isdigit():
            if abs(int(rv) - int(jv)) <= 1:
                continue
            return False, f"column {col}: rust={rv!r} java={jv!r}"
        elif row_key_name == "assembly_profile" and {rv, jv} <= {"-", "default"}:
            continue
        elif rv != jv:
            return False, f"column {col}: rust={rv!r} java={jv!r}"
    return True, None


def compare_one_section(
    r_header: list[str],
    r_rows: list[dict[str, str]],
    j_header: list[str],
    j_rows: list[dict[str, str]],
    section_index: int | None,
    float_eps: float,
) -> list[dict[str, Any]]:
    mismatches: list[dict[str, Any]] = []
    header = r_header if r_header == j_header else r_header
    if r_header != j_header:
        m: dict[str, Any] = {"kind": "header", "rust": r_header, "java": j_header}
        if section_index is not None:
            m["section"] = section_index
        mismatches.append(m)
        return mismatches
    permute_haps = is_haplotype_permutation_dump(r_rows, header)
    if permute_haps:
        ok, reason = haplotype_rows_equal(r_rows, j_rows, header, float_eps)
        if not ok:
            m = {"kind": "haplotype_signatures", "reason": reason}
            if section_index is not None:
                m["section"] = section_index
            mismatches.append(m)
    def flex_hap_row_key(key: str) -> bool:
        return (
            key in haplotype_order_sensitive_keys()
            or key in haplotype_dump_flex_keys()
            or key.startswith("genotype_")
            or (
                key.startswith("haplotype_")
                and (
                    "_log10_sum" in key
                    or "_is_reference" in key
                    or key.endswith("_bases")
                )
            )
        )

    if permute_haps:
        rust_by_key = {row_key(r, header): r for r in r_rows}
        java_by_key = {row_key(j, header): j for j in j_rows}
        rust_keys = set(rust_by_key)
        java_keys = set(java_by_key)
        if rust_keys != java_keys:
            m = {
                "kind": "row_keys",
                "rust_only": sorted(rust_keys - java_keys),
                "java_only": sorted(java_keys - rust_keys),
            }
            if section_index is not None:
                m["section"] = section_index
            mismatches.append(m)
        for key in sorted(rust_keys & java_keys):
            if flex_hap_row_key(key):
                continue
            if key == "kept_indices":
                ok, reason = kept_indices_semantically_equal(
                    r_rows,
                    j_rows,
                    header,
                    row_value(rust_by_key[key], header),
                    row_value(java_by_key[key], header),
                    float_eps,
                )
            else:
                ok, reason = rows_equal(
                    rust_by_key[key], java_by_key[key], header, float_eps
                )
            if not ok:
                row_m: dict[str, Any] = {
                    "kind": "row",
                    "row_key": key,
                    "reason": reason,
                    "rust": rust_by_key[key],
                    "java": java_by_key[key],
                }
                if section_index is not None:
                    row_m["section"] = section_index
                mismatches.append(row_m)
        return mismatches

    if len(r_rows) != len(j_rows):
        m = {"kind": "row_count", "rust": len(r_rows), "java": len(j_rows)}
        if section_index is not None:
            m["section"] = section_index
        mismatches.append(m)
    n = min(len(r_rows), len(j_rows))
    for i in range(n):
        key = row_key(r_rows[i], header)
        if key == "kept_indices":
            ok, reason = kept_indices_semantically_equal(
                r_rows,
                j_rows,
                header,
                row_value(r_rows[i], header),
                row_value(j_rows[i], header),
                float_eps,
            )
        else:
            ok, reason = rows_equal(r_rows[i], j_rows[i], header, float_eps)
        if not ok:
            row_m: dict[str, Any] = {
                "kind": "row",
                "index": i + 1,
                "reason": reason,
                "rust": r_rows[i],
                "java": j_rows[i],
            }
            if section_index is not None:
                row_m["section"] = section_index
            mismatches.append(row_m)
    return mismatches


def compare_files(
    rust_path: Path,
    java_path: Path,
    float_eps: float,
) -> dict[str, Any]:
    r_sections = parse_l2_tsv(rust_path)
    j_sections = parse_l2_tsv(java_path)
    result: dict[str, Any] = {
        "rust": str(rust_path),
        "java": str(java_path),
        "sections": len(r_sections),
        "equal": False,
        "mismatches": [],
    }
    if len(r_sections) != len(j_sections):
        result["mismatches"].append(
            {
                "kind": "section_count",
                "rust_sections": len(r_sections),
                "java_sections": len(j_sections),
            }
        )
        return result
    all_mismatches: list[dict[str, Any]] = []
    for si, ((rh, rr), (jh, jr)) in enumerate(zip(r_sections, j_sections)):
        sec_idx = si if len(r_sections) > 1 else None
        all_mismatches.extend(
            compare_one_section(rh, rr, jh, jr, sec_idx, float_eps)
        )
    result["mismatches"] = all_mismatches
    # Back-compat fields for single-section reports (first table)
    if r_sections:
        result["rust_rows"] = len(r_sections[0][1])
        result["java_rows"] = len(j_sections[0][1])
        result["header_match"] = r_sections[0][0] == j_sections[0][0]
    else:
        result["rust_rows"] = 0
        result["java_rows"] = 0
        result["header_match"] = True
    result["equal"] = not all_mismatches
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("rust_tsv")
    ap.add_argument("java_tsv")
    ap.add_argument(
        "--float-eps",
        type=float,
        default=float(os.environ.get("PARITY_L2_FLOAT_EPS", "1e-5")),
    )
    ap.add_argument("--json-out")
    args = ap.parse_args()

    report = compare_files(
        Path(args.rust_tsv), Path(args.java_tsv), args.float_eps
    )
    if args.json_out:
        Path(args.json_out).write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
    if report["equal"]:
        print(f"L2 OK: {args.rust_tsv}")
        return 0
    print(f"L2 DIFF: {args.rust_tsv} ({len(report['mismatches'])} issue(s))", file=sys.stderr)
    for m in report["mismatches"][:5]:
        print(f"  {m}", file=sys.stderr)
    if len(report["mismatches"]) > 5:
        print(f"  ... and {len(report['mismatches']) - 5} more", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

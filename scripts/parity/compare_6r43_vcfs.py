#!/usr/bin/env python3
"""Compare Java vs Rust HC VCFs for 6R.43 holdouts.

Classification:
  A = equivalent (alleles + FORMAT + INFO QUAL/MLE/QD within tolerance)
  B = representation-only (same alleles/likelihood meaning; formatting/order/float print)
  C = proven implementation divergence (allele set, GT/PL, or INFO that changes meaning)
  D = Java behavior unavailable/ambiguous

Does not patch production code.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

QUAL_TOL = 0.05
QD_TOL = 0.02
AF_TOL = 0.02

INFO_KEYS = ("AC", "AF", "AN", "MLEAC", "MLEAF", "QD")


def parse_vcf_records(path: Path) -> list[dict]:
    rows = []
    if not path.is_file():
        return rows
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 10:
            continue
        chrom, pos, _id, ref, alt, qual, _filt, info, fmt, sample = f[:10]
        info_d = {}
        for part in info.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                info_d[k] = v
        fmt_keys = fmt.split(":")
        fmt_vals = sample.split(":")
        sample_d = dict(zip(fmt_keys, fmt_vals))
        alts = alt.split(",")
        rows.append(
            {
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alts[0],
                "alts": alts,
                "qual": float(qual) if qual not in (".", "") else None,
                "info": info_d,
                "gt": sample_d.get("GT"),
                "ad": sample_d.get("AD"),
                "dp": sample_d.get("DP"),
                "gq": sample_d.get("GQ"),
                "pl": sample_d.get("PL"),
            }
        )
    return rows


def key(r: dict) -> tuple:
    return (r["chrom"], r["pos"], r["ref"], r["alt"])


def close(a, b, tol) -> bool:
    if a is None or b is None:
        return a is b
    try:
        return abs(float(a) - float(b)) <= tol
    except ValueError:
        return str(a) == str(b)


def classify_pair(j: dict, r: dict) -> tuple[str, list[str], str]:
    """Return (class, mismatch notes, first_divergent_stage)."""
    notes = []
    # Alleles already keyed; FORMAT
    for field in ("gt", "ad", "dp", "gq", "pl"):
        jv, rv = j.get(field), r.get(field)
        if jv != rv:
            notes.append(f"{field.upper()} java={jv} rust={rv}")
    if notes:
        return "C", notes, "H"
    # QUAL / MLE / QD
    if not close(j.get("qual"), r.get("qual"), QUAL_TOL):
        notes.append(f"QUAL java={j.get('qual')} rust={r.get('qual')}")
        return "C", notes, "I"
    for k in INFO_KEYS:
        jv = j["info"].get(k)
        rv = r["info"].get(k)
        if k in ("AF", "MLEAF", "QD") or k == "QUAL":
            if not close(jv, rv, QD_TOL if k == "QD" else AF_TOL):
                notes.append(f"{k} java={jv} rust={rv}")
        else:
            if jv != rv:
                notes.append(f"{k} java={jv} rust={rv}")
    if notes:
        return "C", notes, "I"
    return "A", [], "J"


def compare_region(java_vcf: Path, rust_vcf: Path) -> dict:
    jr = parse_vcf_records(java_vcf)
    rr = parse_vcf_records(rust_vcf)
    jmap = {key(x): x for x in jr}
    rmap = {key(x): x for x in rr}
    only_j = sorted(set(jmap) - set(rmap))
    only_r = sorted(set(rmap) - set(jmap))
    shared = sorted(set(jmap) & set(rmap))
    site_rows = []
    classes = []
    first_div = "J" if jr or rr else "J"
    for k in shared:
        cls, notes, stage = classify_pair(jmap[k], rmap[k])
        classes.append(cls)
        if cls == "C" and first_div == "J":
            first_div = stage
        site_rows.append(
            {
                "site": f"{k[0]}:{k[1]} {k[2]}/{k[3]}",
                "class": cls,
                "notes": notes,
                "java_qd": jmap[k]["info"].get("QD"),
                "rust_qd": rmap[k]["info"].get("QD"),
                "java_qual": jmap[k].get("qual"),
                "rust_qual": rmap[k].get("qual"),
                "java_gt": jmap[k].get("gt"),
                "rust_gt": rmap[k].get("gt"),
                "java_pl": jmap[k].get("pl"),
                "rust_pl": rmap[k].get("pl"),
            }
        )
    if only_j or only_r:
        first_div = "F"
        overall = "C"
    elif "C" in classes:
        overall = "C"
    elif classes:
        overall = "A"
    elif not jr and not rr:
        overall = "A"  # both empty
        first_div = "J"
    else:
        overall = "D"

    return {
        "java_n": len(jr),
        "rust_n": len(rr),
        "shared_n": len(shared),
        "java_only": [f"{c}:{p} {r}/{a}" for c, p, r, a in only_j],
        "rust_only": [f"{c}:{p} {r}/{a}" for c, p, r, a in only_r],
        "overall": overall,
        "first_divergent_stage": first_div if overall == "C" else ("J" if overall == "A" else "?"),
        "sites": site_rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--json-out", required=True)
    args = ap.parse_args()
    panel = json.loads(Path(args.panel).read_text())
    out_root = Path(args.out_root)
    summary = {"regions": [], "aggregate": {}}
    counts = {"A": 0, "B": 0, "C": 0, "D": 0, "vcf_mismatch": 0, "internal_only": 0}
    for spec in panel["regions"]:
        rid = spec["id"]
        j = out_root / rid / "java.vcf"
        r = out_root / rid / "rust.vcf"
        row = {
            "id": rid,
            "interval": spec["interval"],
            "role": spec["role"],
            "why": spec["why"],
            "java_vcf": str(j) if j.is_file() else None,
            "rust_vcf": str(r) if r.is_file() else None,
        }
        if not j.is_file() and not r.is_file():
            row["overall"] = "D"
            row["note"] = "missing both VCFs"
            counts["D"] += 1
        elif not j.is_file():
            row["overall"] = "D"
            row["note"] = "missing java.vcf"
            counts["D"] += 1
        elif not r.is_file():
            row["overall"] = "D"
            row["note"] = "missing rust.vcf"
            counts["D"] += 1
        else:
            cmp_ = compare_region(j, r)
            row.update(cmp_)
            counts[cmp_["overall"]] = counts.get(cmp_["overall"], 0) + 1
            if cmp_["overall"] == "C":
                counts["vcf_mismatch"] += 1
        summary["regions"].append(row)
    n = len(summary["regions"])
    summary["aggregate"] = {
        "regions_tested": n,
        "fully_equivalent": counts["A"],
        "representation_only": counts["B"],
        "divergent": counts["C"],
        "unknown": counts["D"],
        "final_vcf_mismatches": counts["vcf_mismatch"],
        "internal_only_mismatches": counts["internal_only"],
    }
    Path(args.json_out).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["aggregate"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

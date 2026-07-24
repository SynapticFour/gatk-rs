#!/usr/bin/env python3
"""Compare Java vs Rust GenotypeGVCFs callsets (site identity + GT + rough QUAL)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_body(path: Path) -> Dict[Tuple[str, int], dict]:
    rows: Dict[Tuple[str, int], dict] = {}
    samples: List[str] = []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith("#CHROM"):
            samples = ln.split("\t")[9:]
            continue
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        chrom, pos = f[0], int(f[1])
        ref, alt = f[3], f[4]
        qual = None if f[5] == "." else float(f[5])
        fmt = f[8].split(":") if len(f) > 8 else []
        gts: Dict[str, Optional[str]] = {}
        for i, name in enumerate(samples):
            if i >= len(f) - 9:
                gts[name] = None
                continue
            parts = f[9 + i].split(":")
            fmap = {k: parts[j] if j < len(parts) else "." for j, k in enumerate(fmt)}
            gts[name] = fmap.get("GT")
        rows[(chrom, pos)] = {
            "ref": ref,
            "alt": set(a for a in alt.split(",") if a and a != "."),
            "qual": qual,
            "gt": gts,
        }
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--java", required=True)
    ap.add_argument("--rust", required=True)
    ap.add_argument("--label", default="genotype-gvcfs")
    ap.add_argument("--qual-tol", type=float, default=20.0)
    args = ap.parse_args()

    ja = parse_body(Path(args.java))
    rb = parse_body(Path(args.rust))
    keys = sorted(set(ja) | set(rb))
    mismatches = 0
    for key in keys:
        if key not in ja:
            print(f"[{args.label}] only-rust {key}")
            mismatches += 1
            continue
        if key not in rb:
            print(f"[{args.label}] only-java {key}")
            mismatches += 1
            continue
        a, b = ja[key], rb[key]
        if a["ref"] != b["ref"] or a["alt"] != b["alt"]:
            print(f"[{args.label}] alleles {key}: java={a['ref']}/{sorted(a['alt'])} rust={b['ref']}/{sorted(b['alt'])}")
            mismatches += 1
        common = set(a["gt"]) & set(b["gt"])
        for s in sorted(common):
            if a["gt"][s] != b["gt"][s]:
                # Normalize ././. vs .
                ag = a["gt"][s] or "."
                bg = b["gt"][s] or "."
                if ag.replace(".", "") == "" and bg.replace(".", "") == "":
                    continue
                print(f"[{args.label}] GT {key} {s}: java={a['gt'][s]} rust={b['gt'][s]}")
                mismatches += 1
        if a["qual"] is not None and b["qual"] is not None:
            if abs(a["qual"] - b["qual"]) > args.qual_tol:
                print(
                    f"[{args.label}] QUAL {key}: java={a['qual']:.2f} rust={b['qual']:.2f} "
                    f"(tol={args.qual_tol})"
                )
                mismatches += 1

    if mismatches:
        print(f"[{args.label}] FAIL mismatches={mismatches} java={len(ja)} rust={len(rb)}")
        return 1
    print(f"[{args.label}] OK sites={len(keys)} (alleles/GT/QUAL±{args.qual_tol})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

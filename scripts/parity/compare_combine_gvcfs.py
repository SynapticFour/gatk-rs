#!/usr/bin/env python3
"""Compare Java vs Rust CombineGVCFs outputs on allele/PL semantics.

Not a full header-identity check. Compares body rows keyed by (CHROM, POS):
REF, ALT set, and per-sample PL vectors when both sides emit the site.
"""
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
        fmt = f[8].split(":") if len(f) > 8 else []
        sample_cols = f[9:] if len(f) > 9 else []
        pl_by_sample: Dict[str, Optional[str]] = {}
        for i, name in enumerate(samples):
            if i >= len(sample_cols):
                pl_by_sample[name] = None
                continue
            parts = sample_cols[i].split(":")
            fmap = {k: parts[j] if j < len(parts) else "." for j, k in enumerate(fmt)}
            pl_by_sample[name] = fmap.get("PL")
        end = None
        if f[7] != ".":
            for item in f[7].split(";"):
                if item.startswith("END="):
                    end = int(item.split("=", 1)[1])
        rows[(chrom, pos)] = {
            "ref": ref,
            "alt": set(alt.split(",")) if alt != "." else set(),
            "end": end,
            "pl": pl_by_sample,
            "samples": samples,
        }
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--java", required=True)
    ap.add_argument("--rust", required=True)
    ap.add_argument("--label", default="combine-gvcfs")
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
        if a["ref"] != b["ref"]:
            print(f"[{args.label}] REF mismatch {key}: java={a['ref']} rust={b['ref']}")
            mismatches += 1
        if a["alt"] != b["alt"]:
            print(f"[{args.label}] ALT mismatch {key}: java={sorted(a['alt'])} rust={sorted(b['alt'])}")
            mismatches += 1
        # Compare PL for samples present on both sides (name intersection).
        common = set(a["pl"]) & set(b["pl"])
        for s in sorted(common):
            if a["pl"][s] != b["pl"][s]:
                print(
                    f"[{args.label}] PL mismatch {key} sample={s}: "
                    f"java={a['pl'][s]} rust={b['pl'][s]}"
                )
                mismatches += 1

    if mismatches:
        print(f"[{args.label}] FAIL mismatches={mismatches} java_sites={len(ja)} rust_sites={len(rb)}")
        return 1
    print(f"[{args.label}] OK sites={len(keys)} (REF/ALT/PL)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

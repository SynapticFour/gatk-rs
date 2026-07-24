#!/usr/bin/env python3
"""Build nightly trio equivalence region BEDs + GATK interval lists.

Produces:
  - full chr20 / chr21 (hs37d5 numeric contigs)
  - 3–5 hard regions carved from GIAB stratification BEDs (segdups, tandem
    repeats, alldifficult, MHC), each capped by a base-pair budget

Outputs under --out-dir:
  manifest.json, *.bed (0-based for samtools -L), *.intervals (1-based for -L)
"""
from __future__ import annotations

import argparse
import gzip
import json
import pathlib
from typing import Iterable


# hs37d5 lengths (numeric contigs, matching parity/realworld/assets/hs37d5.simple.fa)
CHROM_LEN = {
    "20": 63_025_520,
    "21": 48_129_895,
    "6": 171_115_067,
}


def norm_chrom(raw: str) -> str:
    c = raw.strip()
    if c.startswith("chr"):
        c = c[3:]
    return c


def open_bed(path: pathlib.Path):
    if path.suffix == ".gz" or str(path).endswith(".bed.gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def load_bed_intervals(
    path: pathlib.Path, chroms: set[str]
) -> list[tuple[str, int, int]]:
    out: list[tuple[str, int, int]] = []
    with open_bed(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom = norm_chrom(cols[0])
            if chrom not in chroms:
                continue
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            if end > start:
                out.append((chrom, start, end))
    out.sort()
    return out


def merge(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    if not intervals:
        return []
    merged = [intervals[0]]
    for chrom, start, end in intervals[1:]:
        pc, ps, pe = merged[-1]
        if chrom == pc and start <= pe:
            merged[-1] = (pc, ps, max(pe, end))
        else:
            merged.append((chrom, start, end))
    return merged


def take_budget(
    intervals: Iterable[tuple[str, int, int]], budget_bp: int
) -> list[tuple[str, int, int]]:
    """Prefer largest intervals first until budget is exhausted."""
    ranked = sorted(intervals, key=lambda x: -(x[2] - x[1]))
    chosen: list[tuple[str, int, int]] = []
    used = 0
    for chrom, start, end in ranked:
        length = end - start
        if length <= 0:
            continue
        if used >= budget_bp:
            break
        if used + length > budget_bp:
            end = start + (budget_bp - used)
            length = end - start
        chosen.append((chrom, start, end))
        used += length
    return merge(sorted(chosen))


def write_bed(path: pathlib.Path, intervals: list[tuple[str, int, int]]) -> int:
    bp = 0
    with path.open("w", encoding="utf-8") as fh:
        for chrom, start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\n")
            bp += end - start
    return bp


def write_intervals(path: pathlib.Path, intervals: list[tuple[str, int, int]]) -> None:
    """1-based inclusive GATK-style interval strings (one per line)."""
    with path.open("w", encoding="utf-8") as fh:
        for chrom, start, end in intervals:
            # BED half-open → closed 1-based
            fh.write(f"{chrom}:{start + 1}-{end}\n")


def add_region(
    regions: list[dict],
    out_dir: pathlib.Path,
    name: str,
    kind: str,
    intervals: list[tuple[str, int, int]],
) -> None:
    if not intervals:
        return
    bed = out_dir / f"{name}.bed"
    iv = out_dir / f"{name}.intervals"
    bp = write_bed(bed, intervals)
    write_intervals(iv, intervals)
    regions.append(
        {
            "name": name,
            "kind": kind,
            "bed": str(bed),
            "intervals_file": str(iv),
            "n_intervals": len(intervals),
            "span_bp": bp,
        }
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--strat-root", type=pathlib.Path, required=True)
    ap.add_argument("--out-dir", type=pathlib.Path, required=True)
    ap.add_argument(
        "--hard-budget-bp",
        type=int,
        default=2_000_000,
        help="Max bases per hard-region BED (default 2 Mb)",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    regions: list[dict] = []

    # Full chromosomes (medium regions)
    add_region(
        regions,
        args.out_dir,
        "chr20",
        "chromosome",
        [("20", 0, CHROM_LEN["20"])],
    )
    add_region(
        regions,
        args.out_dir,
        "chr21",
        "chromosome",
        [("21", 0, CHROM_LEN["21"])],
    )

    strat = {
        "segdups": args.strat_root / "GRCh37_segdups.bed.gz",
        "tandem_repeats": args.strat_root
        / "GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
        "alldifficult": args.strat_root / "GRCh37_alldifficultregions.bed.gz",
        "mhc": args.strat_root / "GRCh37_MHC.bed.gz",
    }

    hard_specs = [
        ("hard_segdups_chr20", "segdups", {"20"}),
        ("hard_tr_chr20", "tandem_repeats", {"20"}),
        ("hard_alldifficult_chr21", "alldifficult", {"21"}),
        ("hard_segdups_chr21", "segdups", {"21"}),
        ("hard_mhc", "mhc", {"6"}),
    ]

    for name, key, chroms in hard_specs:
        path = strat[key]
        if not path.is_file():
            print(f"[nightly-regions] skip {name}: missing {path}")
            continue
        ivals = take_budget(load_bed_intervals(path, chroms), args.hard_budget_bp)
        add_region(regions, args.out_dir, name, "hard", ivals)

    manifest = {
        "assembly": "hs37d5",
        "contig_style": "numeric",
        "hard_budget_bp": args.hard_budget_bp,
        "regions": regions,
    }
    (args.out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[nightly-regions] wrote {len(regions)} regions → {args.out_dir}")
    for r in regions:
        print(f"  - {r['name']}: {r['span_bp']:,} bp ({r['n_intervals']} intervals)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

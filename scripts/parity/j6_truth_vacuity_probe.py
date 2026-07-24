#!/usr/bin/env python3
"""Probe whether a BAM has non-vacuous depth on GIAB truth sites in an interval (R3).

Exit 0 when at least ``--min-covered-truth`` GIAB truth variants in the eval
interval ∩ high-confidence BED have ``samtools view -c`` depth ≥ ``--min-depth``.

Usage:
  python3 scripts/parity/j6_truth_vacuity_probe.py \\
    --bam parity/realworld/na12878_giab_window_b37/NA12878_giab_window.b37.bam \\
    --truth-vcf parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \\
    --regions-bed parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.bed \\
    --eval-interval 20:10000000-10050000 \\
    --min-covered-truth 5 --min-depth 5
"""
from __future__ import annotations

import argparse
import gzip
import json
import pathlib
import subprocess
import sys

# Reuse interval/BED helpers from p13.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from p13_truth_eval import (  # noqa: E402
    canon_contig,
    in_eval_interval,
    in_regions,
    load_regions,
    parse_eval_interval,
)


def open_text(path: pathlib.Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def bam_depth(bam: pathlib.Path, chrom: str, pos1: int) -> int:
    region = f"{chrom}:{pos1}-{pos1}"
    out = subprocess.check_output(["samtools", "view", "-c", str(bam), region], text=True)
    return int(out.strip())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True, type=pathlib.Path)
    ap.add_argument("--truth-vcf", required=True, type=pathlib.Path)
    ap.add_argument("--regions-bed", type=pathlib.Path, default=None)
    ap.add_argument("--eval-interval", required=True)
    ap.add_argument("--min-depth", type=int, default=5)
    ap.add_argument("--min-covered-truth", type=int, default=5)
    ap.add_argument("--max-sites-scan", type=int, default=5000)
    ap.add_argument("--json-out", type=pathlib.Path, default=None)
    args = ap.parse_args()

    interval = parse_eval_interval(args.eval_interval)
    if interval is None:
        print(f"[vacuity] bad interval: {args.eval_interval}", file=sys.stderr)
        return 2
    regions = load_regions(args.regions_bed)

    truth_in_scope = 0
    covered = 0
    examples: list[dict] = []
    with open_text(args.truth_vcf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom = canon_contig(cols[0])
            try:
                pos1 = int(cols[1])
            except ValueError:
                continue
            if not in_eval_interval(chrom, pos1, interval):
                continue
            if not in_regions(chrom, pos1, regions):
                continue
            truth_in_scope += 1
            if truth_in_scope > args.max_sites_scan:
                break
            # Prefer numeric contig style matching BAM (b37 has "20" not "chr20").
            d = bam_depth(args.bam, chrom, pos1)
            if d >= args.min_depth:
                covered += 1
                if len(examples) < 8:
                    examples.append(
                        {
                            "chrom": chrom,
                            "pos": pos1,
                            "ref": cols[3],
                            "alt": cols[4],
                            "depth": d,
                        }
                    )

    vacuous = covered < args.min_covered_truth
    payload = {
        "label": "j6-truth-vacuity-probe",
        "bam": str(args.bam),
        "eval_interval": args.eval_interval,
        "truth_in_scope": truth_in_scope,
        "covered_truth_sites": covered,
        "min_depth": args.min_depth,
        "min_covered_truth": args.min_covered_truth,
        "vacuous": vacuous,
        "examples": examples,
    }
    text = json.dumps(payload, indent=2) + "\n"
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")
    if vacuous:
        print(
            f"[vacuity] VACUOUS: covered={covered} < min={args.min_covered_truth}",
            file=sys.stderr,
        )
        return 1
    print(f"[vacuity] NON-VACUOUS: covered={covered}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

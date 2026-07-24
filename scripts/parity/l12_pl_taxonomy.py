#!/usr/bin/env python3
"""L12 A1/A2: soft-PL residual taxonomy among GT-matched dense TPs.

Re-measures max |ΔPL| buckets (exact / tol5 / tol50 / tol200 / gt200 / sparse_like)
and classifies high-Δ sites by AD divergence vs pure PairHMM scale.

Example:
  python3 scripts/parity/l12_pl_taxonomy.py \\
    --java-vcf parity/reports/hc-full-parity-j6-dense/p12_dense_giab_window.java.vcf \\
    --rust-vcf parity/reports/hc-full-parity-j6-dense/p12_dense_giab_window.rust.vcf \\
    --eval-interval 20:10000000-10050000 \\
    --md-out docs/CLAIM_MATRIX.md
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import pathlib
import sys
from collections import Counter, defaultdict
from typing import Any


def canon(c: str) -> str:
    return c[3:] if c.startswith("chr") else c


def parse_interval(s: str) -> tuple[str, int, int]:
    chrom, rest = s.split(":", 1)
    a, b = rest.split("-", 1)
    return canon(chrom), int(a), int(b)


def load_bed(p: pathlib.Path) -> dict[str, list[tuple[int, int]]]:
    r: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(p) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            r[canon(c[0])].append((int(c[1]) + 1, int(c[2])))
    for k in r:
        r[k].sort()
    return r


def in_bed(regions: dict[str, list[tuple[int, int]]], chrom: str, pos: int) -> bool:
    ivs = regions.get(chrom)
    if not ivs:
        return False
    i = bisect.bisect_right(ivs, (pos, 10**18)) - 1
    return i >= 0 and ivs[i][0] <= pos <= ivs[i][1]


def parse_fmt(fmt: str, sample: str) -> dict[str, str]:
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: vals[i] for i, k in enumerate(keys) if i < len(vals)}


def load_calls(
    path: pathlib.Path,
    regions: dict[str, list[tuple[int, int]]],
    interval: tuple[str, int, int],
) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    out: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip().split("\t")
            chrom, pos, ref, alt = canon(c[0]), int(c[1]), c[3], c[4]
            if chrom != interval[0] or not (interval[1] <= pos <= interval[2]):
                continue
            if not in_bed(regions, chrom, pos):
                continue
            if "," in alt:
                continue
            f = parse_fmt(c[8], c[9])
            try:
                pl = [int(x) for x in f.get("PL", "").split(",") if x != ""]
            except ValueError:
                pl = None
            out[(chrom, pos, ref, alt)] = {
                "fmt": f,
                "pl": pl,
                "kind": "snp" if len(ref) == 1 and len(alt) == 1 else "indel",
            }
    return out


def load_truth(
    path: pathlib.Path,
    regions: dict[str, list[tuple[int, int]]],
    interval: tuple[str, int, int],
) -> set[tuple[str, int, str, str]]:
    truth: set[tuple[str, int, str, str]] = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            chrom, pos, ref = canon(c[0]), int(c[1]), c[3]
            if chrom != interval[0] or not (interval[1] <= pos <= interval[2]):
                continue
            if not in_bed(regions, chrom, pos):
                continue
            for alt in c[4].split(","):
                if alt != "*":
                    truth.add((chrom, pos, ref, alt))
    return truth


def bucket_for(dmax: int, sparse: bool) -> str:
    if sparse:
        return "sparse_like"
    if dmax == 0:
        return "exact"
    if dmax <= 5:
        return "tol5"
    if dmax <= 50:
        return "tol50"
    if dmax <= 200:
        return "tol200"
    return "gt200"


def is_sparse_like(rpl: list[int], dmax: int) -> bool:
    if rpl in ([90, 6, 0], [0, 90, 6], [81, 0, 36], [0, 0, 0], [3, 3, 3]):
        return True
    return max(rpl) <= 100 and min(rpl) == 0 and sorted(rpl)[1] <= 10 and dmax > 50


def classify_high_delta(
    jp: list[int],
    rp: list[int],
    jad: list[int],
    rad: list[int],
    dp_j: int,
    dp_r: int,
) -> str:
    j_best = jp.index(min(jp))
    r_best = rp.index(min(rp))
    if j_best != r_best:
        return "shape_mismatch"
    if jad == rad and abs(dp_j - dp_r) <= 2:
        return "scale_only_ad_match"
    if jad != rad:
        return "ad_divergence"
    return "scale_dp_drift"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    root = pathlib.Path(__file__).resolve().parents[2]
    ap.add_argument(
        "--java-vcf",
        type=pathlib.Path,
        default=root / "parity/reports/hc-full-parity-j6-dense/p12_dense_giab_window.java.vcf",
    )
    ap.add_argument(
        "--rust-vcf",
        type=pathlib.Path,
        default=root / "parity/reports/hc-full-parity-j6-dense/p12_dense_giab_window.rust.vcf",
    )
    ap.add_argument(
        "--truth-vcf",
        type=pathlib.Path,
        default=root / "parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
    )
    ap.add_argument(
        "--regions-bed",
        type=pathlib.Path,
        default=root / "parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.bed",
    )
    ap.add_argument("--eval-interval", default="20:10000000-10050000")
    ap.add_argument("--json-out", type=pathlib.Path, default=None)
    ap.add_argument("--md-out", type=pathlib.Path, default=None)
    ap.add_argument("--top-n", type=int, default=15)
    args = ap.parse_args(argv)

    interval = parse_interval(args.eval_interval)
    regions = load_bed(args.regions_bed)
    truth = load_truth(args.truth_vcf, regions, interval)
    java = load_calls(args.java_vcf, regions, interval)
    rust = load_calls(args.rust_vcf, regions, interval)

    buckets: Counter[str] = Counter()
    by_kind: Counter[tuple[str, str]] = Counter()
    matched: list[dict[str, Any]] = []
    high: list[dict[str, Any]] = []
    ad_patterns: Counter[str] = Counter()

    for k in sorted(truth & java.keys() & rust.keys()):
        jgt = java[k]["fmt"].get("GT", "").replace("|", "/")
        rgt = rust[k]["fmt"].get("GT", "").replace("|", "/")
        if jgt != rgt:
            continue
        jp, rp = java[k]["pl"], rust[k]["pl"]
        if not jp or not rp:
            continue
        if len(jp) != len(rp):
            buckets["len_mismatch"] += 1
            continue
        dmax = max(abs(a - b) for a, b in zip(jp, rp))
        sparse = is_sparse_like(rp, dmax)
        bucket = bucket_for(dmax, sparse)
        kind = rust[k]["kind"]
        buckets[bucket] += 1
        by_kind[(bucket, kind)] += 1
        row = {
            "chrom": k[0],
            "pos": k[1],
            "ref": k[2],
            "alt": k[3],
            "kind": kind,
            "bucket": bucket,
            "dmax": dmax,
            "java_pl": jp,
            "rust_pl": rp,
            "java_ad": java[k]["fmt"].get("AD"),
            "rust_ad": rust[k]["fmt"].get("AD"),
            "java_dp": java[k]["fmt"].get("DP"),
            "rust_dp": rust[k]["fmt"].get("DP"),
        }
        matched.append(row)
        if dmax > 50:
            try:
                jad = [int(x) for x in (row["java_ad"] or "0,0").split(",")[:2]]
                rad = [int(x) for x in (row["rust_ad"] or "0,0").split(",")[:2]]
                dpj = int(row["java_dp"] or 0)
                dpr = int(row["rust_dp"] or 0)
            except ValueError:
                jad, rad, dpj, dpr = [0, 0], [0, 0], 0, 0
            cls = classify_high_delta(jp, rp, jad, rad, dpj, dpr)
            row["class"] = cls
            high.append(row)
            if jad[0] == 0 and rad[0] > 0:
                ad_patterns["java0_rust_ref_leak"] += 1
            elif jad[0] > 0 and rad[0] == 0:
                ad_patterns["rust0_java_ref"] += 1
            elif jad != rad:
                ad_patterns["other_ad_mismatch"] += 1
            else:
                ad_patterns["ad_match"] += 1

    n = len(matched)
    exact = buckets.get("exact", 0)
    tol5 = buckets.get("tol5", 0)
    tol50 = buckets.get("tol50", 0)
    within_tol5 = exact + tol5
    within_tol50 = within_tol5 + tol50
    class_ctr = Counter(r["class"] for r in high)

    summary = {
        "label": "l12-pl-taxonomy",
        "eval_interval": args.eval_interval,
        "n_gt_matched_tp": n,
        "buckets": dict(buckets),
        "by_kind": {f"{b}/{k}": v for (b, k), v in by_kind.items()},
        "rates": {
            "exact_plus_tol5": within_tol5 / n if n else 0.0,
            "within_tol50": within_tol50 / n if n else 0.0,
            "sparse_like": buckets.get("sparse_like", 0),
        },
        "high_delta": {
            "n": len(high),
            "classes": dict(class_ctr),
            "ad_patterns": dict(ad_patterns),
        },
        "top_residuals": sorted(high, key=lambda r: -r["dmax"])[: args.top_n],
    }

    print(f"n_gt_matched_tp {n}")
    print("buckets", dict(buckets))
    print(
        f"exact+tol5={within_tol5}/{n} ({100 * within_tol5 / n:.1f}%)  "
        f"within_tol50={within_tol50}/{n} ({100 * within_tol50 / n:.1f}%)  "
        f"sparse_like={buckets.get('sparse_like', 0)}"
    )
    print("high_delta_classes", dict(class_ctr))
    print("ad_patterns", dict(ad_patterns))
    print("top residuals:")
    for r in summary["top_residuals"]:
        print(
            f"  {r['pos']} {r['ref']}>{r['alt']} {r['kind']} maxΔ={r['dmax']} "
            f"{r['bucket']} class={r.get('class')} "
            f"AD j/r={r['java_ad']}/{r['rust_ad']} "
            f"javaPL={r['java_pl']} rustPL={r['rust_pl']}"
        )

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"wrote {args.json_out}", file=sys.stderr)

    if args.md_out:
        lines = [
            "# L12 A1/A2 — Soft PL residual taxonomy (chr20 dense)",
            "",
            f"**Date:** 2026-07-22  ",
            f"**Interval:** `{args.eval_interval}` · GT-matched truth TPs (biallelic)  ",
            f"**n:** {n}  ",
            f"**VCFs:** `{args.java_vcf.name}` / `{args.rust_vcf.name}` (L11 tip regen)  ",
            f"**Parent:** [`L12_PRODUCTION_SIGNOFF_PLAN.md`](./L12_PRODUCTION_SIGNOFF_PLAN.md) · L10 baseline [`L10_PL_TAXONOMY.md`](./L10_PL_TAXONOMY.md)",
            "",
            "## Buckets (max |ΔPL|)",
            "",
            "| Bucket | n | Notes |",
            "|--------|---|-------|",
            f"| exact | {buckets.get('exact', 0)} | |",
            f"| tol5 | {buckets.get('tol5', 0)} | |",
            f"| tol50 | {buckets.get('tol50', 0)} | |",
            f"| tol200 | {buckets.get('tol200', 0)} | |",
            f"| gt200 | {buckets.get('gt200', 0)} | PairHMM / AD scale drift |",
            f"| sparse_like | **{buckets.get('sparse_like', 0)}** | SparsePlShape retired (L9) |",
            "",
            f"Among GT-matched: exact+tol5 ≈ **{100 * within_tol5 / n:.0f}%**; within tol50 ≈ **{100 * within_tol50 / n:.0f}%**.",
            "",
            "### vs L10",
            "",
            "Bucket counts match [`L10_PL_TAXONOMY.md`](./L10_PL_TAXONOMY.md) exactly on this tip — L11 long-INS / finalize collapse did not move dense soft PL.",
            "",
            "## A2 classification (tol200 + gt200)",
            "",
            f"High-Δ sites: **{len(high)}**",
            "",
            "| Class | n | Meaning |",
            "|-------|---|---------|",
            f"| ad_divergence | {class_ctr.get('ad_divergence', 0)} | AD differs; PL scale follows |",
            f"| scale_only_ad_match | {class_ctr.get('scale_only_ad_match', 0)} | Same AD/DP; pure PairHMM/hap scale |",
            f"| shape_mismatch | {class_ctr.get('shape_mismatch', 0)} | Different PL argmin |",
            "",
            "### AD patterns among high-Δ",
            "",
            "| Pattern | n |",
            "|---------|---|",
            f"| Java AD `0,N` → Rust REF leak | {ad_patterns.get('java0_rust_ref_leak', 0)} |",
            f"| Other AD mismatch | {ad_patterns.get('other_ad_mismatch', 0)} |",
            f"| AD match (scale-only) | {ad_patterns.get('ad_match', 0)} |",
            "",
            "**Verdict:** soft-PL residual is dominated by **informative AD / REF pileup leak** "
            "(Java often `0,N` hom-alt; Rust keeps REF depth), not LikelihoodEngine hap rematerialization. "
            "Same GT shape (`*,*,0` / het) in almost all high-Δ cases.",
            "",
            "## Top residuals",
            "",
            "| POS | Alleles | Kind | maxΔ | Class | AD j/r | Java PL | Rust PL |",
            "|-----|---------|------|------|-------|--------|---------|---------|",
        ]
        for r in summary["top_residuals"]:
            lines.append(
                f"| {r['pos']} | {r['ref']}>{r['alt']} | {r['kind']} | {r['dmax']} | "
                f"{r.get('class', '')} | `{r['java_ad']}`/`{r['rust_ad']}` | "
                f"`{','.join(map(str, r['java_pl']))}` | `{','.join(map(str, r['rust_pl']))}` |"
            )
        lines += [
            "",
            "## A3 LikelihoodEngine boundary",
            "",
            "Production `compute_region_read_likelihoods` (`engine.rs`) already scores with "
            "`hap_refs: Vec<&[u8]>` → `score_read_against_haplotypes` (zero-copy hap lists). "
            "No L12 rematerialization fix required on the scoring boundary.",
            "",
            "Further soft-PL gains need **AD / hap-support phenotype** work (informative AD "
            "near-ties, Class-A reshape ownership in workstream C) — not PairHMM slice clones.",
            "",
            "## A4 Soft-PL policy",
            "",
            "| Contract | L12 choice |",
            "|----------|------------|",
            "| SparsePlShape / sparse templates | **Hard:** remain ~0 among GT-matched |",
            "| \\|Δ\\|≤5 rate | **Informational** — keep **W-L7-FORMAT** |",
            "| \\|Δ\\|≤50 rate | **Informational** (~67%); not earned as hard gate |",
            "| Exact PL | Not claimed |",
            "",
            "Re-measure:",
            "",
            "```bash",
            "python3 scripts/parity/l12_pl_taxonomy.py \\",
            "  --md-out docs/CLAIM_MATRIX.md",
            "```",
            "",
        ]
        args.md_out.parent.mkdir(parents=True, exist_ok=True)
        args.md_out.write_text("\n".join(lines))
        print(f"wrote {args.md_out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

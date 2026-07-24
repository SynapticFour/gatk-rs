#!/usr/bin/env python3
"""P13: compare Java/Rust callsets to GIAB truth with optional chrom/BED/interval scoping."""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import pathlib
import re

INTERVAL_RE = re.compile(
    r"^\s*(?P<chr>[^:]+)\s*:\s*(?P<s>\d+)\s*-\s*(?P<e>\d+)\s*$"
)


def canon_contig(name: str) -> str:
    n = name.strip()
    if n.startswith("chr"):
        n = n[3:]
    return n


def parse_eval_interval(spec: str | None) -> tuple[str, int, int] | None:
    if not spec or not spec.strip():
        return None
    m = INTERVAL_RE.match(spec.strip())
    if not m:
        return None
    chrom = canon_contig(m.group("chr"))
    start = int(m.group("s"))
    end = int(m.group("e"))
    if end < start:
        return None
    return (chrom, start, end)


def load_regions(path: pathlib.Path | None):
    if path is None or not path.exists():
        return None
    regions: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            chrom = canon_contig(c[0])
            try:
                s = int(c[1]) + 1
                e = int(c[2])
            except ValueError:
                continue
            regions.setdefault(chrom, []).append((s, e))
    for chrom in list(regions.keys()):
        regions[chrom].sort()
    return regions


def in_regions(chrom: str, pos1: int, regions) -> bool:
    """True if 1-based pos1 lies in a BED interval (intervals sorted by start)."""
    if regions is None:
        return True
    ivs = regions.get(chrom)
    if not ivs:
        return False
    # Rightmost interval with start <= pos1 (GIAB high-confidence BED is non-overlapping).
    i = bisect.bisect_right(ivs, (pos1, 10**18)) - 1
    if i < 0:
        return False
    s, e = ivs[i]
    return s <= pos1 <= e


def in_eval_interval(chrom: str, pos1: int, interval: tuple[str, int, int] | None) -> bool:
    if interval is None:
        return True
    ichr, lo, hi = interval
    return chrom == ichr and lo <= pos1 <= hi


def variant_kind(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    return "indel"


def load_vcf(
    path: pathlib.Path,
    chrom_filter: str,
    regions_bed,
    eval_interval: tuple[str, int, int] | None,
) -> set[tuple[str, str, str, str]]:
    if not path.exists():
        return set()
    out: set[tuple[str, str, str, str]] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    ichr, lo, hi = (None, 0, 0)
    if eval_interval is not None:
        ichr, lo, hi = eval_interval
    seen_eval_contig = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 5:
                continue
            chrom = canon_contig(c[0])
            if chrom_filter and chrom != canon_contig(chrom_filter):
                continue
            # Eval interval: skip other contigs and positions outside the window early.
            if ichr is not None:
                if chrom != ichr:
                    # Contig-sorted VCFs: stop once we leave the eval contig.
                    if seen_eval_contig:
                        break
                    continue
                seen_eval_contig = True
                try:
                    pos1 = int(c[1])
                except ValueError:
                    continue
                if pos1 < lo:
                    continue
                if pos1 > hi:
                    break
            else:
                try:
                    pos1 = int(c[1])
                except ValueError:
                    continue
            if not in_regions(chrom, pos1, regions_bed):
                continue
            if not in_eval_interval(chrom, pos1, eval_interval):
                continue
            ref = c[3].strip()
            alt_field = c[4].strip()
            if alt_field in (".", "<NON_REF>"):
                continue
            alts = [
                a
                for a in alt_field.split(",")
                if a and a != "<NON_REF>" and a != "*"
            ]
            if not alts:
                continue
            # L10: expand multi-allelic ALT lists (was first-alt only — under-counted
            # holdout STR sites such as 20:15031984).
            for alt in alts:
                out.add((chrom, str(pos1), ref, alt))
    return out


def metrics(callset: set, truthset: set) -> dict:
    tp = len(callset & truthset)
    fp = len(callset - truthset)
    fn = len(truthset - callset)
    p = tp / (tp + fp) if (tp + fp) else 0.0
    r = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "precision": p, "recall": r, "f1": f1}


def stratified_metrics(
    java: set[tuple[str, str, str, str]],
    rust: set[tuple[str, str, str, str]],
    truth: set[tuple[str, str, str, str]],
) -> dict:
    out: dict = {}
    for kind in ("snp", "indel"):
        jk = {v for v in java if variant_kind(v[2], v[3]) == kind}
        rk = {v for v in rust if variant_kind(v[2], v[3]) == kind}
        tk = {v for v in truth if variant_kind(v[2], v[3]) == kind}
        out[kind] = {
            "truth_count": len(tk),
            "java": metrics(jk, tk),
            "rust": metrics(rk, tk),
        }
    return out


def check_threshold_pair(
    label: str,
    java_m: dict,
    rust_m: dict,
    cfg: dict,
    failures: list[str],
) -> None:
    java_f1 = java_m.get("f1", 0.0)
    rust_f1 = rust_m.get("f1", 0.0)
    min_frac = cfg.get("rust_min_f1_fraction_of_java")
    if min_frac is not None and java_f1 > 0.0 and rust_f1 < java_f1 * float(min_frac):
        failures.append(
            f"{label}: rust_f1 {rust_f1:.6f} < {min_frac} * java_f1 {java_f1:.6f}"
        )
    max_delta = cfg.get("rust_java_f1_max_delta")
    if max_delta is not None and rust_f1 < java_f1 - float(max_delta):
        failures.append(
            f"{label}: rust_f1 {rust_f1:.6f} < java_f1 {java_f1:.6f} - {max_delta}"
        )
    min_rust = cfg.get("rust_min_f1")
    if min_rust is not None and rust_f1 < float(min_rust):
        failures.append(f"{label}: rust_f1 {rust_f1:.6f} < min {min_rust}")


def evaluate_gate(payload: dict, thresholds: dict | None) -> tuple[str, list[str]]:
    if not thresholds:
        return "skipped", []
    if payload.get("status") == "truth_missing":
        return "skipped", ["truth_missing"]
    if payload.get("status") == "truth_empty":
        return "fail", ["truth_empty_in_eval_scope"]
    truth_n = payload.get("truth_variant_count", 0)
    min_truth = int(thresholds.get("min_truth_variants", 1))
    if truth_n < min_truth:
        return "fail", [f"truth_variant_count {truth_n} < min {min_truth}"]

    failures: list[str] = []
    check_threshold_pair("overall", payload["java"], payload["rust"], thresholds, failures)

    strat_cfg = thresholds.get("stratified") or {}
    stratified = payload.get("stratified") or {}
    for kind, cfg in strat_cfg.items():
        block = stratified.get(kind)
        if not block or block.get("truth_count", 0) == 0:
            continue
        check_threshold_pair(
            kind,
            block["java"],
            block["rust"],
            {**thresholds, **cfg},
            failures,
        )
    return ("pass" if not failures else "fail"), failures


def write_truth_missing(json_out: pathlib.Path, md_out: pathlib.Path) -> None:
    payload = {
        "label": "phase13-realworld-truth-eval",
        "status": "truth_missing",
        "gate_status": "skipped",
        "notes": "Set P13_TRUTH_VCF to a GIAB truth VCF to evaluate Java and Rust callsets against external truth.",
    }
    json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_out.write_text(
        "# P13 Real-world Truth Eval\n\n- status: **truth_missing**\n- next: export `P13_TRUTH_VCF=/path/to/giab_truth.vcf.gz` and rerun\n",
        encoding="utf-8",
    )
    print("[p13-truth] truth VCF missing; evaluation skipped")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--java-vcf", type=pathlib.Path, required=True)
    p.add_argument("--rust-vcf", type=pathlib.Path, required=True)
    p.add_argument("--truth-vcf", type=pathlib.Path, required=True)
    p.add_argument("--json-out", type=pathlib.Path, required=True)
    p.add_argument("--md-out", type=pathlib.Path, required=True)
    p.add_argument("--chrom-filter", default="")
    p.add_argument("--regions-bed", type=pathlib.Path, default=None)
    p.add_argument(
        "--eval-interval",
        default="",
        help="GATK-style chrom:start-end (1-based inclusive), e.g. 20:100-50000. "
        "When set, truth and calls are restricted to this window (after BED filter).",
    )
    p.add_argument(
        "--thresholds-json",
        type=pathlib.Path,
        default=None,
        help="L6 gate thresholds (see parity/fixtures/hc-full-parity/j6/thresholds.json)",
    )
    p.add_argument(
        "--strict-gate",
        action="store_true",
        help="Exit non-zero when gate_status is fail",
    )
    args = p.parse_args()

    regions = load_regions(args.regions_bed)
    eval_iv = parse_eval_interval(args.eval_interval.strip() or None)

    truth = load_vcf(args.truth_vcf, args.chrom_filter.strip(), regions, eval_iv)
    java = load_vcf(args.java_vcf, args.chrom_filter.strip(), regions, eval_iv)
    rust = load_vcf(args.rust_vcf, args.chrom_filter.strip(), regions, eval_iv)

    status = "pass" if truth else "truth_empty"
    stratified = stratified_metrics(java, rust, truth)
    thresholds = None
    if args.thresholds_json and args.thresholds_json.is_file():
        thresholds = json.loads(args.thresholds_json.read_text(encoding="utf-8"))

    payload = {
        "label": "phase13-realworld-truth-eval",
        "status": status,
        "truth_variant_count": len(truth),
        "eval_interval": args.eval_interval.strip() or None,
        "java": metrics(java, truth),
        "rust": metrics(rust, truth),
        "stratified": stratified,
        "java_vcf": str(args.java_vcf),
        "rust_vcf": str(args.rust_vcf),
        "truth_vcf": str(args.truth_vcf),
        "regions_bed": str(args.regions_bed) if args.regions_bed else None,
        "thresholds_json": str(args.thresholds_json) if args.thresholds_json else None,
    }
    gate_status, gate_failures = evaluate_gate(payload, thresholds)
    payload["gate_status"] = gate_status
    payload["gate_failures"] = gate_failures

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    iv_note = payload["eval_interval"] or "(none — whole chromosome in BED scope)"
    md_lines = [
        "# P13 Real-world Truth Eval",
        "",
        f"- status: **{status}**",
        f"- gate_status: **{gate_status}**",
        f"- eval_interval: `{iv_note}`",
        f"- truth variants (eval scope): `{len(truth)}`",
        f"- java f1: `{payload['java']['f1']:.6f}` (P={payload['java']['precision']:.6f}, R={payload['java']['recall']:.6f})",
        f"- rust f1: `{payload['rust']['f1']:.6f}` (P={payload['rust']['precision']:.6f}, R={payload['rust']['recall']:.6f})",
    ]
    for kind in ("snp", "indel"):
        block = stratified[kind]
        j = block["java"]
        r = block["rust"]
        md_lines.append(
            f"- {kind}: truth `{block['truth_count']}` — "
            f"java F1 `{j['f1']:.6f}` rust F1 `{r['f1']:.6f}`"
        )
    if gate_failures:
        md_lines.append(f"- gate_failures: `{'; '.join(gate_failures)}`")
    args.md_out.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(
        f"[p13-truth] status={status} gate={gate_status} "
        f"java_f1={payload['java']['f1']:.6f} rust_f1={payload['rust']['f1']:.6f}"
    )
    if args.strict_gate and gate_status == "fail":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

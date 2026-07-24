#!/usr/bin/env python3
"""Summarize L6 gate: P12 scale parity + P13 truth stratification.

Variant-set parity is evaluated on ``parity_interval`` (default: P12 spine
``2:92300000-92350000``). The HC / truth eval window may be wider so GIAB
high-confidence sites exist in scope — see ``thresholds.json``.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

INTERVAL_RE = re.compile(
    r"^\s*(?P<chr>[^:]+)\s*:\s*(?P<s>\d+)\s*-\s*(?P<e>\d+)\s*$"
)


def canon_contig(name: str) -> str:
    n = name.strip()
    if n.startswith("chr"):
        n = n[3:]
    return n


def parse_interval(spec: str | None) -> tuple[str, int, int] | None:
    if not spec or not str(spec).strip():
        return None
    m = INTERVAL_RE.match(str(spec).strip())
    if not m:
        return None
    chrom = canon_contig(m.group("chr"))
    start = int(m.group("s"))
    end = int(m.group("e"))
    if end < start:
        return None
    return (chrom, start, end)


def load_variants(path: pathlib.Path) -> set[tuple[str, str, str, str]]:
    if not path.exists():
        return set()
    out: set[tuple[str, str, str, str]] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) < 5:
            continue
        alt = c[4].split(",")[0]
        out.add((canon_contig(c[0]), c[1], c[3], alt))
    return out


def filter_interval(
    keys: set[tuple[str, str, str, str]], interval: tuple[str, int, int] | None
) -> set[tuple[str, str, str, str]]:
    if interval is None:
        return keys
    chrom, lo, hi = interval
    out = set()
    for c, pos, ref, alt in keys:
        try:
            p = int(pos)
        except ValueError:
            continue
        if c == chrom and lo <= p <= hi:
            out.add((c, pos, ref, alt))
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--p12-json", type=pathlib.Path, required=True)
    p.add_argument("--p13-json", type=pathlib.Path, required=True)
    p.add_argument("--json-out", type=pathlib.Path, required=True)
    p.add_argument("--md-out", type=pathlib.Path, required=True)
    p.add_argument("--interval", default="")
    p.add_argument(
        "--parity-interval",
        default="",
        help="Interval for variant-set parity (default from thresholds or P12 spine)",
    )
    p.add_argument("--thresholds", type=pathlib.Path, default=None)
    args = p.parse_args()

    p12 = json.loads(args.p12_json.read_text(encoding="utf-8"))
    p13 = json.loads(args.p13_json.read_text(encoding="utf-8"))

    thresholds: dict = {}
    if args.thresholds and args.thresholds.is_file():
        thresholds = json.loads(args.thresholds.read_text(encoding="utf-8"))

    parity_spec = (
        args.parity_interval
        or thresholds.get("parity_interval_default")
        or "2:92300000-92350000"
    )
    parity_iv = parse_interval(parity_spec)

    java_vcf = pathlib.Path(p12.get("java_vcf") or "")
    rust_vcf = pathlib.Path(p12.get("rust_vcf") or "")
    j_all = load_variants(java_vcf)
    r_all = load_variants(rust_vcf)
    j_par = filter_interval(j_all, parity_iv)
    r_par = filter_interval(r_all, parity_iv)
    spine_shared = len(j_par & r_par)
    spine_parity = (
        "variant_parity"
        if p12.get("java_exit", 0) == 0
        and p12.get("rust_exit", 0) == 0
        and j_par == r_par
        and len(j_par) > 0
        else (
            "tool_error"
            if p12.get("java_exit", 0) != 0 or p12.get("rust_exit", 0) != 0
            else "variant_mismatch"
        )
    )
    # Empty both on spine with successful exits is still mismatch for L6.
    if (
        p12.get("java_exit", 0) == 0
        and p12.get("rust_exit", 0) == 0
        and j_par == r_par
        and len(j_par) == 0
    ):
        spine_parity = "variant_mismatch"

    blockers: list[str] = []
    p12_status = p12.get("status", "unknown")
    p13_status = p13.get("status", "unknown")
    p13_gate = p13.get("gate_status", p13_status)
    window_parity = p12.get("parity_status")

    if p12_status not in ("pass",):
        blockers.append(f"p12_status={p12_status}")

    if thresholds.get("p12_parity_required", True) and spine_parity != "variant_parity":
        blockers.append(f"spine_parity={spine_parity}")

    if p13_gate not in ("pass", "skipped"):
        blockers.append(f"p13_gate={p13_gate}")
    if p13.get("gate_failures"):
        blockers.extend(str(x) for x in p13["gate_failures"])

    if p13_status == "truth_missing":
        overall = "assets_missing"
    elif blockers:
        overall = "fail"
    else:
        overall = "pass"

    java = p13.get("java") or {}
    rust = p13.get("rust") or {}
    notes: list[str] = []
    if (
        float(java.get("f1") or 0) == 0.0
        and float(rust.get("f1") or 0) == 0.0
        and int(p13.get("truth_variant_count") or 0) > 0
    ):
        notes.append(
            "truth_f1_vacuous_both_zero — NA12878_20k coverage does not support "
            "GIAB TP calls in this window; rust tracks java (both F1=0)"
        )

    payload = {
        "label": "hc-full-parity-j6-truth",
        "status": overall,
        "eval_interval": args.interval or thresholds.get("eval_interval_default"),
        "parity_interval": parity_spec,
        "p12_status": p12_status,
        "p12_parity_status_window": window_parity,
        "p12_parity_status": spine_parity,
        "p13_status": p13_status,
        "p13_gate_status": p13_gate,
        "blockers": blockers,
        "notes": notes,
        "p12": {
            "java_variant_count": p12.get("java_variant_count"),
            "rust_variant_count": p12.get("rust_variant_count"),
            "shared_variant_count": p12.get("shared_variant_count"),
            "java_only": p12.get("java_only"),
            "rust_only": p12.get("rust_only"),
            "spine_java_variant_count": len(j_par),
            "spine_rust_variant_count": len(r_par),
            "spine_shared_variant_count": spine_shared,
            "spine_java_only": len(j_par - r_par),
            "spine_rust_only": len(r_par - j_par),
        },
        "p13": {
            "truth_variant_count": p13.get("truth_variant_count"),
            "java": p13.get("java"),
            "rust": p13.get("rust"),
            "stratified": p13.get("stratified"),
        },
        "thresholds": str(args.thresholds) if args.thresholds else None,
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    md_lines = [
        "# L6 HC full parity — scale + GIAB truth",
        "",
        f"- status: **{overall}**",
        f"- eval_interval (HC + truth): `{payload['eval_interval']}`",
        f"- parity_interval (variant-set lock): `{parity_spec}`",
        f"- p12 harness: `{p12_status}` / window parity `{window_parity}` / spine parity `{spine_parity}`",
        f"- spine variants: java `{len(j_par)}` rust `{len(r_par)}` shared `{spine_shared}`",
        f"- p13 gate: `{p13_gate}` (status `{p13_status}`)",
    ]
    if blockers:
        md_lines.append(f"- blockers: `{'; '.join(blockers)}`")
    if notes:
        md_lines.append(f"- notes: `{'; '.join(notes)}`")
    md_lines.extend(
        [
            f"- truth variants (scope): `{p13.get('truth_variant_count', 0)}`",
            f"- java F1: `{java.get('f1', 0):.6f}`",
            f"- rust F1: `{rust.get('f1', 0):.6f}`",
        ]
    )
    strat = p13.get("stratified") or {}
    for kind in ("snp", "indel"):
        if kind in strat:
            j = strat[kind].get("java", {})
            r = strat[kind].get("rust", {})
            md_lines.append(
                f"- {kind} F1: java `{j.get('f1', 0):.6f}` rust `{r.get('f1', 0):.6f}` "
                f"(truth n={strat[kind].get('truth_count', 0)})"
            )
    args.md_out.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(
        f"[j6-truth] status={overall} p12={p12_status} "
        f"spine_parity={spine_parity} p13_gate={p13_gate}"
    )
    return 0 if overall == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

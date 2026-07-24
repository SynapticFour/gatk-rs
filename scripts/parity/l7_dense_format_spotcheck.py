#!/usr/bin/env python3
"""L7-A4: FORMAT spot-check (GT/GQ/DP/AD/PL) on dense true-positive sites.

Compares Rust vs Java on sites where both callsets match truth (first-ALT keys)
inside the eval interval and optional high-confidence BED.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
from typing import Iterable

INTERVAL_RE = re.compile(r"^\s*(?P<chr>[^:]+)\s*:\s*(?P<s>\d+)\s*-\s*(?P<e>\d+)\s*$")


def canon_contig(name: str) -> str:
    n = name.strip()
    return n[3:] if n.startswith("chr") else n


def parse_interval(spec: str) -> tuple[str, int, int]:
    m = INTERVAL_RE.match(spec.strip())
    if not m:
        raise SystemExit(f"bad interval: {spec!r}")
    return canon_contig(m.group("chr")), int(m.group("s")), int(m.group("e"))


def load_bed(path: pathlib.Path | None) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    if path is None or not path.is_file():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#") or line.startswith("track"):
            continue
        c = line.split()
        if len(c) < 3:
            continue
        chrom = canon_contig(c[0])
        out.setdefault(chrom, []).append((int(c[1]), int(c[2])))
    for chrom in out:
        out[chrom].sort()
    return out


def in_bed(bed: dict[str, list[tuple[int, int]]], chrom: str, pos1: int) -> bool:
    ivs = bed.get(chrom)
    if not ivs:
        return True
    pos0 = pos1 - 1
    lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos0 < s:
            hi = mid - 1
        elif pos0 >= e:
            lo = mid + 1
        else:
            return True
    return False


def parse_format(fmt: str, sample: str) -> dict[str, str]:
    keys = fmt.split(":")
    vals = sample.split(":")
    if len(keys) != len(vals):
        return {}
    return dict(zip(keys, vals))


def load_records(
    path: pathlib.Path,
    interval: tuple[str, int, int],
    bed: dict[str, list[tuple[int, int]]],
) -> dict[tuple[str, int, str, str], dict]:
    chrom, lo, hi = interval
    out: dict[tuple[str, int, str, str], dict] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 10:
            continue
        c = canon_contig(cols[0])
        pos = int(cols[1])
        if c != chrom or pos < lo or pos > hi:
            continue
        if not in_bed(bed, c, pos):
            continue
        ref = cols[3]
        alt0 = cols[4].split(",")[0]
        key = (c, pos, ref, alt0)
        out[key] = {
            "alts": cols[4],
            "qual": cols[5],
            "filter": cols[6],
            "info": cols[7],
            "format": parse_format(cols[8], cols[9]),
            "raw_format": cols[8],
            "raw_sample": cols[9],
        }
    return out


def load_truth_keys(
    path: pathlib.Path,
    interval: tuple[str, int, int],
    bed: dict[str, list[tuple[int, int]]],
) -> set[tuple[str, int, str, str]]:
    import gzip

    chrom, lo, hi = interval
    opener = gzip.open if str(path).endswith(".gz") else open
    keys: set[tuple[str, int, str, str]] = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            c = canon_contig(cols[0])
            pos = int(cols[1])
            if c != chrom or pos < lo or pos > hi:
                continue
            if not in_bed(bed, c, pos):
                continue
            ref = cols[3]
            alt0 = cols[4].split(",")[0]
            keys.add((c, pos, ref, alt0))
    return keys


def _cap_gq(v: str | None) -> str | None:
    if v is None:
        return None
    try:
        return str(min(int(v), 99))
    except ValueError:
        return v


def _ad_l1(rust_ad: str | None, java_ad: str) -> int | None:
    try:
        r = [int(x) for x in (rust_ad or "").split(",") if x != ""]
        j = [int(x) for x in java_ad.split(",") if x != ""]
    except ValueError:
        return None
    if len(r) < 2 or len(j) < 2:
        return None
    return abs(r[0] - j[0]) + abs(r[1] - j[1])


def compare_format(
    rust: dict[str, str],
    java: dict[str, str],
    pl_tol: int,
    *,
    hard_keys: tuple[str, ...],
    soft_keys: tuple[str, ...] = ("AD", "DP", "PL"),
    ad_l1_tol: int = 2,
    dp_tol: int = 2,
) -> tuple[list[str], list[str]]:
    """Return (hard_errors, soft_errors). Soft keys default to AD/DP/PL."""
    hard: list[str] = []
    soft: list[str] = []
    soft_set = set(soft_keys)
    if "GT" in hard_keys and "GT" in java and rust.get("GT") != java.get("GT"):
        hard.append(f"GT: rust={rust.get('GT')!r} java={java.get('GT')!r}")
    if "GQ" in hard_keys and "GQ" in java:
        rg = _cap_gq(rust.get("GQ"))
        jg = _cap_gq(java.get("GQ"))
        if rg != jg:
            hard.append(f"GQ: rust={rust.get('GQ')!r}->{rg!r} java={java.get('GQ')!r}")
    if "DP" in soft_set and "DP" in java:
        try:
            rd = int(rust.get("DP", ""))
            jd = int(java["DP"])
            if abs(rd - jd) > dp_tol:
                soft.append(f"DP: rust={rust.get('DP')!r} java={java.get('DP')!r}")
        except ValueError:
            if rust.get("DP") != java.get("DP"):
                soft.append(f"DP: rust={rust.get('DP')!r} java={java.get('DP')!r}")
    if "AD" in soft_set and "AD" in java:
        l1 = _ad_l1(rust.get("AD"), java["AD"])
        if l1 is None or l1 > ad_l1_tol:
            soft.append(f"AD: rust={rust.get('AD')!r} java={java['AD']!r}")
    if "PL" in soft_set and "PL" in java:
        try:
            r_pl = [int(x) for x in rust.get("PL", "").split(",") if x != ""]
            j_pl = [int(x) for x in java["PL"].split(",") if x != ""]
        except ValueError:
            soft.append("PL: parse error")
        else:
            if len(r_pl) != len(j_pl):
                soft.append(f"PL length: rust={len(r_pl)} java={len(j_pl)}")
            else:
                for i, (rp, jp) in enumerate(zip(r_pl, j_pl)):
                    if abs(rp - jp) > pl_tol:
                        soft.append(f"PL[{i}]: rust={rp} java={jp}")
    return hard, soft


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--java-vcf", type=pathlib.Path, required=True)
    ap.add_argument("--rust-vcf", type=pathlib.Path, required=True)
    ap.add_argument("--truth-vcf", type=pathlib.Path, required=True)
    ap.add_argument("--regions-bed", type=pathlib.Path, default=None)
    ap.add_argument("--eval-interval", required=True)
    ap.add_argument("--json-out", type=pathlib.Path, required=True)
    ap.add_argument("--md-out", type=pathlib.Path, required=True)
    ap.add_argument("--max-sites", type=int, default=40)
    ap.add_argument("--pl-tol", type=int, default=5, help="Soft PL absolute tolerance (L8)")
    ap.add_argument("--ad-l1-tol", type=int, default=2, help="Soft AD L1 tolerance (L8)")
    ap.add_argument("--dp-tol", type=int, default=2, help="Soft DP absolute tolerance (L8)")
    ap.add_argument(
        "--max-hard-mismatch-rate",
        type=float,
        default=0.15,
        help="Fail if hard (GT + GQ≤99) mismatch rate among checked TPs exceeds this",
    )
    ap.add_argument(
        "--max-soft-mismatch-rate",
        type=float,
        default=None,
        help=(
            "L8: if set, fail when soft mismatch rate among GT-matched TPs exceeds this "
            "(soft keys default AD,DP,PL; use --soft-keys AD,DP to exclude PL)"
        ),
    )
    ap.add_argument(
        "--hard-keys",
        default="GT",
        help="Comma-separated hard FORMAT keys (L7 default: GT; GQ/AD/PL residual tracked soft)",
    )
    ap.add_argument(
        "--soft-keys",
        default="AD,DP,PL",
        help="Comma-separated soft FORMAT keys gated by tolerances (L8 AD/DP gate: AD,DP)",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    interval = parse_interval(args.eval_interval)
    bed = load_bed(args.regions_bed)
    truth = load_truth_keys(args.truth_vcf, interval, bed)
    java = load_records(args.java_vcf, interval, bed)
    rust = load_records(args.rust_vcf, interval, bed)
    tps = sorted(truth & set(java) & set(rust))
    checked = tps[: max(0, args.max_sites)]
    hard_keys = tuple(k.strip() for k in args.hard_keys.split(",") if k.strip())
    soft_keys = tuple(k.strip() for k in args.soft_keys.split(",") if k.strip())

    hard_mismatches: list[dict] = []
    soft_mismatches: list[dict] = []
    hard_matches = 0
    soft_matches = 0
    for key in checked:
        chrom, pos, ref, alt = key
        hard, soft = compare_format(
            rust[key]["format"],
            java[key]["format"],
            args.pl_tol,
            hard_keys=hard_keys,
            soft_keys=soft_keys,
            ad_l1_tol=args.ad_l1_tol,
            dp_tol=args.dp_tol,
        )
        site = f"{chrom}:{pos}:{ref}>{alt}"
        if hard:
            hard_mismatches.append(
                {
                    "site": site,
                    "errors": hard,
                    "soft": soft,
                    "rust": rust[key]["format"],
                    "java": java[key]["format"],
                }
            )
        else:
            hard_matches += 1
            # Soft residual is tracked among GT-matched TPs only (L8).
            if soft:
                soft_mismatches.append({"site": site, "errors": soft})
            else:
                soft_matches += 1

    n = len(checked)
    hard_rate = (len(hard_mismatches) / n) if n else 1.0
    soft_denom = max(hard_matches, 1)
    soft_rate = (len(soft_mismatches) / soft_denom) if hard_matches else 1.0
    hard_ok = n > 0 and hard_rate <= args.max_hard_mismatch_rate
    soft_ok = True
    if args.max_soft_mismatch_rate is not None:
        soft_ok = hard_matches > 0 and soft_rate <= args.max_soft_mismatch_rate
    status = "pass" if hard_ok and soft_ok else "fail"

    summary = {
        "label": "l7-a4-dense-format-spotcheck",
        "status": status,
        "eval_interval": args.eval_interval,
        "truth_tp_shared": len(tps),
        "checked": n,
        "hard_keys": list(hard_keys),
        "soft_keys": list(soft_keys),
        "hard_match": hard_matches,
        "hard_mismatch": len(hard_mismatches),
        "hard_mismatch_rate": hard_rate,
        "max_hard_mismatch_rate": args.max_hard_mismatch_rate,
        "soft_match": soft_matches,
        "soft_mismatch": len(soft_mismatches),
        "soft_mismatch_rate": soft_rate,
        "max_soft_mismatch_rate": args.max_soft_mismatch_rate,
        "pl_tol": args.pl_tol,
        "ad_l1_tol": args.ad_l1_tol,
        "dp_tol": args.dp_tol,
        "hard_mismatches": hard_mismatches[:20],
        "soft_mismatch_examples": soft_mismatches[:10],
        "notes": (
            "Hard gate: GT (+ optional GQ). Soft keys among GT-matched TPs only; "
            "L8 tolerances AD L1≤ad_l1_tol, |DP|≤dp_tol, |PL[i]|≤pl_tol. "
            "Exact PL identity remains W-L7-FORMAT / L9 when soft-keys omit PL."
        ),
        "java_vcf": str(args.java_vcf),
        "rust_vcf": str(args.rust_vcf),
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# L7-A4 dense FORMAT spot-check",
        "",
        f"- status: `{status}`",
        f"- interval: `{args.eval_interval}`",
        f"- shared truth TPs (first-ALT): `{len(tps)}`",
        f"- checked: `{n}` (cap `{args.max_sites}`)",
        f"- hard keys: `{','.join(hard_keys)}` (GQ capped at 99)",
        f"- hard match: `{hard_matches}` mismatch `{len(hard_mismatches)}` "
        f"(rate `{hard_rate:.3f}`, max `{args.max_hard_mismatch_rate}`)",
        f"- soft keys: `{','.join(soft_keys)}`",
        f"- soft among GT-matched: `{soft_matches}` mismatch `{len(soft_mismatches)}` "
        f"(rate `{soft_rate:.3f}`"
        + (
            f", max `{args.max_soft_mismatch_rate}`)"
            if args.max_soft_mismatch_rate is not None
            else ", informational)"
        ),
        f"- soft tols: AD L1≤`{args.ad_l1_tol}`, |DP|≤`{args.dp_tol}`, |PL[i]|≤`{args.pl_tol}`",
        "",
    ]
    if hard_mismatches:
        lines.append("## Hard mismatch examples")
        for row in hard_mismatches[:10]:
            lines.append(f"- `{row['site']}`: {'; '.join(row['errors'])}")
        lines.append("")
    if soft_mismatches:
        lines.append("## Soft mismatch examples (AD/DP/PL)")
        for row in soft_mismatches[:5]:
            lines.append(f"- `{row['site']}`: {'; '.join(row['errors'])}")
        lines.append("")
    args.md_out.write_text("\n".join(lines), encoding="utf-8")

    soft_gate = (
        f" soft_max={args.max_soft_mismatch_rate}"
        if args.max_soft_mismatch_rate is not None
        else ""
    )
    print(
        f"[l7-a4-format] status={status} checked={n} "
        f"hard_match={hard_matches} hard_mismatch={len(hard_mismatches)} "
        f"hard_rate={hard_rate:.3f} soft_rate={soft_rate:.3f} "
        f"soft_keys={','.join(soft_keys)}{soft_gate}"
    )
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())

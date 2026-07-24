#!/usr/bin/env python3
"""Standalone Java↔Rust callset comparison for external pilots.

No dependency on gatk-rs CI, Docker helpers, or internal parity fixtures.
Works with callsets you already produced (your Java pipeline + gatk-rs).

Modes
-----
1. **Direct diff (no truth)** — site / allele / GT identity + FORMAT drift classes.
2. **Truth eval** — optional hap.py (preferred) or RTG vcfeval against your truth
   VCF (+ confident BED), then reports Rust−Java ΔF1.

Exit codes
----------
  0  hard gates passed (soft FORMAT drift may still be listed)
  1  hard mismatch (missing sites, GT/allele/FILTER disagreement beyond tolerances)
  2  configuration / tool error (e.g. truth mode requested but hap.py/rtg missing)

Examples
--------
  # Pure Rust-vs-Java (no truth set):
  python3 scripts/pilot/compare_callsets.py \\
    --java java.genotyped.vcf --rust rust.genotyped.vcf \\
    --out pilot_out

  # With your truth + confident regions:
  python3 scripts/pilot/compare_callsets.py \\
    --java java.vcf --rust rust.vcf \\
    --reference hs37d5.fa \\
    --truth HG001_benchmark.vcf.gz --confident HG001_benchmark.bed \\
    --out pilot_out --engine auto
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple


# ---------------------------------------------------------------------------
# VCF I/O
# ---------------------------------------------------------------------------


def open_text(path: Path):
    name = path.name.lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def canon_chrom(c: str) -> str:
    t = c.strip()
    return t[3:] if t.startswith("chr") or t.startswith("CHR") else t


def normalize_gt(gt: str) -> str:
    return gt.replace("|", "/")


@dataclass
class SampleFmt:
    gt: str
    fields: Dict[str, str] = field(default_factory=dict)


@dataclass
class Site:
    chrom: str
    pos: int
    ref: str
    alts: List[str]
    qual: Optional[float]
    filt: str
    samples: Dict[str, SampleFmt]


Key = Tuple[str, int]  # chrom, pos


def parse_vcf(path: Path) -> Tuple[List[str], Dict[Key, Site]]:
    samples: List[str] = []
    sites: Dict[Key, Site] = {}
    with open_text(path) as fh:
        for ln in fh:
            if not ln or ln.startswith("##"):
                continue
            if ln.startswith("#CHROM"):
                samples = ln.rstrip("\n").split("\t")[9:]
                continue
            if ln.startswith("#"):
                continue
            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom = canon_chrom(cols[0])
            pos = int(cols[1])
            ref = cols[3].upper()
            alts = [a.upper() for a in cols[4].split(",") if a and a != "."]
            qual: Optional[float]
            try:
                qual = None if cols[5] == "." else float(cols[5])
            except ValueError:
                qual = None
            filt = cols[6]
            fmt_keys = cols[8].split(":") if len(cols) > 8 else []
            smap: Dict[str, SampleFmt] = {}
            for i, name in enumerate(samples):
                if 9 + i >= len(cols):
                    break
                parts = cols[9 + i].split(":")
                fmap = {
                    k: (parts[j] if j < len(parts) else ".")
                    for j, k in enumerate(fmt_keys)
                }
                smap[name] = SampleFmt(gt=normalize_gt(fmap.get("GT", ".")), fields=fmap)
            sites[(chrom, pos)] = Site(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alts=alts,
                qual=qual,
                filt=filt,
                samples=smap,
            )
    return samples, sites


def allele_set(site: Site) -> Set[str]:
    return set(site.alts)


def called_alleles(site: Site, sample: str) -> Optional[frozenset]:
    """Unordered bases implied by GT (handles ALT-order differences)."""
    sf = site.samples.get(sample)
    if sf is None:
        return None
    gt = sf.gt
    if gt in (".", "./.", ".|."):
        return frozenset()
    alleles = [site.ref] + site.alts
    out = []
    for tok in gt.split("/"):
        if tok in (".", ""):
            continue
        try:
            idx = int(tok)
        except ValueError:
            return None
        if idx < 0 or idx >= len(alleles):
            return None
        out.append(alleles[idx])
    return frozenset(out)


# ---------------------------------------------------------------------------
# Direct compare
# ---------------------------------------------------------------------------


def _parse_int_list(s: str) -> Optional[List[int]]:
    if not s or s == ".":
        return None
    try:
        return [int(x) for x in s.split(",")]
    except ValueError:
        return None


def format_soft_ok(
    java_f: Dict[str, str],
    rust_f: Dict[str, str],
    *,
    ad_tol: float,
    dp_tol: int,
    gq_tol: int,
    pl_tol: int,
) -> Tuple[bool, List[str]]:
    """Return (ok, reasons). Soft FORMAT drift that pilots may ignore."""
    reasons: List[str] = []
    # AD (relative)
    jad, rad = _parse_int_list(java_f.get("AD", ".")), _parse_int_list(rust_f.get("AD", "."))
    if jad is not None and rad is not None and len(jad) == len(rad):
        for a, b in zip(jad, rad):
            denom = max(a, b, 1)
            if abs(a - b) / denom > ad_tol and abs(a - b) > 1:
                reasons.append(f"AD java={jad} rust={rad}")
                break
    elif jad != rad and (jad is not None or rad is not None):
        # Length/presence mismatch — still soft if both missing-ish
        if not (jad is None and rad is None):
            reasons.append(f"AD java={java_f.get('AD')} rust={rust_f.get('AD')}")

    def _int(v: str) -> Optional[int]:
        try:
            return int(v)
        except (TypeError, ValueError):
            return None

    jdp, rdp = _int(java_f.get("DP", ".")), _int(rust_f.get("DP", "."))
    if jdp is not None and rdp is not None and abs(jdp - rdp) > dp_tol:
        reasons.append(f"DP java={jdp} rust={rdp}")

    jgq, rgq = _int(java_f.get("GQ", ".")), _int(rust_f.get("GQ", "."))
    if jgq is not None and rgq is not None and abs(jgq - rgq) > gq_tol:
        reasons.append(f"GQ java={jgq} rust={rgq}")

    jpl, rpl = _parse_int_list(java_f.get("PL", ".")), _parse_int_list(rust_f.get("PL", "."))
    if jpl is not None and rpl is not None:
        if len(jpl) != len(rpl):
            reasons.append(f"PL_len java={len(jpl)} rust={len(rpl)}")
        else:
            for a, b in zip(jpl, rpl):
                if abs(a - b) > pl_tol:
                    reasons.append(f"PL java={jpl} rust={rpl}")
                    break
    return (len(reasons) == 0), reasons


def compare_direct(
    java_path: Path,
    rust_path: Path,
    *,
    qual_tol: float,
    ad_tol: float,
    dp_tol: int,
    gq_tol: int,
    pl_tol: int,
    ignore_filter: bool,
) -> dict:
    j_samples, j_sites = parse_vcf(java_path)
    r_samples, r_sites = parse_vcf(rust_path)
    keys = sorted(set(j_sites) | set(r_sites))

    only_java = 0
    only_rust = 0
    allele_mismatch = 0
    gt_mismatch = 0
    filter_mismatch = 0
    qual_mismatch = 0
    format_soft = 0
    exact_match = 0
    examples: Dict[str, List[str]] = defaultdict(list)

    def note(kind: str, msg: str, limit: int = 25) -> None:
        if len(examples[kind]) < limit:
            examples[kind].append(msg)

    common_samples = sorted(set(j_samples) & set(r_samples))
    if not common_samples and (j_samples or r_samples):
        # Fall back to positional names if headers differ
        n = min(len(j_samples), len(r_samples))
        common_samples = [f"__idx{i}__" for i in range(n)]

    for key in keys:
        if key not in j_sites:
            only_rust += 1
            note("only_rust", f"{key[0]}:{key[1]}")
            continue
        if key not in r_sites:
            only_java += 1
            note("only_java", f"{key[0]}:{key[1]}")
            continue
        js, rs = j_sites[key], r_sites[key]
        hard_here = False
        if js.ref != rs.ref or allele_set(js) != allele_set(rs):
            allele_mismatch += 1
            hard_here = True
            note(
                "allele",
                f"{key[0]}:{key[1]} java={js.ref}/{sorted(allele_set(js))} "
                f"rust={rs.ref}/{sorted(allele_set(rs))}",
            )
        if not ignore_filter:
            jf = set(js.filt.replace(";", ",").split(",")) - {".", "PASS", ""}
            rf = set(rs.filt.replace(";", ",").split(",")) - {".", "PASS", ""}
            # PASS vs . treated as equal for soft-filter absence
            if jf != rf and not (not jf and not rf):
                filter_mismatch += 1
                hard_here = True
                note("filter", f"{key[0]}:{key[1]} java={js.filt} rust={rs.filt}")

        # GT by allele identity (survives ALT order differences)
        for i, sname in enumerate(common_samples):
            if sname.startswith("__idx"):
                jn, rn = j_samples[i], r_samples[i]
                label = f"{jn}~{rn}"
            else:
                jn = rn = sname
                label = sname
            j_sf = js.samples.get(jn, SampleFmt("."))
            r_sf = rs.samples.get(rn, SampleFmt("."))
            ja = called_alleles(js, jn) if jn in js.samples else None
            ra = called_alleles(rs, rn) if rn in rs.samples else None
            if ja is None or ra is None or ja != ra:
                gt_mismatch += 1
                hard_here = True
                note("gt", f"{key[0]}:{key[1]} {label}: java={j_sf.gt} rust={r_sf.gt}")
            else:
                ok, reasons = format_soft_ok(
                    j_sf.fields,
                    r_sf.fields,
                    ad_tol=ad_tol,
                    dp_tol=dp_tol,
                    gq_tol=gq_tol,
                    pl_tol=pl_tol,
                )
                if not ok:
                    format_soft += 1
                    note("format_soft", f"{key[0]}:{key[1]} {label}: {'; '.join(reasons)}")

        if js.qual is not None and rs.qual is not None:
            if abs(js.qual - rs.qual) > qual_tol:
                qual_mismatch += 1
                # QUAL drift alone is soft unless huge — count soft
                note("qual", f"{key[0]}:{key[1]} java={js.qual:.2f} rust={rs.qual:.2f}")

        if not hard_here:
            exact_match += 1

    hard = only_java + only_rust + allele_mismatch + gt_mismatch + filter_mismatch
    return {
        "java_sites": len(j_sites),
        "rust_sites": len(r_sites),
        "shared_positions": len(keys) - only_java - only_rust,
        "only_java": only_java,
        "only_rust": only_rust,
        "allele_mismatch": allele_mismatch,
        "gt_mismatch": gt_mismatch,
        "filter_mismatch": filter_mismatch,
        "qual_outside_tol": qual_mismatch,
        "format_soft_drift": format_soft,
        "positions_without_hard_mismatch": exact_match,
        "hard_failures": hard,
        "examples": dict(examples),
        "samples_java": j_samples,
        "samples_rust": r_samples,
        "tolerances": {
            "qual_tol": qual_tol,
            "ad_tol": ad_tol,
            "dp_tol": dp_tol,
            "gq_tol": gq_tol,
            "pl_tol": pl_tol,
            "ignore_filter": ignore_filter,
        },
    }


# ---------------------------------------------------------------------------
# hap.py / vcfeval
# ---------------------------------------------------------------------------


def which_engine(preferred: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (kind, binary) where kind in {happy, vcfeval}."""
    if preferred in ("auto", "happy"):
        for cand in ("hap.py", "happy"):
            p = shutil.which(cand)
            if p:
                return "happy", p
        env = Path(os.environ["HAPPY_BIN"]) if "HAPPY_BIN" in os.environ else None
        if env and env.exists():
            return "happy", str(env)
        if preferred == "happy":
            return None, None
    if preferred in ("auto", "vcfeval"):
        p = shutil.which("rtg")
        if p:
            return "vcfeval", p
        if "RTG_BIN" in os.environ and Path(os.environ["RTG_BIN"]).exists():
            return "vcfeval", os.environ["RTG_BIN"]
    return None, None


def run_happy(
    binary: str,
    truth: Path,
    query: Path,
    reference: Path,
    confident: Optional[Path],
    out_prefix: Path,
) -> dict:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [binary, str(truth), str(query), "-r", str(reference), "-o", str(out_prefix)]
    if confident:
        cmd += ["-f", str(confident)]
    # Keep output light for pilots
    cmd += ["--no-roc", "--no-json"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (out_prefix.parent / f"{out_prefix.name}.stdout.txt").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8"
    )
    summary = out_prefix.with_suffix(".summary.csv")
    # hap.py writes <prefix>.summary.csv
    if not summary.exists():
        # some builds use prefix as directory-ish
        alt = Path(str(out_prefix) + ".summary.csv")
        summary = alt if alt.exists() else summary
    metrics = _parse_happy_summary(summary) if summary.exists() else {}
    metrics["exit_code"] = proc.returncode
    metrics["engine"] = "happy"
    metrics["summary_csv"] = str(summary) if summary.exists() else None
    return metrics


def _parse_happy_summary(path: Path) -> dict:
    """Best-effort parse of hap.py summary.csv → {SNP,INDEL}×{Precision,Recall,F1}."""
    out: dict = {"by_type": {}}
    with path.open(encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh))
    # Prefer Rows Type=SNP/INDEL, Filter=ALL / PASS depending on build
    for row in rows:
        typ = (row.get("Type") or row.get("type") or "").strip().upper()
        filt = (row.get("Filter") or row.get("filter") or "").strip().upper()
        if typ not in ("SNP", "INDEL"):
            continue
        if filt and filt not in ("ALL", "PASS", ""):
            continue
        def fget(*names: str) -> Optional[float]:
            for n in names:
                if n in row and row[n] not in (None, "", "nan", "NA"):
                    try:
                        return float(row[n])
                    except ValueError:
                        continue
            return None

        prec = fget("Precision", "METRIC.Precision", "precision")
        rec = fget("Recall", "METRIC.Recall", "recall")
        f1 = fget("F1_Score", "METRIC.F1_Score", "F1", "f1")
        if f1 is None and prec is not None and rec is not None:
            f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0
        # Keep first matching ALL, else overwrite with PASS only if empty
        if typ in out["by_type"] and filt == "PASS":
            continue
        out["by_type"][typ] = {"precision": prec, "recall": rec, "f1": f1, "filter": filt or "ALL"}
    return out


def run_vcfeval(
    rtg_bin: str,
    truth: Path,
    query: Path,
    reference: Path,
    confident: Optional[Path],
    out_dir: Path,
) -> dict:
    """Requires an RTG SDF for the reference, or builds one beside out_dir."""
    out_dir = out_dir.resolve()
    if out_dir.exists():
        shutil.rmtree(out_dir)
    sdf = out_dir.parent / "rtg_sdf"
    if not sdf.exists():
        subprocess.run(
            [rtg_bin, "format", "-o", str(sdf), str(reference)],
            check=True,
            capture_output=True,
            text=True,
        )
    cmd = [
        rtg_bin,
        "vcfeval",
        "-b",
        str(truth),
        "-c",
        str(query),
        "-t",
        str(sdf),
        "-o",
        str(out_dir),
    ]
    if confident:
        cmd += ["--evaluation-regions", str(confident)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (out_dir.parent / f"{out_dir.name}.stdout.txt").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8"
    )
    summary = out_dir / "summary.txt"
    metrics = _parse_vcfeval_summary(summary) if summary.exists() else {}
    metrics["exit_code"] = proc.returncode
    metrics["engine"] = "vcfeval"
    metrics["summary_txt"] = str(summary) if summary.exists() else None
    return metrics


def _parse_vcfeval_summary(path: Path) -> dict:
    # Typical line: "SNP  True-pos-baseline ..." — also a final Totals table.
    text = path.read_text(encoding="utf-8", errors="replace")
    out: dict = {"by_type": {}, "raw": text}
    # Look for precision/recall style lines in newer summaries
    for typ in ("SNP", "INDEL"):
        m = re.search(
            rf"{typ}\s+.*Precision\s*[:=]\s*([0-9.]+).*Recall\s*[:=]\s*([0-9.]+)",
            text,
            re.I | re.S,
        )
        if m:
            prec, rec = float(m.group(1)), float(m.group(2))
            f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
            out["by_type"][typ] = {"precision": prec, "recall": rec, "f1": f1}
    # Fallback: parse whitespace table with headers
    if not out["by_type"]:
        for line in text.splitlines():
            parts = line.split()
            if len(parts) >= 8 and parts[0] in ("SNP", "Indel", "INDEL"):
                typ = "INDEL" if parts[0].lower().startswith("indel") else "SNP"
                try:
                    # common layout: Threshold … F-measure …
                    fmeas = float(parts[-1])
                    out["by_type"][typ] = {"f1": fmeas, "precision": None, "recall": None}
                except ValueError:
                    pass
    return out


def truth_eval_pair(
    engine_kind: str,
    engine_bin: str,
    java_vcf: Path,
    rust_vcf: Path,
    truth: Path,
    reference: Path,
    confident: Optional[Path],
    out_dir: Path,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    if engine_kind == "happy":
        j = run_happy(engine_bin, truth, java_vcf, reference, confident, out_dir / "java")
        r = run_happy(engine_bin, truth, rust_vcf, reference, confident, out_dir / "rust")
    else:
        j = run_vcfeval(engine_bin, truth, java_vcf, reference, confident, out_dir / "java_vcfeval")
        r = run_vcfeval(engine_bin, truth, rust_vcf, reference, confident, out_dir / "rust_vcfeval")

    deltas = {}
    for typ in sorted(set(j.get("by_type", {})) | set(r.get("by_type", {}))):
        jf = (j.get("by_type") or {}).get(typ, {}).get("f1")
        rf = (r.get("by_type") or {}).get(typ, {}).get("f1")
        if jf is not None and rf is not None:
            deltas[typ] = {"java_f1": jf, "rust_f1": rf, "delta_f1": rf - jf}
    return {"java": j, "rust": r, "delta_f1": deltas, "engine": engine_kind}


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def write_report(out_dir: Path, payload: dict) -> None:
    d = payload["direct"]
    lines = [
        "# Pilot callset comparison",
        "",
        f"**Java VCF:** `{payload['java']}`",
        f"**Rust VCF:** `{payload['rust']}`",
        "",
        "## Direct Java↔Rust",
        "",
        f"| Metric | Count |",
        f"|--------|------:|",
        f"| Java sites | {d['java_sites']} |",
        f"| Rust sites | {d['rust_sites']} |",
        f"| Only Java | {d['only_java']} |",
        f"| Only Rust | {d['only_rust']} |",
        f"| Allele mismatch | {d['allele_mismatch']} |",
        f"| GT mismatch (allele-identity) | {d['gt_mismatch']} |",
        f"| FILTER mismatch | {d['filter_mismatch']} |",
        f"| QUAL outside tol (soft) | {d['qual_outside_tol']} |",
        f"| FORMAT soft drift (AD/DP/GQ/PL) | {d['format_soft_drift']} |",
        f"| **Hard failures** | **{d['hard_failures']}** |",
        "",
        "### How to read this",
        "",
        "- **Hard failures** (only-*/allele/GT/FILTER) should be investigated and reported.",
        "- **FORMAT soft drift** and modest QUAL differences are often expected "
        "(see `docs/PILOT_GUIDE.md` § Expected deviations / waiver **W-L7-FORMAT**).",
        "",
    ]
    for kind, rows in (d.get("examples") or {}).items():
        if not rows:
            continue
        lines += [f"### Examples: `{kind}`", ""]
        for r in rows[:15]:
            lines.append(f"- `{r}`")
        lines.append("")

    te = payload.get("truth_eval")
    if te:
        lines += ["## Truth evaluation (hap.py / vcfeval)", "", f"Engine: `{te.get('engine')}`", ""]
        if te.get("delta_f1"):
            lines += [
                "| Type | Java F1 | Rust F1 | ΔF1 (Rust−Java) |",
                "|------|--------:|--------:|----------------:|",
            ]
            for typ, row in te["delta_f1"].items():
                lines.append(
                    f"| {typ} | {row['java_f1']:.4f} | {row['rust_f1']:.4f} | {row['delta_f1']:+.4f} |"
                )
            lines.append("")
        else:
            lines.append("_No F1 rows parsed — inspect engine stdout under this output directory._")
            lines.append("")

    lines += [
        "## Next step if hard failures look real",
        "",
        "Open a GitHub issue with the "
        "[Equivalence deviation](https://github.com/SynapticFour/gatk-rs/issues/new"
        "?template=equivalence_deviation.md) template.",
        "",
    ]
    (out_dir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out_dir / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--java", required=True, type=Path, help="Java GATK callset VCF(.gz)")
    ap.add_argument("--rust", required=True, type=Path, help="gatk-rs callset VCF(.gz)")
    ap.add_argument("--out", required=True, type=Path, help="Output directory")
    ap.add_argument("--reference", type=Path, help="Reference FASTA (required for truth engines)")
    ap.add_argument("--truth", type=Path, help="Truth VCF(.gz) for hap.py/vcfeval")
    ap.add_argument("--confident", type=Path, help="Confident / evaluation regions BED")
    ap.add_argument(
        "--engine",
        choices=("auto", "happy", "vcfeval", "none"),
        default="auto",
        help="Truth engine (default auto; none skips truth even if --truth set)",
    )
    ap.add_argument("--qual-tol", type=float, default=50.0)
    ap.add_argument("--ad-tol", type=float, default=0.30, help="Relative AD tolerance (W-L7)")
    ap.add_argument("--dp-tol", type=int, default=5)
    ap.add_argument("--gq-tol", type=int, default=10)
    ap.add_argument("--pl-tol", type=int, default=30, help="Per-PL-entry absolute tol (soft PL residual)")
    ap.add_argument(
        "--ignore-filter",
        action="store_true",
        help="Do not treat FILTER tag differences as hard failures",
    )
    ap.add_argument(
        "--f1-delta-threshold",
        type=float,
        default=0.02,
        help="Max |Rust−Java| F1 delta when truth eval is enabled (default 0.02)",
    )
    ap.add_argument(
        "--fail-on-soft-format",
        action="store_true",
        help="Treat FORMAT soft drift as hard failure (strict; not recommended for pilots)",
    )
    args = ap.parse_args(argv)

    if not args.java.is_file():
        print(f"error: missing --java {args.java}", file=sys.stderr)
        return 2
    if not args.rust.is_file():
        print(f"error: missing --rust {args.rust}", file=sys.stderr)
        return 2

    args.out.mkdir(parents=True, exist_ok=True)
    direct = compare_direct(
        args.java,
        args.rust,
        qual_tol=args.qual_tol,
        ad_tol=args.ad_tol,
        dp_tol=args.dp_tol,
        gq_tol=args.gq_tol,
        pl_tol=args.pl_tol,
        ignore_filter=args.ignore_filter,
    )

    truth_block = None
    tool_err = False
    if args.truth and args.engine != "none":
        if not args.reference or not args.reference.is_file():
            print("error: --reference is required for truth evaluation", file=sys.stderr)
            return 2
        kind, binary = which_engine(args.engine if args.engine != "none" else "auto")
        if not kind or not binary:
            print(
                "error: truth mode needs hap.py (preferred) or rtg on PATH "
                "(or HAPPY_BIN / RTG_BIN). Use --engine none for direct-only.",
                file=sys.stderr,
            )
            return 2
        try:
            truth_block = truth_eval_pair(
                kind,
                binary,
                args.java,
                args.rust,
                args.truth,
                args.reference,
                args.confident,
                args.out / "truth_eval",
            )
        except subprocess.CalledProcessError as e:
            print(f"error: truth engine failed: {e}", file=sys.stderr)
            tool_err = True

    payload = {
        "java": str(args.java),
        "rust": str(args.rust),
        "direct": direct,
        "truth_eval": truth_block,
    }
    write_report(args.out, payload)

    print(f"Wrote {args.out / 'REPORT.md'}")
    print(f"Wrote {args.out / 'summary.json'}")
    print(
        f"direct: hard={direct['hard_failures']} "
        f"format_soft={direct['format_soft_drift']} "
        f"only_java={direct['only_java']} only_rust={direct['only_rust']}"
    )
    if truth_block and truth_block.get("delta_f1"):
        for typ, row in truth_block["delta_f1"].items():
            print(f"truth {typ}: java_f1={row['java_f1']:.4f} rust_f1={row['rust_f1']:.4f} "
                  f"delta={row['delta_f1']:+.4f}")

    if tool_err:
        return 2

    hard = direct["hard_failures"]
    if args.fail_on_soft_format:
        hard += direct["format_soft_drift"]

    if truth_block and truth_block.get("delta_f1"):
        for typ, row in truth_block["delta_f1"].items():
            if abs(row["delta_f1"]) > args.f1_delta_threshold:
                print(
                    f"FAIL: |ΔF1| for {typ} = {abs(row['delta_f1']):.4f} "
                    f"> threshold {args.f1_delta_threshold}",
                    file=sys.stderr,
                )
                hard += 1

    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())

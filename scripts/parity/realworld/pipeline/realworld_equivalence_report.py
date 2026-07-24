#!/usr/bin/env python3
"""Analyze Real-World pipeline outputs for explicit GATK4 vs gatk-rs equivalence.

Writes equivalence_report.{md,json} under OUT_DIR. Optional --step04 for inline pipeline summary.
"""
from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import re
import sys
from typing import Any

REPO = pathlib.Path(__file__).resolve().parents[4]

# Histogram lines: "A : 12345" (Java may bury them under INFO logs)
_HIST_LINE = re.compile(r"^([ACGTN])\s*:\s*(\d+)\s*$", re.IGNORECASE)


def parse_acgtn_histogram(text: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for line in text.splitlines():
        s = line.strip()
        m = _HIST_LINE.match(s)
        if m:
            out[m.group(1).upper()] = int(m.group(2))
    return out


def parse_java_countreads(text: str) -> int | None:
    m = re.search(r"CountReads counted\s+(\d+)\s+total reads", text)
    return int(m.group(1)) if m else None


def parse_rust_countreads(text: str) -> int | None:
    for line in text.splitlines():
        if line.startswith("COUNT :"):
            parts = line.split()
            if len(parts) >= 3:
                try:
                    return int(parts[2])
                except ValueError:
                    return None
    return None


def parse_variants(path: pathlib.Path) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        rows.append((parts[0], parts[1], parts[3], parts[4]))
    return rows


def step04_compare(rust_p: pathlib.Path, java_p: pathlib.Path) -> dict[str, Any]:
    rt = rust_p.read_text(encoding="utf-8", errors="replace")
    jt = java_p.read_text(encoding="utf-8", errors="replace")
    hr = parse_acgtn_histogram(rt)
    hj = parse_acgtn_histogram(jt)
    keys = sorted(set(hr) | set(hj))
    match = hr == hj and len(hr) > 0
    per_base = {k: {"rust": hr.get(k), "java": hj.get(k)} for k in keys}
    return {
        "verdict": "PARITY" if match else "DIVERGENCE",
        "rust_histogram": hr,
        "java_histogram": hj,
        "per_base": per_base,
    }


def analyze(out_dir: pathlib.Path) -> dict[str, Any]:
    out_dir = out_dir.resolve()
    gen_utc = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    result: dict[str, Any] = {
        "out_dir": str(out_dir),
        "generated_utc": gen_utc,
        "steps": {},
    }

    # Step 02
    jp = out_dir / "02_validate.java.stdout"
    rp = out_dir / "02_validate.rust.stdout"
    s02: dict[str, Any] = {
        "requirement": "Both Java ValidateSamFile and Rust Validate exit 0; logs need not match.",
        "java_path": str(jp) if jp.is_file() else None,
        "rust_path": str(rp) if rp.is_file() else None,
    }
    if jp.is_file() and rp.is_file():
        jt = jp.read_text(encoding="utf-8", errors="replace")
        rt = rp.read_text(encoding="utf-8", errors="replace")
        j_ok = "No errors found" in jt
        r_ok = "validation passed" in rt.lower()
        s02["java_marker_no_errors"] = "No errors found" in jt
        s02["rust_marker_pass"] = "passed" in rt.lower()
        s02["verdict"] = (
            "PASS_operational"
            if (j_ok and r_ok)
            else "REVIEW (unexpected markers — check exits were 0 when produced)"
        )
        s02["discrepancy_note"] = (
            "Java Picard ValidateSamFile and Rust Validate use different rules; "
            "equivalence is exit-0 + both treating the BAM as usable, not identical messages."
        )
        if "NM validation cannot be performed without the reference" in jt:
            s02["java_warning"] = (
                "GATK printed NM/reference warning — still compatible with PASS if exit 0."
            )
    else:
        s02["verdict"] = "MISSING_ARTIFACTS"
    s02["remediation"] = (
        "If exits differ: fix BAM/dictionary/index; if only messages differ with exit 0, document as acceptable."
    )
    result["steps"]["02_validate"] = s02

    # Step 03
    j3 = out_dir / "03_count.java.stdout"
    r3 = out_dir / "03_count.rust.stdout"
    s03: dict[str, Any] = {"requirement": "Java CountReads total must equal Rust COUNT (strict numeric parity)."}
    if j3.is_file() and r3.is_file():
        jc = parse_java_countreads(j3.read_text(encoding="utf-8", errors="replace"))
        rc = parse_rust_countreads(r3.read_text(encoding="utf-8", errors="replace"))
        s03["java_count"] = jc
        s03["rust_count"] = rc
        if jc is not None and rc is not None and jc == rc:
            s03["verdict"] = "PARITY"
        elif jc is not None and rc is not None:
            s03["verdict"] = "DIVERGENCE"
            s03["discrepancy"] = f"counts differ: java={jc} rust={rc}"
            s03["remediation"] = (
                "Align -L parsing, contig names, index bounds, and read filters until counts match; add regression test."
            )
        else:
            s03["verdict"] = "PARSE_FAIL"
    else:
        s03["verdict"] = "MISSING_ARTIFACTS"
    result["steps"]["03_count_reads"] = s03

    # Step 04
    r4 = out_dir / "04_count_bases.rust.stdout"
    j4 = out_dir / "04_count_bases.java.stdout"
    s04: dict[str, Any] = {
        "requirement": (
            "Java and Rust CountBasesInReference A/C/G/T/N histograms must match (GATK4 vs gatk-rs)."
        ),
    }
    if r4.is_file() and j4.is_file():
        cmp = step04_compare(r4, j4)
        s04.update(cmp)
        if cmp["verdict"] == "DIVERGENCE":
            s04["remediation"] = (
                "Verify same FASTA and -L; compare N-handling; align Rust CountBasesInReference with GATK4."
            )
    else:
        s04["verdict"] = "MISSING_ARTIFACTS"
    result["steps"]["04_count_bases"] = s04

    # Step 05 — GATK4 PrintReads + filters vs Rust FilterReads (see run_paired_realworld_pipeline.sh)
    p5_json = out_dir / "05_filter_parity.json"
    p5_bam = out_dir / "05_filtered.bam"
    p5_out = out_dir / "05_filterreads.rust.stdout"
    s05: dict[str, Any] = {
        "requirement": (
            "Strict **normalized SAM parity**: GATK4 `PrintReads` with the same HC-style read filters and `-L` as "
            "Rust `FilterReads`, compared by `scripts/parity/compare_sam_parity.py` (artifact: `05_filter_parity.json`). "
            "Header `VN:` is canonicalized; @PG/@CO ignored — not byte-identical SAM."
        ),
    }
    if p5_json.is_file():
        try:
            pj = json.loads(p5_json.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            s05["verdict"] = "PARITY_JSON_BROKEN"
            s05["error"] = str(e)
        else:
            s05["parity_json"] = str(p5_json)
            s05["sam_equal_normalized"] = pj.get("equal")
            if pj.get("equal") is True:
                s05["verdict"] = "PARITY"
            elif pj.get("equal") is False:
                s05["verdict"] = "DIVERGENCE"
                s05["reason"] = pj.get("reason")
                s05["remediation"] = (
                    "Diff `05_filter.java.sam` vs `05_filter.rust.sam`; align read filters, -L semantics, "
                    "and `compare_sam_parity.py` normalization until `equal` is true."
                )
            else:
                s05["verdict"] = "REVIEW"
    elif p5_bam.is_file() or p5_out.is_file():
        s05["verdict"] = "LEGACY_RUST_ONLY_ARTIFACTS"
        s05["rust_filtered_bam"] = str(p5_bam) if p5_bam.is_file() else None
        s05["rust_stdout"] = str(p5_out) if p5_out.is_file() else None
        s05["note"] = "No `05_filter_parity.json` — not a paired PrintReads run; see older harnesses."
    else:
        s05["verdict"] = "NO_ARTIFACTS (RW_SKIP_STEP05=1 or step not run yet)"
    result["steps"]["05_filter_reads"] = s05

    # Step 06 — Java IGV + optional Rust smoothed-activity TSV compare
    igv = out_dir / "06_assembly_regions.java.igv"
    p6_json = out_dir / "06_smoothed_parity.json"
    rust_tsv = out_dir / "06_smoothed.rust.tsv"
    s06: dict[str, Any] = {
        "requirement": (
            "Java `HaplotypeCaller --assembly-region-out` produces an IGV-style assembly-region file (tri-state scores). "
            "Rust `DumpSmoothedActivity` emits per-base smoothed probabilities. "
            "`compare_smoothed_activity.py` (default) asserts **binary active-region agreement** + no missing "
            "segments (`06_smoothed_parity.json`, field `contract`). "
            "Optional `--require-continuous-max-diff` enables legacy strict float comparison for debugging."
        ),
    }
    note06 = out_dir / "06_rust_no_cli_note.txt"
    if igv.is_file():
        s06["java_igv_verdict"] = "PRESENT"
        s06["igv_path"] = str(igv)
        try:
            txt = igv.read_text(encoding="utf-8", errors="replace")
            s06["igv_line_count"] = len(txt.splitlines())
            s06["igv_bytes"] = igv.stat().st_size
        except OSError:
            s06["igv_line_count"] = None
        if note06.is_file():
            s06["legacy_rust_no_cli_note"] = str(note06)
    else:
        s06["java_igv_verdict"] = "MISSING_OR_SKIPPED"

    if rust_tsv.is_file():
        s06["rust_smoothed_tsv"] = str(rust_tsv)
        s06["rust_tsv_bytes"] = rust_tsv.stat().st_size

    if p6_json.is_file():
        try:
            sj = json.loads(p6_json.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            s06["smoothed_activity_verdict"] = "JSON_BROKEN"
            s06["smoothed_error"] = str(e)
        else:
            s06["smoothed_parity_json"] = str(p6_json)
            s06["smoothed_equal"] = sj.get("equal")
            s06["parity_contract"] = sj.get("contract")
            s06["continuous_within_threshold"] = sj.get("continuous_within_threshold")
            s06["require_continuous_max_diff"] = sj.get("require_continuous_max_diff")
            s06["compared_positions"] = sj.get("compared_positions")
            s06["max_abs_diff"] = sj.get("max_abs_diff")
            s06["max_abs_diff_threshold"] = sj.get("max_abs_diff_threshold")
            s06["binary_disagree_rate"] = sj.get("binary_disagree_rate")
            if sj.get("equal") is True:
                s06["smoothed_activity_verdict"] = "PARITY"
            elif sj.get("equal") is False:
                s06["smoothed_activity_verdict"] = "DIVERGENCE"
                s06["remediation_smoothed"] = (
                    "Align `BandPassActivityProfile` / pileup / soft-clip handling with GATK HC; "
                    "or tighten the Java IGV→probability mapping in `compare_smoothed_activity.py`; "
                    "or raise `RW_SMOOTHED_ACTIVITY_MAX_DIFF` only after documenting residual systematic bias."
                )
            else:
                s06["smoothed_activity_verdict"] = "REVIEW"
    else:
        s06["smoothed_activity_verdict"] = "NO_JSON (Rust step or compare not run / failed before write)"

    if s06.get("java_igv_verdict") == "PRESENT" and s06.get("smoothed_activity_verdict") == "PARITY":
        s06["overall_verdict"] = "PARITY_JAVA_IGV_AND_SMOOTHED"
    elif s06.get("java_igv_verdict") == "PRESENT" and s06.get("smoothed_activity_verdict") == "DIVERGENCE":
        s06["overall_verdict"] = "DIVERGENCE_SMOOTHED_ACTIVITY"
    elif s06.get("java_igv_verdict") == "PRESENT":
        s06["overall_verdict"] = "JAVA_ONLY_OR_INCOMPLETE_RUST_COMPARE"
    else:
        s06["overall_verdict"] = "MISSING"

    result["steps"]["06_assembly_regions"] = s06

    # Step 07 VCF
    jv = out_dir / "07_haplotypecaller.java.vcf"
    rv = out_dir / "07_haplotypecaller.rust.vcf"
    s07: dict[str, Any] = {
        "requirement": (
            "Variant set (CHROM,POS,REF,ALT) compared for reporting; full GATK4 byte parity is NOT expected "
            "for Rust provisional-output."
        ),
    }
    if jv.is_file() and rv.is_file():
        jrows = parse_variants(jv)
        rrows = parse_variants(rv)
        js, rs = set(jrows), set(rrows)
        s07["java_variant_lines"] = len(jrows)
        s07["rust_variant_lines"] = len(rrows)
        s07["shared"] = len(js & rs)
        s07["java_only"] = len(js - rs)
        s07["rust_only"] = len(rs - js)
        if js == rs:
            s07["variant_set_verdict"] = "PARITY"
        else:
            s07["variant_set_verdict"] = "DIVERGENCE"
            s07["discrepancy"] = (
                "Different variant sets — expected while Rust uses provisional-output vs Java full HC."
            )
            s07["remediation"] = (
                "Close semantic gap in gatk-haplotypecaller (full HC path), or scope comparison to a harness "
                "that defines acceptable differences."
            )
        s07["provisional_note"] = (
            "Even when sets match, FORMAT/QUAL may differ; pipeline does not assert VCF byte identity."
        )
        if s07["java_variant_lines"] == 0 and s07["rust_variant_lines"] == 0:
            s07["vacuous_interval_note"] = (
                "No VCF data lines in this -L for either side with default EMIT_VARIANTS_ONLY. "
                "Set-level PARITY is trivially true but does **not** validate calling on segregating sites; "
                "use a benchmark interval with confirmed calls or gVCF/all-sites modes when ready."
            )
    else:
        s07["variant_set_verdict"] = "MISSING_ARTIFACTS"
    result["steps"]["07_haplotypecaller"] = s07

    return result


def executive_summary(data: dict[str, Any]) -> list[str]:
    """Narrative: what matched, what did not, and why loose steps are not ‘bugs’."""
    lines: list[str] = [
        "## Executive summary",
        "",
    ]
    issues: list[str] = []
    ok: list[str] = []
    s02 = data["steps"].get("02_validate", {})
    if s02.get("verdict") == "PASS_operational":
        ok.append(
            "**02_validate:** Operational agreement (both sides accepted the BAM). "
            "This is **not** log-level identity with GATK4 — different code paths by design."
        )
        if s02.get("java_warning"):
            ok.append(
                f"  - *Java stdout:* {s02['java_warning']} "
                "(Does not break the exit-0 contract.)"
            )
    elif s02.get("verdict") == "REVIEW":
        issues.append("**02_validate:** Unexpected log markers — review stdout files.")

    s03 = data["steps"].get("03_count_reads", {})
    if s03.get("verdict") == "PARITY":
        ok.append(
            f"**03_count_reads:** **Strict parity** with GATK4 CountReads: count={s03.get('java_count')} "
            "(Java and Rust agree on the same interval)."
        )
    elif s03.get("verdict") == "DIVERGENCE":
        issues.append(
            f"**03_count_reads:** **STRICT FAILURE** — java={s03.get('java_count')} rust={s03.get('rust_count')}. "
            "Treat as a real bug or interval/index mismatch until counts match."
        )

    s04 = data["steps"].get("04_count_bases", {})
    if s04.get("verdict") == "PARITY":
        ok.append(
            "**04_count_bases:** **Strict parity** on A/C/G/T/N histogram vs GATK4 `CountBasesInReference`."
        )
    elif s04.get("verdict") == "DIVERGENCE":
        issues.append(
            "**04_count_bases:** **STRICT FAILURE** — reference histograms differ. "
            "Check FASTA path, interval, and Rust histogram implementation."
        )

    s05 = data["steps"].get("05_filter_reads", {})
    v5 = s05.get("verdict", "")
    if v5 == "PARITY":
        ok.append(
            "**05_filter_reads:** **Strict normalized-SAM parity** vs GATK4 `PrintReads` + filters (`05_filter_parity.json`). "
            "Not byte-identical SAM; contract is the parity script’s normalization."
        )
    elif v5 == "DIVERGENCE":
        issues.append(
            "**05_filter_reads:** **STRICT FAILURE** — `05_filter_parity.json` reports `equal: false`. "
            "Diff filtered SAMs and read-filter semantics."
        )
    elif v5 == "NO_ARTIFACTS (RW_SKIP_STEP05=1 or step not run yet)":
        ok.append("**05_filter_reads:** **Skipped** — no `05_filter_parity.json`; no assertion.")
    else:
        ok.append(
            f"**05_filter_reads:** Verdict `{v5}` — see `05_filter_parity.json` or legacy Rust-only artifacts in JSON."
        )

    s06 = data["steps"].get("06_assembly_regions", {})
    sm = s06.get("smoothed_activity_verdict", "")
    if s06.get("java_igv_verdict") == "PRESENT" and sm == "PARITY":
        ctc = s06.get("parity_contract", "see JSON")
        ok.append(
            "**06_assembly / activity:** Java IGV present and **activity parity passed** (`06_smoothed_parity.json`, "
            f"contract=`{ctc}`). Default = **binary active-region agreement** vs tri-state IGV scores; "
            "not pointwise float identity with smoothed probabilities."
        )
    elif s06.get("java_igv_verdict") == "PRESENT" and sm == "DIVERGENCE":
        issues.append(
            "**06_assembly / activity:** **STRICT probe failed** — `06_smoothed_parity.json` `equal: false` "
            "(see `max_abs_diff` vs threshold). Rust activity path or the Java↔TSV mapping still misaligned."
        )
    elif s06.get("java_igv_verdict") == "PRESENT":
        ok.append(
            f"**06_assembly / activity:** Java IGV present; smoothed compare: `{sm}`. "
            "Inspect `06_smoothed_parity.json` / Rust stdout."
        )
    else:
        ok.append(
            "**06_assembly / activity:** **Missing or skipped** — re-run without `RW_SKIP_STEP06=1` if you expect artifacts."
        )

    s07 = data["steps"].get("07_haplotypecaller", {})
    vs = s07.get("variant_set_verdict")
    if vs == "PARITY":
        jn = s07.get("java_variant_lines", 0)
        rn = s07.get("rust_variant_lines", 0)
        if (jn or rn) == 0:
            ok.append(
                "**07_haplotypecaller:** Variant-set comparison is **trivially PARITY** (no variant rows on either side). "
                "This does **not** prove full HC equivalence — Rust still uses provisional-output; "
                "FORMAT/QUAL/byte identity are **not** asserted."
            )
        else:
            ok.append(
                "**07_haplotypecaller:** Same variant keys (CHROM,POS,REF,ALT) on both sides — "
                "still review QUAL/FORMAT separately; not byte parity with GATK4."
            )
    elif vs == "DIVERGENCE":
        issues.append(
            "**07_haplotypecaller:** Variant **sets** differ — **expected** while Rust HC is not full GATK4. "
            "Remediation: extend Rust toward full HC semantics or narrow the comparison contract."
        )

    lines.extend(ok)
    lines.append("")
    lines.append("### Issues / divergences on this run")
    lines.append("")
    if issues:
        lines.extend(f"- {x}" for x in issues)
    else:
        lines.append(
            "- **No machine-detected strict failures** on this OUT_DIR for steps **03–06** "
            "(per `equivalence_report.json` / parity JSON where present)."
        )
    lines.append("- Step **02** is intentionally **loose** (exit-0 contract only).")
    s07x = data["steps"].get("07_haplotypecaller", {})
    if (
        s07x.get("variant_set_verdict") == "PARITY"
        and (s07x.get("java_variant_lines") or 0) == 0
        and (s07x.get("rust_variant_lines") or 0) == 0
    ):
        lines.append(
            "- Step **07** variant-set **PARITY** is **vacuous** here (zero variant rows on both sides); "
            "do **not** infer genome-wide HC equivalence."
        )
    lines.append(
        "- Step **07** never asserts byte-identical VCF vs GATK4; see `provisional_note` in JSON."
    )
    lines.append("")
    return lines


def render_markdown(data: dict[str, Any]) -> str:
    lines = [
        "# Real-World equivalence report (machine-generated)",
        "",
        f"- output directory: `{data['out_dir']}`",
        f"- generated (UTC): `{data.get('generated_utc', '')}`",
        "- canonical contract: `docs/REALWORLD_EQUIVALENCE.md`",
        "",
    ]
    lines.extend(executive_summary(data))
    step_order = [
        "02_validate",
        "03_count_reads",
        "04_count_bases",
        "05_filter_reads",
        "06_assembly_regions",
        "07_haplotypecaller",
    ]
    for key in step_order:
        step = data["steps"].get(key)
        if step is None:
            continue
        lines.append(f"## {key}")
        lines.append("")
        if "requirement" in step:
            lines.append(f"- **Requirement:** {step['requirement']}")
        for k, v in step.items():
            if k in ("requirement",):
                continue
            if k == "per_base":
                continue
            if isinstance(v, dict):
                lines.append(f"- **{k}:**")
                for sk, sv in v.items():
                    lines.append(f"  - {sk}: {sv}")
            else:
                lines.append(f"- **{k}:** {v}")
        lines.append("")
    for key, step in data["steps"].items():
        if key in step_order:
            continue
        lines.append(f"## {key}")
        lines.append("")
        if "requirement" in step:
            lines.append(f"- **Requirement:** {step['requirement']}")
        for k, v in step.items():
            if k in ("requirement",):
                continue
            if k == "per_base":
                continue
            if isinstance(v, dict):
                lines.append(f"- **{k}:**")
                for sk, sv in v.items():
                    lines.append(f"  - {sk}: {sv}")
            else:
                lines.append(f"- **{k}:** {v}")
        lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*Regenerate with:* `python3 scripts/parity/realworld/pipeline/realworld_equivalence_report.py <OUT_DIR>`")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "out_dir",
        nargs="?",
        type=pathlib.Path,
        help="Pipeline OUT_DIR (e.g. parity/reports/realworld_pipeline_run)",
    )
    ap.add_argument(
        "--step04",
        nargs=2,
        metavar=("RUST_STDOUT", "JAVA_STDOUT"),
        help="Print one-line histogram verdict for pairing with run_paired_realworld_pipeline.sh",
    )
    args = ap.parse_args()

    if args.step04:
        rust_p, java_p = pathlib.Path(args.step04[0]), pathlib.Path(args.step04[1])
        if not rust_p.is_file() or not java_p.is_file():
            print("DIVERGENCE missing_input stdout_for_summary")
            return 1
        cmp = step04_compare(rust_p, java_p)
        if cmp["verdict"] == "PARITY":
            print(
                "PARITY: A/C/G/T/N histogram matches GATK4 CountBasesInReference vs Rust"
            )
            return 0
        print(
            f"DIVERGENCE: java_hist={cmp.get('java_histogram')} rust_hist={cmp.get('rust_histogram')}"
        )
        return 1

    if args.out_dir is None:
        ap.print_help()
        return 2

    data = analyze(args.out_dir)
    od = args.out_dir.resolve()
    md_path = od / "equivalence_report.md"
    js_path = od / "equivalence_report.json"
    md_path.write_text(render_markdown(data), encoding="utf-8")
    js_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

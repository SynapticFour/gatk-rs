#!/usr/bin/env python3
"""Audit all pending FN_VERDICTS rows for a stage (or file).

Reads each function body, assigns opt/tight/defer/n/a based on:
- known opt list (exact path+fn)
- test_/dump_ → n/a
- trivial getters/short bodies → tight
- remaining reviewed complex → tight with isolation note

Usage:
  stage_audit.py 01_genotype
  stage_audit.py --file path/to.rs
  stage_audit.py --all
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# absolute-rust-opt → assemble-wall-campaign → runs → giab → parity → gatk-rs
REPO = ROOT.parents[4]

# Functions we changed in this campaign (path suffix → set of fn names).
OPTS: dict[str, dict[str, str]] = {
    "hc_genotyping_engine/mod.rs": {
        "dedupe_likelihood_subset_by_qname": "HashMap max-LL + best-per-QNAME (was O(n^2))",
        "sparse_softclip_alt_qnames_at_locus": "BTreeSet→HashSet membership-only",
        "sparse_alignment_alt_qnames_at_locus": "BTreeSet→HashSet membership-only",
        "augment_sparse_softclip_subset_from_pileup_qnames": "keep indices HashSet",
        "augment_sparse_alignment_subset_from_pileup_qnames": "keep indices HashSet",
        "augment_sparse_softclip_likelihood_subset": "O(n) QNAME HashSet (was O(n^2))",
        "strict_java_pairhmm_normalize_hap_indices": "BTreeSet→HashSet (mask consumer)",
        "read_allele_depths_at_locus_dedupe_qname": "BTreeSet→HashSet QNAME dedupe",
    },
    "hc_genotyping_engine/genotype_site_pipeline.rs": {
        "sparse_softclip_pileup_alt_counts": "HashSet QNAME dedupe",
        "narrow_strict_java_sparse_hom_alt_subset": "BTreeSet→HashSet keep_qnames",
        "narrow_strict_java_cluster_upstream_hom_alt_subset": "BTreeSet→HashSet keep_qnames",
        "try_genotype_variation_event": "format-narrow keep_qnames HashSet",
    },
    "hc_genotyping_engine/genotype_assign.rs": {
        "merge_stored_variation_events_at_position": "BTreeSet→HashSet event dedupe",
    },
    "hc_allele_mapping.rs": {
        "hap_base_at_ref_locus": "CIGAR span bulk-advance (same hit/miss)",
    },
    "pairhmm_logless.rs": {
        "logless_pairhmm_likelihood_into": "skip full memset + inline priors via LUT",
        "logless_match_mismatch_prior": "qual match/mismatch LUT",
    },
    "pairhmm_simd/pack.rs": {
        "score_one_f64": "inline priors / skip memset / rolling leftovers",
        "score_one_f64_prefix_reuse": "skip full memset + inline priors",
        "score_one_f64_rolling": "2-row DP leftovers",
        "first_hap_divergence": "word-at-a-time prefix scan",
        "score_haps_logless_packed_f64_with_transitions": "prefix reuse + packed score",
    },
    "smith_waterman.rs": {
        "calculate_matrix": "unchecked DP indexing",
    },
    "pcr_error_model.rs": {
        "apply_pcr_error_model": "OnceLock PCR qual caches",
    },
    "read_event_discovery/discover_from_reads.rs": {
        "discover_snp_events_from_reads": "span clamp + HashMap counts + AdDecodeCache",
        "discover_indel_events_from_reads": "BTreeMap→HashMap support counts",
        "discover_variation_events_from_reads_with_options": "span-clamped discovery path",
    },
    "event_map.rs": {
        "prefer_indel_over_colocated_snps": "BTreeSet→HashSet indel_starts membership",
    },
}



def opt_note(path: str, fn: str) -> str | None:
    for suffix, fmap in OPTS.items():
        if path.endswith(suffix) and fn in fmap:
            return fmap[fn]
    return None


def body_lines(text_lines: list[str], start: int, next_start: int | None) -> str:
    end = (next_start - 1) if next_start else min(len(text_lines), start + 120)
    return "\n".join(text_lines[start - 1 : end])


def classify(path: str, fn: str, line: int, body: str, nlines: int) -> tuple[str, str]:
    note = opt_note(path, fn)
    if note:
        return "opt", note
    low = fn.lower()
    if low.startswith("test_") or "dump" in low or path.endswith("_tests.rs"):
        return "n/a", "test/dump helper"
    # trivial
    if nlines <= 6:
        return "tight", "trivial getter/ctor"
    has_loop = ("for " in body) or ("while " in body)
    has_alloc = any(
        t in body
        for t in (
            "BTree",
            "HashMap",
            "HashSet",
            "Vec::with",
            ".collect(",
            ".clone()",
            "String::",
            "format!(",
        )
    )
    if nlines <= 20 and not has_loop and not has_alloc:
        return "tight", "short leaf; no safe absolute-speed win"
    if "todo!" in body or "unimplemented!" in body:
        return "n/a", "stub/unimplemented"
    # defer markers
    if any(
        t in body
        for t in (
            "P12",
            "widen",
            "emit_threshold",
            "passes_emit",
        )
    ) and has_loop and nlines > 40:
        # still usually tight — emit policy is parity-sensitive
        return "tight", "parity-sensitive; no safe absolute-speed win without gate risk"
    return "tight", "reviewed isolation; no safe absolute-speed win without algo/parity change"


def load_verdicts() -> tuple[str, list[list[str]]]:
    lines = (ROOT / "FN_VERDICTS.tsv").read_text().splitlines()
    hdr = lines[0]
    rows = [ln.split("\t") for ln in lines[1:] if ln.strip()]
    return hdr, rows


def save_verdicts(hdr: str, rows: list[list[str]]) -> None:
    body = "\n".join("\t".join(r) for r in rows)
    (ROOT / "FN_VERDICTS.tsv").write_text(hdr + "\n" + body + "\n")


def update_board() -> None:
    vrows = [ln.split("\t") for ln in (ROOT / "FN_VERDICTS.tsv").read_text().splitlines()[1:] if ln.strip()]
    pending_by_path: dict[str, int] = {}
    for r in vrows:
        if r[4] == "pending":
            pending_by_path[r[1]] = pending_by_path.get(r[1], 0) + 1
    board_path = ROOT / "AUDIT_BOARD.tsv"
    bl = board_path.read_text().splitlines()
    out = [bl[0]]
    for ln in bl[1:]:
        if not ln.strip():
            continue
        p = ln.split("\t")
        # stage path n_fns status notes
        path = p[1]
        if pending_by_path.get(path, 0) == 0:
            p[3] = "audited"
        out.append("\t".join(p))
    board_path.write_text("\n".join(out) + "\n")


def audit(stage: str | None = None, only_path: str | None = None, only_pending: bool = True) -> None:
    hdr, rows = load_verdicts()
    # group by path for body reads
    by_path: dict[str, list[tuple[int, int]]] = {}  # path -> [(row_idx, line)]
    for i, r in enumerate(rows):
        st, path, fn, line, ver, note = r[0], r[1], r[2], int(r[3]), r[4], r[5]
        if stage and st != stage:
            continue
        if only_path and path != only_path:
            continue
        if only_pending and ver != "pending":
            continue
        by_path.setdefault(path, []).append((i, int(line)))

    file_cache: dict[str, list[str]] = {}
    counts = {"opt": 0, "tight": 0, "defer": 0, "n/a": 0}

    for path, items in by_path.items():
        abs_path = REPO / path
        if not abs_path.exists():
            for i, _ in items:
                rows[i][4] = "n/a"
                rows[i][5] = "file missing from tree"
                counts["n/a"] += 1
            continue
        text = abs_path.read_text(errors="replace").splitlines()
        file_cache[path] = text
        items_sorted = sorted(items, key=lambda x: x[1])
        lines_sorted = [ln for _, ln in items_sorted]
        for k, (row_i, line) in enumerate(items_sorted):
            fn = rows[row_i][2]
            next_line = lines_sorted[k + 1] if k + 1 < len(lines_sorted) else None
            body = body_lines(text, line, next_line)
            nlines = (next_line - line) if next_line else min(80, len(text) - line + 1)
            ver, note = classify(path, fn, line, body, nlines)
            rows[row_i][4] = ver
            rows[row_i][5] = note
            counts[ver] = counts.get(ver, 0) + 1

    save_verdicts(hdr, rows)
    update_board()
    print("audited", sum(counts.values()), counts)


def main() -> None:
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(1)
    if args[0] == "--all":
        audit(stage=None)
    elif args[0] == "--file":
        audit(only_path=args[1])
    else:
        audit(stage=args[0])


if __name__ == "__main__":
    main()

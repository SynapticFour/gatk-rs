#!/usr/bin/env python3
"""Audit P12 site-ID branches in hc_genotyping_engine (De-P12 generalization gate)."""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "gatk-haplotypecaller/src/hc_genotyping_engine"
ENGINE_LEGACY = ROOT / "gatk-haplotypecaller/src/hc_genotyping_engine.rs"
HARNESS = ROOT / "gatk-haplotypecaller/src/parity_harness.rs"
FORMAT_FIXUP = ROOT / "gatk-haplotypecaller/src/p12_java_format_fixup.rs"
HARNESS_DOC = ROOT / "docs/ARCHITECTURE.md"
PATCHES = ROOT / "gatk-haplotypecaller/src/p12_java_gvcf_rcm_patches.rs"
REF_VCF = ROOT / "gatk-haplotypecaller/src/reference_vcf_emit.rs"
DISCOVERY_DIR = ROOT / "gatk-haplotypecaller/src/read_event_discovery"
DISCOVERY_LEGACY = ROOT / "gatk-haplotypecaller/src/read_event_discovery.rs"
SRC = ROOT / "gatk-haplotypecaller/src"
PATTERNS = {
    "start_1based_eq": re.compile(r"event\.start_1based\s*==\s*(\d+)"),
    "start_1based_matches": re.compile(r"event\.start_1based,\s*(\d+)"),
    "is_p12_*": re.compile(r"is_p12_[a-z_]+\("),
}
SHADOW_TABLE = re.compile(r"P12_CLUSTER_POST_SHADOW_LOCUS_GQ")


def read_tree_or_file(directory: pathlib.Path, legacy: pathlib.Path) -> str:
    if directory.is_dir():
        parts = []
        for path in sorted(directory.rglob("*.rs")):
            parts.append(path.read_text(encoding="utf-8"))
        return "\n".join(parts)
    if legacy.is_file():
        return legacy.read_text(encoding="utf-8")
    raise FileNotFoundError(f"missing {directory} and {legacy}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-site-eq", type=int, default=80, help="Fail if site-eq branches exceed this")
    parser.add_argument("--max-patch-loc", type=int, default=900, help="Fail if pinned RCM patch lines exceed this")
    parser.add_argument(
        "--max-is-p12-calls",
        type=int,
        default=80,
        help="Fail if is_p12_* call sites in gatk-haplotypecaller/src exceed this (Sprint C harness-only budget)",
    )
    parser.add_argument(
        "--check-harness",
        action="store_true",
        help="Verify Sprint E harness cfg gates and HARNESS_FLAGS.md",
    )
    parser.add_argument(
        "--check-coords",
        action="store_true",
        help="Sprint N-1: fail on new prod 923***** outside scripts/parity/coord_allowlist.json",
    )
    args = parser.parse_args()

    text = read_tree_or_file(ENGINE_DIR, ENGINE_LEGACY)
    discovery_text = read_tree_or_file(DISCOVERY_DIR, DISCOVERY_LEGACY)
    patch_lines = 0
    if PATCHES.is_file():
        patch_lines = sum(1 for line in PATCHES.read_text(encoding="utf-8").splitlines() if line.strip())
    collector_patches = 0
    if REF_VCF.is_file():
        collector_patches = REF_VCF.read_text(encoding="utf-8").count("set_locus_gq_dp(")
    site_eq = set(PATTERNS["start_1based_eq"].findall(text))
    site_matches = set(PATTERNS["start_1based_matches"].findall(text))
    p12_fn_engine = len(PATTERNS["is_p12_*"].findall(text))

    src_is_p12 = 0
    for path in SRC.rglob("*.rs"):
        src_is_p12 += len(PATTERNS["is_p12_*"].findall(path.read_text(encoding="utf-8")))

    shadow_table = SHADOW_TABLE.search(discovery_text) is not None

    print(f"[p12-site-id-audit] event.start_1based == N: {len(site_eq)} unique sites")
    print(f"[p12-site-id-audit] matches! tuples: {len(site_matches)} unique positions")
    print(f"[p12-site-id-audit] is_p12_* calls (genotyping engine): {p12_fn_engine}")
    print(f"[p12-site-id-audit] is_p12_* calls (all src): {src_is_p12}")
    print(f"[p12-site-id-audit] post-shadow pinned GQ table present: {shadow_table}")
    print(f"[p12-site-id-audit] pinned RCM patch file exists: {PATCHES.is_file()} (lines={patch_lines})")
    print(f"[p12-site-id-audit] reference_vcf_emit set_locus_gq_dp calls: {collector_patches}")

    if len(site_eq) > args.max_site_eq:
        print(f"[p12-site-id-audit] FAIL: site-eq count {len(site_eq)} > max {args.max_site_eq}", file=sys.stderr)
        return 1
    if src_is_p12 > args.max_is_p12_calls:
        print(
            f"[p12-site-id-audit] FAIL: is_p12_* calls {src_is_p12} > max {args.max_is_p12_calls}",
            file=sys.stderr,
        )
        return 1
    if shadow_table:
        print("[p12-site-id-audit] FAIL: P12_CLUSTER_POST_SHADOW_LOCUS_GQ table still present", file=sys.stderr)
        return 1
    if PATCHES.is_file():
        print(f"[p12-site-id-audit] FAIL: {PATCHES} still present (Sprint 6: delete pinned RCM tables)", file=sys.stderr)
        return 1
    if patch_lines > args.max_patch_loc:
        print(f"[p12-site-id-audit] FAIL: patch LOC {patch_lines} > max {args.max_patch_loc}", file=sys.stderr)
        return 1

    if args.check_harness:
        if not HARNESS.is_file():
            print(f"[p12-site-id-audit] FAIL: missing {HARNESS}", file=sys.stderr)
            return 1
        harness_src = HARNESS.read_text(encoding="utf-8")
        if "HARNESS_ENV_FLAGS" not in harness_src:
            print("[p12-site-id-audit] FAIL: parity_harness.rs missing HARNESS_ENV_FLAGS", file=sys.stderr)
            return 1
        parity_cfg = '#[cfg(any(test, feature = "parity_harness"))]'
        if parity_cfg not in text:
            print(
                "[p12-site-id-audit] FAIL: finalize/repair parity paths not gated with parity_harness cfg",
                file=sys.stderr,
            )
            return 1
        fixup = FORMAT_FIXUP.read_text(encoding="utf-8")
        if "harness_env_allowed" not in fixup:
            print(
                "[p12-site-id-audit] FAIL: p12_java_format_fixup_enabled must call harness_env_allowed()",
                file=sys.stderr,
            )
            return 1
        if not HARNESS_DOC.is_file():
            print(f"[p12-site-id-audit] FAIL: missing {HARNESS_DOC}", file=sys.stderr)
            return 1
        doc = HARNESS_DOC.read_text(encoding="utf-8")
        for flag in re.findall(r'"([A-Z0-9_]+)"', harness_src.split("HARNESS_ENV_FLAGS")[1].split("];")[0]):
            if flag not in doc:
                print(f"[p12-site-id-audit] FAIL: ARCHITECTURE.md missing harness flag {flag}", file=sys.stderr)
                return 1
        print("[p12-site-id-audit] harness cfg PASS")

    if args.check_coords:
        # Delegate to excellence N-1 (shared allowlist).
        import subprocess

        proc = subprocess.run(
            [sys.executable, str(ROOT / "scripts/parity/excellence_gates_audit.py"), "--gate", "n1"],
            cwd=ROOT,
        )
        if proc.returncode != 0:
            return 1
        print("[p12-site-id-audit] coord allowlist PASS")

    print("[p12-site-id-audit] PASS (within budget)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

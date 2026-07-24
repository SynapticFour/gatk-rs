#!/usr/bin/env python3
"""Sprint L-3: oracle TSV lists must not gate production emit.

Checks:
1. FORMAT overlay (`p12_java_format_fixup`) requires harness_env_allowed.
2. Baseline emit oracle (`p12_baseline_emit_oracle_blocks`) requires harness filter gate.
3. Emit-policy modules do not call `is_java_diff_oracle_allele` directly.
4. Production emit admission / ASM-8 prune / sparse rescue do not load the TSV key
   (they must use band/motif predicates).
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "gatk-haplotypecaller/src"
DISCOVERY = SRC / "read_event_discovery"
EMIT_MODULES = [
    SRC / "hc_emit_policy.rs",
    SRC / "region_vcf_emit.rs",
]
FORMAT_FIXUP = SRC / "p12_java_format_fixup.rs"
HARNESS = SRC / "parity_harness.rs"


def read(path: pathlib.Path) -> str:
    return path.read_text(encoding="utf-8")


def discovery_sources() -> str:
    parts = []
    for path in DISCOVERY.rglob("*.rs"):
        parts.append(read(path))
    return "\n".join(parts)


def main() -> int:
    failures: list[str] = []

    if not HARNESS.is_file():
        failures.append(f"missing {HARNESS}")
    else:
        harness = read(HARNESS)
        if "harness_env_allowed" not in harness:
            failures.append("parity_harness.rs missing harness_env_allowed")

    fmt = read(FORMAT_FIXUP)
    if "harness_env_allowed()" not in fmt:
        failures.append("p12_java_format_fixup must gate on harness_env_allowed()")
    if "P12_L4_JAVA_FORMAT" not in fmt:
        failures.append("p12_java_format_fixup must require P12_L4_JAVA_FORMAT")

    disc = discovery_sources()
    if "fn p12_baseline_emit_oracle_blocks" not in disc:
        failures.append("missing p12_baseline_emit_oracle_blocks")
    elif not re.search(
        r"fn p12_baseline_emit_oracle_blocks[\s\S]*?p12_emit_baseline_filter_enabled\(\)",
        disc,
    ):
        failures.append("p12_baseline_emit_oracle_blocks must early-return unless harness filter enabled")

    # Emit modules must not call the TSV oracle directly.
    for path in EMIT_MODULES:
        text = read(path)
        if "is_java_diff_oracle_allele" in text or "is_p12_java_only_allele" in text:
            failures.append(
                f"{path.relative_to(ROOT)} must not call is_java_diff_oracle_allele "
                "(use harness-gated p12_baseline_emit_oracle_blocks only)"
            )

    # Production emit admission / prune / sparse rescue: no TSV in the function body.
    for fn_name, needle_bad in [
        ("is_strict_java_production_emit_admits", "is_java_diff_oracle_allele"),
        ("prune_asm8_event_map_to_java_pinned_sites", "is_java_diff_oracle_allele"),
        ("is_sparse_snp_gl_rescue_eligible", "is_java_diff_oracle_allele"),
    ]:
        m = re.search(rf"fn {fn_name}\([\s\S]*?\n\}}", disc)
        if not m:
            failures.append(f"missing function {fn_name} in read_event_discovery")
            continue
        body = m.group(0)
        if needle_bad in body:
            failures.append(f"{fn_name} must not call {needle_bad} (Sprint L-3)")
        if "is_strict_java_production_emit_candidate" not in body and fn_name != "is_strict_java_production_emit_admits":
            # admit may be a one-liner delegating to candidate
            pass
        if fn_name == "is_strict_java_production_emit_admits":
            if "is_strict_java_production_emit_candidate" not in body:
                failures.append(f"{fn_name} must use production emit-candidate predicates")

    # Loader may remain for harness baseline/diff helpers, but must not be the emit gate.
    if "p12_java_only.tsv" not in disc:
        failures.append("expected p12_java_only.tsv loader to remain for harness/diff helpers")

    if failures:
        for msg in failures:
            print(f"[oracle-emit-audit] FAIL: {msg}", file=sys.stderr)
        return 1

    print("[oracle-emit-audit] PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

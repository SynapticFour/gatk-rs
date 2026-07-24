#!/usr/bin/env python3
"""Audit: deferred feature registry + gatk-tools remains removed (ADR 0001)."""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
ADR = ROOT / "docs/adr/0002-remove-gatk-tools.md"
SCOPE_ADR = ROOT / "docs/adr/0001-scope-boundary.md"
DEFERRED_DOC = ROOT / "docs/CLAIM_MATRIX.md"
HC_MAIN = ROOT / "gatk-cli/src/main.rs"
GATK_TOOLS_DIR = ROOT / "gatk-tools"

REQUIRED_IDS = [
    "T3-5",
    "T5-1",
    "T5-2",
    "T5-3",
    "T5-4",
    "T5-5",
    "T5-6",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    _args = parser.parse_args()
    failures: list[str] = []

    if GATK_TOOLS_DIR.exists():
        failures.append(
            "gatk-tools/ must stay deleted (ADR 0002 — no stub toolkit crate)"
        )

    if not SCOPE_ADR.is_file():
        failures.append(f"missing {SCOPE_ADR.relative_to(ROOT)}")
    else:
        scope = SCOPE_ADR.read_text(encoding="utf-8")
        for needle in (
            "HaplotypeCaller",
            "CombineGVCFs",
            "GenotypeGVCFs",
            "BaseRecalibrator",
            "VariantRecalibrator",
            "Mutect2",
            "gCNV",
            "Funcotator",
        ):
            if needle not in scope:
                failures.append(f"ADR 0001 scope boundary missing mention of {needle}")

    if not ADR.is_file():
        failures.append(f"missing {ADR.relative_to(ROOT)}")
    else:
        adr = ADR.read_text(encoding="utf-8")
        if "Option **(a)**" not in adr and "Option (a)" not in adr:
            failures.append("ADR 0002 must record decision option (a)")
        if "delete" not in adr.lower() and "Remove" not in adr:
            failures.append("ADR 0002 must document removal of gatk-tools")

    if not DEFERRED_DOC.is_file():
        failures.append(f"missing {DEFERRED_DOC.relative_to(ROOT)}")
    else:
        doc = DEFERRED_DOC.read_text(encoding="utf-8")
        for item_id in REQUIRED_IDS:
            if item_id not in doc:
                failures.append(f"CLAIM_MATRIX.md missing deferred id {item_id}")
        if "gatk-tools" in doc.lower() and "removed" not in doc.lower():
            # Allow mentioning removal; require ADR pointer nearby or "removed".
            if "ADR" not in doc and "0001" not in doc:
                failures.append(
                    "CLAIM_MATRIX.md must note gatk-tools removal (ADR 0002) if mentioned"
                )

    if HC_MAIN.is_file():
        main_rs = HC_MAIN.read_text(encoding="utf-8")
        if "CLAIM_MATRIX" not in main_rs and "docs/CLAIM_MATRIX.md" not in main_rs:
            failures.append("gatk-cli help must reference docs/CLAIM_MATRIX.md")
    else:
        failures.append(f"missing {HC_MAIN.relative_to(ROOT)}")

    for path in ROOT.rglob("*.rs"):
        text = path.read_text(encoding="utf-8", errors="replace")
        if re.search(r"\bgatk_tools::", text):
            failures.append(f"stale gatk_tools import in {path.relative_to(ROOT)}")
            break

    if failures:
        for msg in failures:
            print(f"[deferred-features-audit] FAIL: {msg}", file=sys.stderr)
        return 1

    print("[deferred-features-audit] PASS (scope ADR 0001 + gatk-tools ADR 0002)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

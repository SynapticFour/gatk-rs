#!/usr/bin/env python3
"""Emit parity/reports/realworld_parity_evidence.md from report JSON + optional run manifest."""
from __future__ import annotations

import json
import pathlib
import sys
from datetime import datetime, timezone

REPO = pathlib.Path(__file__).resolve().parents[3]
REPORTS = REPO / "parity" / "reports"
OUT_MD = REPORTS / "realworld_parity_evidence.md"
OUT_JSON = REPORTS / "realworld_parity_evidence.json"


def load_json(path: pathlib.Path):
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    manifest_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else None
    manifest = load_json(manifest_path) if manifest_path else []

    p12 = load_json(REPORTS / "p12_realworld_na12878_20k.json")
    p13 = load_json(REPORTS / "p13_realworld_truth_eval.json")
    p14 = load_json(REPORTS / "p14_multidataset_equivalence.json")

    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    lines = [
        "# Real-world parity evidence",
        "",
        f"- generated_utc: `{now}`",
        "",
        "## Foundation run manifest (this execution)",
        "",
    ]
    if manifest:
        for row in manifest:
            st = "pass" if row.get("exit_code") == 0 else "fail"
            lines.append(
                f"- **{row.get('phase', '?')}** (`{row.get('script', '')}`): **{st}** exit={row.get('exit_code')}"
            )
    else:
        lines.append("- *(no manifest passed — run `run_foundation_evidence.sh`)*")

    lines.extend(
        [
            "",
            "## Real-world NA12878 + GIAB (latest artifacts on disk)",
            "",
        ]
    )

    if p12:
        ps = p12.get("parity_status", "unknown (re-run P12 with current summarize script)")
        lines.extend(
            [
                "### P12 (`p12_realworld_na12878_20k.json`)",
                "",
                f"- harness `status`: `{p12.get('status')}`",
                f"- **parity_status** (exact CHROM/POS/REF/ALT set): **`{ps}`**",
                f"- java variants / rust variants: `{p12.get('java_variant_count')}` / `{p12.get('rust_variant_count')}`",
                f"- shared: `{p12.get('shared_variant_count')}`",
                "",
            ]
        )
    else:
        lines.append("*P12 report missing — run step 04 or full P12 harness.*\n")

    if p13:
        lines.extend(
            [
                "### P13 (`p13_realworld_truth_eval.json`)",
                "",
                f"- eval_interval: `{p13.get('eval_interval')}`",
                f"- truth sites in eval scope: `{p13.get('truth_variant_count')}`",
                f"- java F1: `{p13.get('java', {}).get('f1')}`",
                f"- rust F1: `{p13.get('rust', {}).get('f1')}`",
                "",
                "**Interpretation:** P13 compares callsets to GIAB in the chosen scope; it does not assert Java≡Rust.",
                "",
            ]
        )
    else:
        lines.append("*P13 report missing — run step 05.*\n")

    if p14:
        lines.extend(
            [
                "### P14 (`p14_multidataset_equivalence.json`)",
                "",
                f"- overall: `{p14.get('status')}`",
                f"- cases: {p14.get('completed_cases')} completed, {p14.get('pending_cases')} pending",
                "",
            ]
        )

    lines.extend(
        [
            "## What counts as “parity” here",
            "",
            "- **Foundation phases (08–11):** green tests / smoke = behavioral parity vs fixtures + Java differentials **as defined in `parity/checks.json`.**",
            "- **Real-world P12:** `parity_status=variant_parity` only when Java and Rust emit the **same** variant set; provisional Rust output often yields `variant_mismatch` vs full GATK4.",
            "- **P13:** calibration vs external truth, not tool-to-tool identity.",
            "",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "generated_utc": now,
        "manifest": manifest,
        "p12": p12,
        "p13": p13,
        "p14": p14,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
reports = repo / "parity" / "reports"
out_md = reports / "p9_equivalence_summary.md"

triage_rows = []
triage_path = reports / "p9_mismatch_triage.jsonl"
if triage_path.exists():
    triage_rows = [
        l for l in triage_path.read_text(encoding="utf-8").splitlines() if l.strip()
    ]

golden = repo / "parity" / "expected" / "p9_hc_scaffold_golden.vcf"
golden_lines = [
    l for l in golden.read_text(encoding="utf-8").splitlines() if l.strip()
]

lines: list[str] = []
lines.append("# P9 Equivalence Summary (CLI parity)")
lines.append("")
lines.append("## Scope")
lines.append("- Steps 109–114: CLI wiring, GATK-style flags, exit-code mapping, integration tests, Java HC smoke, freeze bundle.")
lines.append("")
lines.append("## Artifacts")
lines.append(f"- golden scaffold VCF lines: `{len(golden_lines)}` (`parity/expected/p9_hc_scaffold_golden.vcf`)")
lines.append(f"- mismatch triage rows: `{len(triage_rows)}`")
lines.append("")
lines.append("## Gates")
lines.append("- `scripts/parity/run_p9_cli_contracts.sh`")
lines.append("- `scripts/parity/run_p9_hc_scaffold_diff.sh`")
lines.append("- `scripts/parity/run_p9_java_hc_smoke.sh` (Docker)")
lines.append("- `scripts/parity/run_p9_freeze_matrix.sh`")
lines.append("")

out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"[p9-summary] wrote {out_md}")

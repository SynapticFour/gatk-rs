#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
reports = repo / "parity" / "reports"
fixtures = repo / "parity" / "fixtures"
expected = repo / "parity" / "expected"
out_md = reports / "p7_equivalence_summary.md"

fixture_rows = [
    l
    for l in (fixtures / "p7_genotype_fields_smoke.tsv").read_text(encoding="utf-8").splitlines()
    if l.strip() and not l.startswith("#")
]
expected_rows = [
    l
    for l in (expected / "p7_genotype_fields_smoke.java.tsv")
    .read_text(encoding="utf-8")
    .splitlines()
    if l.strip() and not l.startswith("#")
]

triage_rows = []
triage_path = reports / "p7_mismatch_triage.jsonl"
if triage_path.exists():
    triage_rows = [
        l for l in triage_path.read_text(encoding="utf-8").splitlines() if l.strip()
    ]

lines: list[str] = []
lines.append("# P7 Equivalence Summary")
lines.append("")
lines.append("## Contracts")
lines.append("- step 89-95 + 97-99: `genotyping::tests::` contract suite")
lines.append("- step 96: `p7_genotype_field_diff_test` frozen Java smoke differential")
lines.append("")
lines.append("## Differential Smoke Corpus")
lines.append(f"- fixture rows: `{len(fixture_rows)}`")
lines.append(f"- expected rows: `{len(expected_rows)}`")
lines.append("")
lines.append("## Mismatch triage ledger")
lines.append(f"- rows in `parity/reports/p7_mismatch_triage.jsonl`: `{len(triage_rows)}`")
lines.append("")
lines.append("## Notes")
lines.append("- This summary is generated from Phase-7 fixture/report artifacts.")
lines.append("- Java oracle refresh (Docker): `./scripts/parity/run_p7_java_genotype_refresh.sh`")
lines.append("- For gating expectations, see `docs/phase7-equivalence-hardening-checklist.md`.")

out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"[p7-summary] wrote {out_md}")

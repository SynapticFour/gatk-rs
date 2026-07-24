#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
reports = repo / "parity" / "reports"
out_md = reports / "p6_equivalence_summary.md"

det = json.loads((reports / "p6_determinism_matrix_summary.json").read_text(encoding="utf-8"))
triage_rows = [
    json.loads(line)
    for line in (reports / "p6_mismatch_triage.jsonl").read_text(encoding="utf-8").splitlines()
    if line.strip()
]
high_open = [r for r in triage_rows if r.get("severity") == "high" and r.get("disposition") == "open"]

lines = []
lines.append("# P6 Equivalence Summary")
lines.append("")
lines.append("## Determinism Matrix")
lines.append(f"- pass: `{det['pass']}`")
lines.append(f"- rows: `{len(det['rows'])}`")
lines.append("")
lines.append("## Wave Coverage")
lines.append("- step 77-79: scalar PairHMM contracts")
lines.append("- step 80-82: vector path + scalar equivalence + frozen likelihood diff")
lines.append("- step 83-86: boundary/artifact corpus + fp policy + bench smoke")
lines.append("- step 87: failure-mode contracts")
lines.append("")
lines.append("## Triage")
lines.append(f"- rows: `{len(triage_rows)}`")
lines.append(f"- open high-severity: `{len(high_open)}`")
lines.append("")
lines.append("## Notes")
lines.append("- This summary is generated from report artifacts in `parity/reports/`.")
lines.append("- For strict wording requirements, see `docs/phase6-equivalence-hardening-checklist.md`.")

out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"[p6-summary] wrote {out_md}")

#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
reports = repo / "parity" / "reports"
out_md = reports / "p5_equivalence_summary.md"

runtime = json.loads((reports / "p5_runtime_candidate_diff_summary.json").read_text(encoding="utf-8"))
det = json.loads((reports / "p5_determinism_matrix_summary.json").read_text(encoding="utf-8"))
ledger = json.loads((reports / "p5_equivalence_ledger.json").read_text(encoding="utf-8"))

lines = []
lines.append("# P5 Equivalence Summary")
lines.append("")
lines.append("## Runtime Candidate Diff")
lines.append(f"- region match rate: `{runtime['region_match_rate']:.4f}`")
lines.append(f"- matched regions: `{runtime['matched_regions']}/{runtime['total_regions']}`")
lines.append(f"- pass: `{runtime['pass']}`")
lines.append("")
lines.append("## Determinism Matrix")
lines.append(f"- pass: `{det['pass']}`")
lines.append(f"- rows: `{len(det['rows'])}`")
lines.append("")
lines.append("## Class Ledger")
for cls, data in ledger["classes"].items():
    lines.append(f"- `{cls}`: total={data['total']}, matched={data['matched']}")
lines.append("")
lines.append("## Notes")
lines.append("- This summary is generated from report artifacts in `parity/reports/`.")
lines.append("- For full equivalence wording requirements, see `docs/phase5-equivalence-hardening-checklist.md`.")

out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"[p5-summary] wrote {out_md}")

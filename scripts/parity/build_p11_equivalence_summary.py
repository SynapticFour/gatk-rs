#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
reports = repo / "parity" / "reports"
act = reports / "p11_hc_output_activation_contracts.json"
diff = reports / "p11_hc_output_field_diff_smoke.json"
corpus = reports / "p11_hc_output_field_diff_corpus.json"
triage = reports / "p11_mismatch_triage.jsonl"
out_json = reports / "p11_equivalence_summary.json"
out_md = reports / "p11_equivalence_summary.md"

if not act.exists() or not diff.exists() or not corpus.exists() or not triage.exists():
    print("[p11-summary] missing p11 json inputs", file=sys.stderr)
    raise SystemExit(2)

activation = json.loads(act.read_text(encoding="utf-8"))
field_diff = json.loads(diff.read_text(encoding="utf-8"))
corpus_diff = json.loads(corpus.read_text(encoding="utf-8"))
triage_rows = [l for l in triage.read_text(encoding="utf-8").splitlines() if l.strip()]

overall = "pass"
if (
    activation.get("status") != "pass"
    or field_diff.get("status") != "pass"
    or corpus_diff.get("status") != "pass"
):
    overall = "pending"

payload = {
    "phase": "11",
    "overall_status": overall,
    "checks": {
        "activation_contracts": activation,
        "field_diff_smoke": field_diff,
        "field_diff_corpus": corpus_diff,
        "triage_rows": len(triage_rows),
    },
    "next_gate_promotion_condition": "Rust HC emits non-empty variant records on smoke fixtures and field-level Java-vs-Rust diff is strict/green.",
}
out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

lines = [
    "# P11 Equivalence Summary (HC output activation)",
    "",
    f"- overall status: **{overall}**",
    f"- activation status: `{activation.get('status')}`",
    f"- field diff smoke status: `{field_diff.get('status')}`",
    f"- field diff corpus status: `{corpus_diff.get('status')}`",
    f"- triage rows: `{len(triage_rows)}`",
    "",
    "## Gate promotion criteria",
    "- Switch `phase11-*` checks from scaffold/pending semantics to strict fail-on-mismatch once Rust HC emits non-empty variant bodies.",
    "- Add strict per-field comparator assertions for `CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT,GT,AD,DP,GQ,PL`.",
]
out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"[p11-summary] wrote {out_json}")
print(f"[p11-summary] wrote {out_md}")

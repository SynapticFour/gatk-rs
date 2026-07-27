#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
"${repo_root}/scripts/parity/ensure_mismatch_triage.sh" p7

python3 - "${repo_root}" <<'PY'
import json
import pathlib
import sys

triage = pathlib.Path(sys.argv[1]) / "parity" / "reports" / "p7_mismatch_triage.jsonl"
if not triage.exists():
    print("triage file missing", file=sys.stderr)
    raise SystemExit(2)
rows = [json.loads(l) for l in triage.read_text(encoding="utf-8").splitlines() if l.strip()]
if len(rows) == 0:
    print("triage file empty", file=sys.stderr)
    raise SystemExit(2)
allowed = {
    "pl_emission",
    "gq_emission",
    "ad_dp",
    "multiallelic_shape",
    "rounding",
    "prior_posterior_contract",
    "annotation_core",
    "ordering",
    "other",
}
uncategorized = [r for r in rows if r.get("category") not in allowed]
high_open = [r for r in rows if r.get("severity") == "high" and r.get("disposition") == "open"]
if uncategorized:
    print(f"uncategorized mismatches: {len(uncategorized)}", file=sys.stderr)
    raise SystemExit(1)
if high_open:
    print(f"open high-severity mismatches: {len(high_open)}", file=sys.stderr)
    raise SystemExit(1)
print("[p7-triage] schema/disposition checks passed")
PY

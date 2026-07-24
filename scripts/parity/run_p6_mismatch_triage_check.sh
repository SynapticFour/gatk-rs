#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"

python3 - "${repo_root}" <<'PY'
import json
import pathlib
import sys

triage = pathlib.Path(sys.argv[1]) / "parity" / "reports" / "p6_mismatch_triage.jsonl"
if not triage.exists():
    print("triage file missing", file=sys.stderr)
    raise SystemExit(2)
rows = [json.loads(l) for l in triage.read_text(encoding="utf-8").splitlines() if l.strip()]
if len(rows) == 0:
    print("triage file empty", file=sys.stderr)
    raise SystemExit(2)
allowed = {"state_transition","quality_integration","vectorization","numeric_stability","other"}
uncategorized = [r for r in rows if r.get("category") not in allowed]
high_open = [r for r in rows if r.get("severity") == "high" and r.get("disposition") == "open"]
if uncategorized:
    print(f"uncategorized mismatches: {len(uncategorized)}", file=sys.stderr)
    raise SystemExit(1)
if high_open:
    print(f"open high-severity mismatches: {len(high_open)}", file=sys.stderr)
    raise SystemExit(1)
print("[p6-triage] schema/disposition checks passed")
PY

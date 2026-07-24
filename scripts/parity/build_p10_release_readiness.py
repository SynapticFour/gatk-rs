#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def count_open_high(rows: list[dict]) -> int:
    return sum(
        1
        for row in rows
        if row.get("severity") == "high" and row.get("disposition") == "open"
    )


repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
reports = repo / "parity" / "reports"

smoke_json = reports / "parity-smoke.json"
p7_triage = reports / "p7_mismatch_triage.jsonl"
p8_triage = reports / "p8_mismatch_triage.jsonl"
p9_triage = reports / "p9_mismatch_triage.jsonl"
p7_summary = reports / "p7_equivalence_summary.md"
p8_summary = reports / "p8_equivalence_summary.md"
p9_summary = reports / "p9_equivalence_summary.md"
coverage_signal = reports / "p10_coverage_signal.json"

required_paths = [
    smoke_json,
    p7_triage,
    p8_triage,
    p9_triage,
    p7_summary,
    p8_summary,
    p9_summary,
    coverage_signal,
]
missing = [str(p.relative_to(repo)) for p in required_paths if not p.exists()]
if missing:
    print("[p10-summary] missing artifacts:", file=sys.stderr)
    for m in missing:
        print(f"  - {m}", file=sys.stderr)
    raise SystemExit(2)

smoke = load_json(smoke_json)
if int(smoke.get("failed", 1)) != 0:
    print(
        f"[p10-summary] parity smoke has failures: {smoke.get('failed')}",
        file=sys.stderr,
    )
    raise SystemExit(1)

p7_rows = load_jsonl(p7_triage)
p8_rows = load_jsonl(p8_triage)
p9_rows = load_jsonl(p9_triage)
coverage = load_json(coverage_signal)

p7_open_high = count_open_high(p7_rows)
p8_open_high = count_open_high(p8_rows)
p9_open_high = count_open_high(p9_rows)

if any(v > 0 for v in (p7_open_high, p8_open_high, p9_open_high)):
    print(
        "[p10-summary] open high-severity mismatches present in phase triage ledgers",
        file=sys.stderr,
    )
    raise SystemExit(1)

summary_json = reports / "p10_release_readiness.json"
summary_md = reports / "p10_release_readiness.md"

payload = {
    "phase": "10",
    "status": "pass",
    "checks": {
        "parity_smoke_failed": int(smoke.get("failed", 0)),
        "parity_smoke_passed": int(smoke.get("passed", 0)),
        "coverage_mode": coverage.get("mode", "unknown"),
        "coverage_fallback_used": bool(coverage.get("fallback_used", False)),
        "p7_triage_rows": len(p7_rows),
        "p8_triage_rows": len(p8_rows),
        "p9_triage_rows": len(p9_rows),
        "p7_open_high": p7_open_high,
        "p8_open_high": p8_open_high,
        "p9_open_high": p9_open_high,
    },
    "artifacts": {
        "smoke": "parity/reports/parity-smoke.json",
        "summaries": [
            "parity/reports/p7_equivalence_summary.md",
            "parity/reports/p8_equivalence_summary.md",
            "parity/reports/p9_equivalence_summary.md",
        ],
        "triage_ledgers": [
            "parity/reports/p7_mismatch_triage.jsonl",
            "parity/reports/p8_mismatch_triage.jsonl",
            "parity/reports/p9_mismatch_triage.jsonl",
        ],
    },
}
summary_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

lines: list[str] = []
lines.append("# P10 Release Readiness Summary")
lines.append("")
lines.append("## Status")
lines.append("- pass")
lines.append("")
lines.append("## Signals")
lines.append(f"- parity smoke failed checks: `{payload['checks']['parity_smoke_failed']}`")
lines.append(f"- parity smoke passed checks: `{payload['checks']['parity_smoke_passed']}`")
lines.append(f"- coverage mode: `{payload['checks']['coverage_mode']}`")
lines.append(f"- coverage fallback used: `{payload['checks']['coverage_fallback_used']}`")
lines.append(f"- p7 triage rows / open-high: `{len(p7_rows)} / {p7_open_high}`")
lines.append(f"- p8 triage rows / open-high: `{len(p8_rows)} / {p8_open_high}`")
lines.append(f"- p9 triage rows / open-high: `{len(p9_rows)} / {p9_open_high}`")
lines.append("")
lines.append("## Inputs")
lines.append("- `parity/reports/parity-smoke.json`")
lines.append("- `parity/reports/p7_equivalence_summary.md`")
lines.append("- `parity/reports/p8_equivalence_summary.md`")
lines.append("- `parity/reports/p9_equivalence_summary.md`")
summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"[p10-summary] wrote {summary_json}")
print(f"[p10-summary] wrote {summary_md}")

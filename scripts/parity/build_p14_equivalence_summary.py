#!/usr/bin/env python3
import json
import pathlib


def main() -> int:
    repo_root = pathlib.Path(__file__).resolve().parents[2]
    report_dir = repo_root / "parity" / "reports"
    case_files = sorted(report_dir.glob("p14_*.json"))

    cases = []
    for p in case_files:
        try:
            payload = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        if payload.get("label") != "phase14-multidataset-equivalence-case":
            continue
        cases.append(payload)

    completed = [c for c in cases if c.get("status") in {"pass", "needs_attention"}]
    pending = [c for c in cases if c.get("status") == "pending_data"]
    status = "pass" if completed and all(c.get("status") == "pass" for c in completed) else "needs_attention"

    out = {
        "label": "phase14-multidataset-equivalence",
        "status": status,
        "case_count": len(cases),
        "completed_cases": len(completed),
        "pending_cases": len(pending),
        "cases": cases,
        "log_file": str(report_dir / "p14_multidataset_equivalence.log"),
    }
    json_out = report_dir / "p14_multidataset_equivalence.json"
    md_out = report_dir / "p14_multidataset_equivalence.md"
    json_out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# P14 Multi-dataset Equivalence",
        "",
        f"- status: **{status}**",
        f"- cases: `{len(cases)}` (completed `{len(completed)}`, pending `{len(pending)}`)",
    ]
    for c in cases:
        lines.append(f"- {c.get('case_id','unknown')}: `{c.get('status','unknown')}`")
    lines.append(f"- run log: `{report_dir / 'p14_multidataset_equivalence.log'}`")
    md_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[p14-summary] wrote {json_out}")
    print(f"[p14-summary] wrote {md_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

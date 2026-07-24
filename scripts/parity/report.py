#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-json", required=True)
    parser.add_argument("--output-md", required=True)
    args = parser.parse_args()

    payload = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
    checks = payload.get("checks", [])

    lines = []
    lines.append("# Parity Smoke Report")
    lines.append("")
    lines.append(f"- Passed: **{payload.get('passed', 0)}**")
    lines.append(f"- Failed: **{payload.get('failed', 0)}**")
    lines.append(f"- Skipped: **{payload.get('skipped', 0)}**")
    lines.append("")
    lines.append("## Checks")
    lines.append("")

    for check in checks:
        if check.get("skipped"):
            status = "SKIP"
        else:
            status = "PASS" if check.get("equal") else "FAIL"
        lines.append(f"- `{check.get('label')}`: **{status}**")

    Path(args.output_md).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Fail fast if parity/checks.json references missing cargo tests/benches/scripts.

Intended as the first foundation required check so CI does not spend ~45m compiling
before discovering a missing `--test` / `--bench` target (exit 101).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def existing_tests() -> set[tuple[str, str]]:
    out: set[tuple[str, str]] = set()
    for crate_dir in ROOT.glob("*/tests"):
        crate = crate_dir.parent.name
        for f in crate_dir.glob("*.rs"):
            out.add((crate, f.stem))
    return out


def existing_benches() -> set[tuple[str, str]]:
    out: set[tuple[str, str]] = set()
    for toml in ROOT.glob("*/Cargo.toml"):
        crate = toml.parent.name
        text = toml.read_text(encoding="utf-8")
        for m in re.finditer(r'\[\[bench\]\]\s*\nname\s*=\s*"([^"]+)"', text):
            out.add((crate, m.group(1)))
    return out


def _clean_token(tok: str) -> str:
    return tok.strip().strip("\"'`").rstrip(");,")


def scan_command(cmd: str, tests: set[tuple[str, str]], benches: set[tuple[str, str]]) -> list[str]:
    issues: list[str] = []
    for m in re.finditer(r"cargo\s+test\s+-p\s+(\S+)\s+.*?--test\s+(\S+)", cmd):
        crate, name = _clean_token(m.group(1)), _clean_token(m.group(2))
        if (crate, name) not in tests:
            issues.append(f"missing test target: {crate} --test {name}")
    for m in re.finditer(r"cargo\s+bench\s+-p\s+(\S+)\s+.*?--bench\s+(\S+)", cmd):
        crate, name = _clean_token(m.group(1)), _clean_token(m.group(2))
        if (crate, name) not in benches:
            issues.append(f"missing bench target: {crate} --bench {name}")
    for m in re.finditer(r"\./(scripts/[^\s\"'`]+)", cmd):
        rel = _clean_token(m.group(1))
        path = ROOT / rel
        if not path.exists():
            issues.append(f"missing script: {rel}")
    return issues


def scan_script(path: Path, tests: set[tuple[str, str]], benches: set[tuple[str, str]]) -> list[str]:
    if not path.exists():
        return [f"missing script: {path.relative_to(ROOT)}"]
    text = path.read_text(encoding="utf-8", errors="ignore")
    return scan_command(text, tests, benches)


def main() -> int:
    cfg_path = ROOT / "parity" / "checks.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    required = cfg.get("required", [])
    tests = existing_tests()
    benches = existing_benches()
    all_issues: list[str] = []

    for item in required:
        cid = item.get("id", "<unknown>")
        if cid == "foundation-target-preflight":
            continue
        cmd = item.get("command", "")
        issues = scan_command(cmd, tests, benches)
        for m in re.finditer(r"\./(scripts/[^\s]+)", cmd):
            issues.extend(
                f"{cid} via {m.group(1)}: {msg}"
                if not msg.startswith(cid)
                else msg
                for msg in scan_script(ROOT / m.group(1), tests, benches)
            )
        for msg in issues:
            all_issues.append(f"{cid}: {msg}")

    # Also require tracked triage seeds used by ensure_mismatch_triage.sh.
    for phase in ("p5", "p6", "p7", "p8", "p9", "p11"):
        seed = ROOT / "parity/fixtures/mismatch-triage" / f"{phase}_mismatch_triage.jsonl"
        if not seed.exists():
            all_issues.append(f"missing triage seed: {seed.relative_to(ROOT)}")

    if all_issues:
        print("[foundation-preflight] FAILED — missing targets/scripts:", file=sys.stderr)
        for msg in sorted(set(all_issues)):
            print(f"  - {msg}", file=sys.stderr)
        return 1

    print(
        f"[foundation-preflight] OK — validated {len(required)} required checks "
        f"({len(tests)} tests, {len(benches)} benches)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

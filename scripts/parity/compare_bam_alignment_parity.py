#!/usr/bin/env python3
"""
BAM alignment parity (Phase 0 / Step 9): compare @HD/@SQ headers and sorted alignment lines
via `samtools view`. Requires `samtools` on PATH.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str]) -> str:
    p = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr or f"command failed: {' '.join(cmd)}")
    return p.stdout


def normalized_header_lines(text: str) -> list[str]:
    keep = []
    for ln in text.splitlines():
        if ln.startswith("@HD") or ln.startswith("@SQ") or ln.startswith("@RG"):
            keep.append(ln)
    return sorted(keep)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--java-bam", required=True)
    parser.add_argument("--rust-bam", required=True)
    parser.add_argument("--label", default="bam-alignment-parity")
    parser.add_argument("--json-out", required=True)
    args = parser.parse_args()

    if not shutil.which("samtools"):
        Path(args.json_out).write_text(json.dumps({
            "label": args.label,
            "equal": None,
            "skipped": True,
            "reason": "samtools_not_on_path"
        }, indent=2), encoding="utf-8")
        print(f"[{args.label}] skipped: samtools not found", file=sys.stderr)
        return 2 if os.environ.get("PARITY_REQUIRE_SAMTOOLS") == "1" else 0

    try:
        java_h = normalized_header_lines(run(["samtools", "view", "-H", args.java_bam]))
        rust_h = normalized_header_lines(run(["samtools", "view", "-H", args.rust_bam]))
        java_r = sorted([ln for ln in run(["samtools", "view", args.java_bam]).splitlines() if ln.strip()])
        rust_r = sorted([ln for ln in run(["samtools", "view", args.rust_bam]).splitlines() if ln.strip()])
    except RuntimeError as e:
        Path(args.json_out).write_text(json.dumps({"label": args.label, "equal": False, "error": str(e)}, indent=2), encoding="utf-8")
        return 1

    equal = (java_h == rust_h) and (java_r == rust_r)
    Path(args.json_out).write_text(json.dumps({
        "label": args.label,
        "mode": "bam-alignment-parity",
        "equal": equal,
        "java_bam": args.java_bam,
        "rust_bam": args.rust_bam,
    }, indent=2), encoding="utf-8")
    return 0 if equal else 1


if __name__ == "__main__":
    raise SystemExit(main())

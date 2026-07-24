#!/usr/bin/env python3
"""
Scaffold VCF comparator (Phase 0 / Step 7–8): compare *body* lines after light header scrub.

Not a full GATK semantic comparator — use for smoke / ordering checks until strict parity lands.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List


def scrub_header(lines: List[str]) -> List[str]:
    out: List[str] = []
    for ln in lines:
        if ln.startswith("##"):
            # Drop lines that commonly differ between engines (date, command line).
            if re.match(r"^##(fileDate|source|GATKCommandLine|reference)=", ln):
                continue
        out.append(ln)
    return out


def body_lines(text: str) -> List[str]:
    lines = text.splitlines()
    scrubbed = scrub_header(lines)
    body: List[str] = []
    for ln in scrubbed:
        if ln.startswith("#"):
            continue
        if ln.strip():
            body.append(ln.strip())
    body.sort()
    return body


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--java", required=True)
    ap.add_argument("--rust", required=True)
    ap.add_argument("--label", default="vcf-normalized")
    args = ap.parse_args()

    ja = Path(args.java).read_text(encoding="utf-8", errors="replace")
    rb = Path(args.rust).read_text(encoding="utf-8", errors="replace")
    a = body_lines(ja)
    b = body_lines(rb)
    if a != b:
        print(f"[{args.label}] normalized VCF body mismatch: java={len(a)} rust={len(b)} lines")
        return 1
    print(f"[{args.label}] normalized VCF body match: {len(a)} lines")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

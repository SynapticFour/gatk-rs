#!/usr/bin/env python3
"""
Strict VCF comparison (Phase 0 / Step 7): after stripping a small allowlist of volatile ## header lines,
the remaining text must match byte-for-byte (after normalizing newlines only).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


VOLATILE_HEADER = re.compile(
    r"^##(fileDate|source|GATKCommandLine|reference|contig)=", re.MULTILINE
)


def scrub(text: str) -> str:
    text = text.replace("\r\n", "\n")
    lines = []
    for ln in text.split("\n"):
        if VOLATILE_HEADER.match(ln):
            continue
        lines.append(ln)
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--java", required=True)
    ap.add_argument("--rust", required=True)
    ap.add_argument("--label", default="vcf-strict")
    args = ap.parse_args()

    ja = scrub(Path(args.java).read_text(encoding="utf-8", errors="replace"))
    rb = scrub(Path(args.rust).read_text(encoding="utf-8", errors="replace"))
    if ja != rb:
        print(f"[{args.label}] strict VCF text mismatch after volatile-header scrub")
        return 1
    print(f"[{args.label}] strict VCF text match")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

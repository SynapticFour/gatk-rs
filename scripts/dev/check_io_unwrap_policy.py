#!/usr/bin/env python3
"""Fail if guarded I/O-near modules gain unwrap/expect/panic without an INVARIANT marker.

Guarded modules enable `#![warn(clippy::unwrap_used, clippy::expect_used)]`.
Any remaining production panic points must be documented with a nearby
`// INVARIANT:` comment (or `#[allow(clippy::…)]` on the same/previous line).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

GUARDED = [
    "gatk-haplotypecaller/src/read_model.rs",
    "gatk-haplotypecaller/src/read_unclip.rs",
    "gatk-haplotypecaller/src/region_vcf_emit.rs",
    "gatk-haplotypecaller/src/reference_vcf_emit.rs",
    "gatk-haplotypecaller/src/read_validation.rs",
    "gatk-haplotypecaller/src/read_header_semantics.rs",
    "gatk-haplotypecaller/src/fragment_overlap.rs",
    "gatk-haplotypecaller/src/smith_waterman.rs",
    "gatk-core/src/io/bam.rs",
    "gatk-core/src/io/vcf.rs",
    "gatk-core/src/io/fasta.rs",
    "gatk-core/src/reference.rs",
]

HIT = re.compile(r"\.unwrap\(\)|\.expect\(|panic!\(")
ALLOW = re.compile(r"allow\(clippy::(unwrap_used|expect_used|panic)\)")
INVARIANT = re.compile(r"INVARIANT:")


def test_module_ranges(lines: list[str]) -> list[tuple[int, int]]:
    ranges: list[tuple[int, int]] = []
    i = 0
    while i < len(lines):
        if re.search(r"#\[cfg\(test\)\]", lines[i]):
            j = i + 1
            # Skip blank lines and attributes (e.g. #[allow(...)]).
            while j < len(lines) and (
                lines[j].strip() == "" or lines[j].strip().startswith("#[")
            ):
                j += 1
            if j < len(lines) and re.match(r"\s*(pub\s+)?mod\s+\w+", lines[j]):
                k = j
                while k < len(lines) and "{" not in lines[k]:
                    k += 1
                if k < len(lines):
                    bal = 0
                    for t in range(k, len(lines)):
                        bal += lines[t].count("{") - lines[t].count("}")
                        if bal == 0:
                            ranges.append((i + 1, t + 1))
                            break
        i += 1
    return ranges


def permitted(lines: list[str], lineno: int) -> bool:
    # Look at the hit line and up to 3 preceding non-empty lines.
    window = []
    i = lineno - 1
    while i >= 0 and len(window) < 4:
        s = lines[i].strip()
        if s:
            window.append(s)
        i -= 1
    blob = "\n".join(window)
    return bool(ALLOW.search(blob) or INVARIANT.search(blob))


def main() -> int:
    violations: list[str] = []
    for rel in GUARDED:
        path = ROOT / rel
        if not path.is_file():
            violations.append(f"missing guarded file: {rel}")
            continue
        text = path.read_text()
        if "clippy::unwrap_used" not in text:
            violations.append(f"{rel}: missing #![warn(clippy::unwrap_used, …)]")
        lines = text.splitlines()
        ranges = test_module_ranges(lines)

        def in_test(n: int) -> bool:
            return any(a <= n <= b for a, b in ranges)

        for n, line in enumerate(lines, 1):
            if in_test(n) or not HIT.search(line):
                continue
            if not permitted(lines, n):
                violations.append(f"{rel}:{n}: unguarded panic point: {line.strip()[:100]}")

    if violations:
        print("I/O unwrap policy violations:", file=sys.stderr)
        for v in violations:
            print(f"  {v}", file=sys.stderr)
        return 1
    print(f"ok: {len(GUARDED)} guarded I/O modules pass unwrap policy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

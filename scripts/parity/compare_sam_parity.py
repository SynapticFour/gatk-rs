#!/usr/bin/env python3
"""
Compare two SAM files for behavioral parity (PrintReads-style round-trip).

- Ignores @PG and @CO (implementation-specific).
- Canonicalizes @HD VN (rust-htslib may echo 1.5 from BAM; GATK PrintWrites 1.6) — not behavioral.
- Compares sorted @HD, @SQ, @RG lines and sorted alignment records.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

# GATK 4 / HTSJDK SAM text uses VN:1.6; BAM round-trips may retain VN:1.5 — ignore for parity.
_CANON_HD_VN = "1.6"


def canonicalize_hd_line(line: str) -> str:
    if not line.startswith("@HD"):
        return line
    parts = line.split("\t")
    out: List[str] = []
    for p in parts:
        if p.startswith("VN:"):
            out.append(f"VN:{_CANON_HD_VN}")
        else:
            out.append(p)
    return "\t".join(out)


def split_sam_lines(text: str) -> List[str]:
    return [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]


def normalize_for_parity(lines: List[str]) -> dict:
    hd = sorted(
        canonicalize_hd_line(ln) for ln in lines if ln.startswith("@HD")
    )
    sq = sorted(ln for ln in lines if ln.startswith("@SQ"))
    rg = sorted(ln for ln in lines if ln.startswith("@RG"))
    align = sorted(ln for ln in lines if ln and not ln.startswith("@"))
    return {"hd": hd, "sq": sq, "rg": rg, "align": align}


def _attach_compact_mismatch(
    result: Dict[str, Any], ja: dict, rb: dict, max_sample: int = 5
) -> None:
    """Attach small diff fields; full normalized dicts only if total line count is modest."""
    for key in ("hd", "sq", "rg", "align"):
        if ja.get(key) != rb.get(key):
            a_list = ja.get(key) or []
            b_list = rb.get(key) or []
            result[f"mismatch_{key}_counts"] = (len(a_list), len(b_list))
            sa = set(a_list)
            sb = set(b_list)
            only_a = sorted(sa - sb)[:max_sample]
            only_b = sorted(sb - sa)[:max_sample]
            if only_a or only_b:
                result[f"mismatch_{key}_only_java_sample"] = only_a
                result[f"mismatch_{key}_only_rust_sample"] = only_b
    total_lines = sum(len(ja.get(k) or []) for k in ("hd", "sq", "rg", "align"))
    if total_lines <= 500:
        result["java_normalized"] = ja
        result["rust_normalized"] = rb
    else:
        result["note"] = (
            "Omitted full java_normalized/rust_normalized (too many lines); "
            "use mismatch_*_sample and SAM files to debug."
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--java-sam", required=True)
    parser.add_argument("--rust-sam", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--json-out", required=True)
    args = parser.parse_args()

    java_text = Path(args.java_sam).read_text(encoding="utf-8", errors="replace")
    rust_text = Path(args.rust_sam).read_text(encoding="utf-8", errors="replace")

    java_lines = split_sam_lines(java_text)
    rust_lines = split_sam_lines(rust_text)

    # Drop noisy header lines not meant for cross-tool parity
    drop_prefixes = ("@PG", "@CO")
    java_f = [ln for ln in java_lines if not ln.startswith(drop_prefixes)]
    rust_f = [ln for ln in rust_lines if not ln.startswith(drop_prefixes)]

    ja = normalize_for_parity(java_f)
    rb = normalize_for_parity(rust_f)
    equal = ja == rb

    result = {
        "label": args.label,
        "mode": "sam-file-parity",
        "equal": equal,
        "java_sam": args.java_sam,
        "rust_sam": args.rust_sam,
    }
    if not equal:
        result["reason"] = "sam_normalized_mismatch"
        # Avoid multi‑MB JSON when whole-genome SAMs diverge: keep a compact diff hint.
        _attach_compact_mismatch(result, ja, rb)

    Path(args.json_out).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if equal else 1


if __name__ == "__main__":
    raise SystemExit(main())

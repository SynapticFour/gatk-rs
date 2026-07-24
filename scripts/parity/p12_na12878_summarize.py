#!/usr/bin/env python3
"""Emit P12 NA12878 real-world JSON/MD with harness vs variant-parity fields."""
from __future__ import annotations

import json
import pathlib
import sys


def variants(path: pathlib.Path) -> list[tuple[str, str, str, str]]:
    if not path.exists():
        return []
    out: list[tuple[str, str, str, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) >= 5:
            out.append((c[0], c[1], c[3], c[4]))
    return out


def parity_status(jset: set, rset: set, java_exit: int, rust_exit: int) -> str:
    if java_exit != 0 or rust_exit != 0:
        return "tool_error"
    if jset == rset:
        return "variant_parity"
    return "variant_mismatch"


def main() -> None:
    if len(sys.argv) < 8:
        print(
            "usage: p12_na12878_summarize.py <json_out> <md_out> <java_vcf> <rust_vcf> "
            "<java_exit> <rust_exit> <notes_json>",
            file=sys.stderr,
        )
        sys.exit(2)
    json_out = pathlib.Path(sys.argv[1])
    md_out = pathlib.Path(sys.argv[2])
    java_vcf = pathlib.Path(sys.argv[3])
    rust_vcf = pathlib.Path(sys.argv[4])
    java_exit = int(sys.argv[5])
    rust_exit = int(sys.argv[6])
    notes = json.loads(sys.argv[7]) if sys.argv[7].strip() else {}
    diff_dir: pathlib.Path | None = None
    if len(sys.argv) >= 9 and sys.argv[8].strip():
        diff_dir = pathlib.Path(sys.argv[8])

    j = variants(java_vcf)
    r = variants(rust_vcf)
    jset = set(j)
    rset = set(r)
    shared = len(jset & rset)
    harness_ok = java_exit == 0 and rust_exit == 0
    status = "pass" if harness_ok else "tool_error"
    pstat = parity_status(jset, rset, java_exit, rust_exit)
    mode = notes.get("mode", "")
    read_augment = notes.get("read_augment", "")
    augment_line = (
        f"- read augment: `{read_augment}`"
        if read_augment
        else ""
    )
    if mode == "rust_only_reuse_java_vcf":
        harness_line = f"- status (harness): **{status}** — Rust run; Java VCF reused (not re-run)"
    else:
        harness_line = f"- status (harness): **{status}** — both tools exited successfully"

    payload = {
        "label": "phase12-realworld-na12878-20k",
        "status": status,
        "parity_status": pstat,
        "java_exit": java_exit,
        "rust_exit": rust_exit,
        "java_variant_count": len(j),
        "rust_variant_count": len(r),
        "shared_variant_count": shared,
        "java_only": len(jset - rset),
        "rust_only": len(rset - jset),
        "java_vcf": str(java_vcf),
        "rust_vcf": str(rust_vcf),
        "notes": notes,
    }
    if diff_dir is not None:
        diff_dir.mkdir(parents=True, exist_ok=True)
        java_only_path = diff_dir / "p12_java_only_diff.tsv"
        rust_only_path = diff_dir / "p12_rust_only.tsv"
        shared_path = diff_dir / "p12_shared.tsv"

        def write_keys(path: pathlib.Path, keys: set[tuple[str, str, str, str]]) -> None:
            lines = ["chrom\tpos\tref\talt"]
            for chrom, pos, ref, alt in sorted(keys, key=lambda k: (k[0], int(k[1]), k[2], k[3])):
                lines.append(f"{chrom}\t{pos}\t{ref}\t{alt}")
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        write_keys(java_only_path, jset - rset)
        write_keys(rust_only_path, rset - jset)
        write_keys(shared_path, jset & rset)
        payload["diff_tsv_dir"] = str(diff_dir)
        payload["java_only_tsv"] = str(java_only_path)
        payload["rust_only_tsv"] = str(rust_only_path)

    json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_out.write_text(
        "\n".join(
            [
                "# P12 Real-world NA12878 20k",
                "",
                harness_line,
                f"- parity_status: **{pstat}** — exact set match of (CHROM,POS,REF,ALT)",
                f"- java exit / rust exit: `{java_exit} / {rust_exit}`",
                f"- java variants: `{len(j)}`",
                f"- rust variants: `{len(r)}`",
                f"- shared variants: `{shared}`",
                f"- java-only / rust-only: `{len(jset - rset)} / {len(rset - jset)}`",
            ]
            + ([augment_line] if augment_line else [])
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        f"[p12-realworld] status={status} parity={pstat} java={len(j)} rust={len(r)} shared={shared}",
        flush=True,
    )
    if diff_dir is not None:
        print(
            f"[p12-realworld] diff: java_only={len(jset - rset)} rust_only={len(rset - jset)} "
            f"→ {diff_dir}/p12_{{java_only,rust_only,shared}}.tsv",
            flush=True,
        )


if __name__ == "__main__":
    main()

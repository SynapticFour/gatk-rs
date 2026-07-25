#!/usr/bin/env python3
"""Sprint N — excellence quality gates for gatk-haplotypecaller.

Gates:
  N-1  Coordinate allowlist — no new prod ``923*****`` outside frozen set
  N-2  Env allowlist — ``std::env::var`` only in runtime_config / parity_harness / dumps
  N-3  Module size — prod ``.rs`` ≤ 2500 lines (grandfather ceilings for known giants)
  N-4  Trailing-bool API — ban legacy double-bool sync/filter call shapes
  N-5  Production unwrap — no ``.unwrap()`` / ``.expect(`` in I/O modules outside tests
  N-6  Doc-claim — deferred + harness + oracle audits + CLAIM_MATRIX link check
  N-7  P12 band freeze — START/END inclusive widths must not widen vs freeze JSON

Usage:
  python3 scripts/parity/excellence_gates_audit.py
  python3 scripts/parity/excellence_gates_audit.py --gate n1
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "gatk-haplotypecaller/src"
COORD_ALLOWLIST = ROOT / "scripts/parity/coord_allowlist.json"
BAND_FREEZE = ROOT / "scripts/parity/p12_band_freeze.json"
CLAIM_MATRIX = ROOT / "docs/CLAIM_MATRIX.md"
HUB = ROOT / "docs/ARCHITECTURE.md"
ROOT_README = ROOT / "README.md"
SITE_SEMANTICS = SRC / "compatibility/java_hc_site_semantics.rs"
CONST_U64_RE = re.compile(r"pub const ([A-Z0-9_]+):\s*u64\s*=\s*(\d+);")

COORD_RE = re.compile(r"923[0-9]{5}")
ENV_VAR_RE = re.compile(r"\bstd::env::var(?:_os)?\s*\(")
UNWRAP_RE = re.compile(r"\.(?:unwrap|expect)\s*\(")
# Legacy shapes retired in Sprint I/K (must not return).
TRAILING_BOOL_PATTERNS = [
    re.compile(
        r"sync_assembly_events_from_haplotype_cigars_with_harvest\s*\([^;]*?,\s*(?:true|false)\s*,\s*(?:true|false)\s*\)",
        re.S,
    ),
    re.compile(
        r"filter_assembly_and_likelihoods\s*\([^;]*?,\s*(?:true|false)\s*,",
        re.S,
    ),
]

# N-3: L14-C emptied — all production modules under default 2500.
SIZE_GRANDFATHER: dict[str, int] = {}
SIZE_MAX_DEFAULT = 2500

# N-2 allowlisted relative paths (under src/) for direct env reads.
ENV_ALLOW_PREFIXES = (
    "parity_harness.rs",
    "runtime_config.rs",
)
ENV_ALLOW_SUFFIX = "_dump.rs"

# N-5 production I/O modules (strip cfg(test) / mod tests before scanning).
IO_MODULES = (
    "run.rs",
    "engine.rs",
    "region_vcf_emit.rs",
    "reference_vcf_emit.rs",
    "gvcf_writer.rs",
)


def rel(path: pathlib.Path) -> str:
    return str(path.relative_to(SRC))


def strip_test_modules(text: str) -> str:
    """Remove #[cfg(test)] mod … { … } (attrs allowed between cfg and mod)."""
    # Drop cfg(test) modules with balanced braces. Production modules often look like:
    #   #[cfg(test)]
    #   #[allow(...)]
    #   mod tests { ... }
    out = []
    i = 0
    while i < len(text):
        m = re.search(
            r"#\[cfg\(test\)\](?:\s*#\[[^\]]*\]|\s*//[^\n]*)*\s*mod\s+\w+\s*\{",
            text[i:],
        )
        if not m:
            out.append(text[i:])
            break
        start = i + m.start()
        out.append(text[i:start])
        brace_start = i + m.end() - 1
        depth = 0
        j = brace_start
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    j += 1
                    break
            j += 1
        i = j
    return "".join(out)


def gate_n1() -> list[str]:
    failures: list[str] = []
    if not COORD_ALLOWLIST.is_file():
        return [f"missing {COORD_ALLOWLIST}"]
    data = json.loads(COORD_ALLOWLIST.read_text(encoding="utf-8"))
    allowed: dict[str, set[str]] = {
        k: set(v) for k, v in data.get("files", {}).items()
    }
    for path in sorted(SRC.rglob("*.rs")):
        text = path.read_text(encoding="utf-8", errors="replace")
        found = set(COORD_RE.findall(text))
        if not found:
            continue
        key = rel(path)
        if key not in allowed:
            # New file with P12 coords: only allowed under compatibility/
            if key.startswith("compatibility/"):
                failures.append(
                    f"N-1: new coords in {key} — add to coord_allowlist.json with rationale "
                    f"({sorted(found)[:8]}{'…' if len(found) > 8 else ''})"
                )
            else:
                failures.append(
                    f"N-1: {key} introduces 923***** outside compatibility/ "
                    f"({sorted(found)[:8]})"
                )
            continue
        extra = found - allowed[key]
        if extra:
            failures.append(
                f"N-1: {key} has new coords not in allowlist: {sorted(extra)[:12]}"
            )
    return failures


def gate_n2() -> list[str]:
    failures: list[str] = []
    for path in sorted(SRC.rglob("*.rs")):
        key = rel(path)
        if key in ENV_ALLOW_PREFIXES or key.endswith(ENV_ALLOW_SUFFIX):
            continue
        # Doc comments mentioning the API are fine.
        for i, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            stripped = line.lstrip()
            if stripped.startswith("//") or stripped.startswith("*"):
                continue
            if ENV_VAR_RE.search(line):
                failures.append(f"N-2: {key}:{i} direct std::env::var — use runtime_config/parity_harness")
    return failures


def gate_n3() -> list[str]:
    failures: list[str] = []
    for path in sorted(SRC.rglob("*.rs")):
        key = rel(path)
        n = sum(1 for _ in path.open(encoding="utf-8", errors="replace"))
        ceiling = SIZE_GRANDFATHER.get(key, SIZE_MAX_DEFAULT)
        if n > ceiling:
            failures.append(f"N-3: {key} has {n} lines > ceiling {ceiling}")
    return failures


def gate_n4() -> list[str]:
    failures: list[str] = []
    for path in sorted(SRC.rglob("*.rs")):
        text = path.read_text(encoding="utf-8", errors="replace")
        # Skip unit-test bodies pulled via #[path] under tests/ — only scan src.
        for pat in TRAILING_BOOL_PATTERNS:
            if pat.search(text):
                failures.append(f"N-4: {rel(path)} matches banned trailing-bool call shape ({pat.pattern[:48]}…)")
    return failures


def gate_n5() -> list[str]:
    failures: list[str] = []
    for rel_path in IO_MODULES:
        path = SRC / rel_path
        if not path.is_file():
            failures.append(f"N-5: missing I/O module {rel_path}")
            continue
        text = strip_test_modules(path.read_text(encoding="utf-8", errors="replace"))
        for i, line in enumerate(text.splitlines(), 1):
            stripped = line.lstrip()
            if stripped.startswith("//"):
                continue
            if UNWRAP_RE.search(line):
                failures.append(f"N-5: {rel_path}:{i} unwrap/expect in production I/O path")
    return failures


def _check_md_links(md_path: pathlib.Path) -> list[str]:
    failures: list[str] = []
    if not md_path.is_file():
        return [f"missing {md_path}"]
    text = md_path.read_text(encoding="utf-8")
    # [text](relative.md) or [text](./foo.md#anchor)
    for m in re.finditer(r"\[[^\]]+\]\(([^)]+)\)", text):
        target = m.group(1).strip()
        if target.startswith(("http://", "https://", "mailto:")):
            continue
        path_part = target.split("#", 1)[0]
        if not path_part or path_part.startswith("/"):
            continue
        resolved = (md_path.parent / path_part).resolve()
        if not resolved.exists():
            failures.append(f"N-6: broken link in {md_path.relative_to(ROOT)} → {target}")
    return failures


def gate_n6() -> list[str]:
    failures: list[str] = []
    # Sub-audits
    for script, label in [
        ("scripts/parity/deferred_features_audit.py", "deferred-features"),
        ("scripts/parity/oracle_emit_audit.py", "oracle-emit"),
    ]:
        proc = subprocess.run(
            [sys.executable, str(ROOT / script)],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            failures.append(f"N-6: {label} audit failed:\n{proc.stderr or proc.stdout}")
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts/parity/p12_site_id_audit.py"), "--check-harness"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        failures.append(f"N-6: harness/site-id audit failed:\n{proc.stderr or proc.stdout}")

    if not CLAIM_MATRIX.is_file():
        failures.append("N-6: missing CLAIM_MATRIX.md")
    else:
        failures.extend(_check_md_links(CLAIM_MATRIX))
        claim_text = CLAIM_MATRIX.read_text(encoding="utf-8")
        for needle in ("W-H1", "W-H3", "genome-wide", "L6"):
            if needle not in claim_text:
                failures.append(f"N-6: CLAIM_MATRIX.md missing expected term {needle!r}")

    for path, needle in [
        (HUB, "CLAIM_MATRIX.md"),
        (ROOT_README, "CLAIM_MATRIX.md"),
        (ROOT_README, "Why does this exist"),
    ]:
        if not path.is_file():
            failures.append(f"N-6: missing {path}")
            continue
        if needle not in path.read_text(encoding="utf-8"):
            failures.append(f"N-6: {path.relative_to(ROOT)} must reference {needle!r}")
    return failures


def gate_n7() -> list[str]:
    """P12 band freeze: inclusive widths must not grow (generalize rule #4)."""
    failures: list[str] = []
    if not BAND_FREEZE.is_file():
        return ["N-7: missing scripts/parity/p12_band_freeze.json"]
    if not SITE_SEMANTICS.is_file():
        return ["N-7: missing java_hc_site_semantics.rs"]
    freeze = json.loads(BAND_FREEZE.read_text(encoding="utf-8"))
    consts = {
        m.group(1): int(m.group(2))
        for m in CONST_U64_RE.finditer(SITE_SEMANTICS.read_text(encoding="utf-8"))
    }
    for band in freeze.get("bands", []):
        sc = band["start_const"]
        ec = band["end_const"]
        if sc not in consts or ec not in consts:
            failures.append(f"N-7: freeze band missing const {sc}/{ec}")
            continue
        start, end = consts[sc], consts[ec]
        width = end - start + 1
        frozen_w = int(band["inclusive_width"])
        if width > frozen_w:
            failures.append(
                f"N-7: band widened {sc}..{ec}: width {width} > frozen {frozen_w} "
                f"(never widen P12 bands; see GENERALIZE_WITHOUT_OVERFIT.md)"
            )
        if start != int(band["start"]) or end != int(band["end"]):
            # Allow shrink (start↑ or end↓) without failing when width ≤ frozen.
            if width > frozen_w:
                continue
            # Record intentional shrinks as OK; reject expansions of either edge
            # that keep width equal via shifting (still a form of widening coverage).
            if start < int(band["start"]) or end > int(band["end"]):
                failures.append(
                    f"N-7: band edge expanded {sc}..{ec}: "
                    f"[{start},{end}] vs frozen [{band['start']},{band['end']}]"
                )
    return failures


GATES = {
    "n1": ("coordinate allowlist", gate_n1),
    "n2": ("env allowlist", gate_n2),
    "n3": ("module size", gate_n3),
    "n4": ("trailing-bool API", gate_n4),
    "n5": ("production unwrap", gate_n5),
    "n6": ("doc-claim CI", gate_n6),
    "n7": ("p12 band freeze", gate_n7),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gate",
        choices=sorted(GATES) + ["all"],
        default="all",
        help="Run a single gate or all (default)",
    )
    args = parser.parse_args()
    selected = list(GATES) if args.gate == "all" else [args.gate]
    all_failures: list[str] = []
    for gid in selected:
        name, fn = GATES[gid]
        fails = fn()
        if fails:
            print(f"[excellence-gates] {gid.upper()} FAIL ({name})", file=sys.stderr)
            for f in fails:
                print(f"  {f}", file=sys.stderr)
            all_failures.extend(fails)
        else:
            print(f"[excellence-gates] {gid.upper()} PASS ({name})")
    if all_failures:
        print(f"[excellence-gates] FAIL ({len(all_failures)} issue(s))", file=sys.stderr)
        return 1
    print("[excellence-gates] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

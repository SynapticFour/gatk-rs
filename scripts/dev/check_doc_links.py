#!/usr/bin/env python3
"""Fail if .rs / .md files reference missing repo paths under docs/ or *.md.

Catches:
  - Markdown links ``[text](relative/or/docs/path)``
  - Path-like tokens ``docs/...`` (repo-root relative)
  - Repo-style ``*.md`` paths (``parity/...``, ``gatk-*/...``, ``scripts/...``, …)

Skips http(s)/mailto, anchors-only, incomplete template stubs, and generated trees.
A reference is OK if it resolves from the repo root **or** relative to the source file.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

MD_LINK = re.compile(r"(?<!!)\[[^\]]*\]\(([^)]+)\)")
DOCS_PATH = re.compile(r"(?<![A-Za-z0-9_./-])(docs/[A-Za-z0-9_./@+-]+)")
MD_FILE = re.compile(
    r"(?<![A-Za-z0-9_./-])((?:(?:\.\./)+|\./)?(?:[A-Za-z0-9_.-]+/)+[A-Za-z0-9_.-]+\.md)"
)

SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
REPO_PREFIXES = (
    "docs/",
    "parity/",
    "scripts/",
    "tools/",
    "fuzz/",
    "gatk-",
    ".github/",
    ".quality-gates/",
)
FILE_SUFFIXES = (".md", ".env", ".txt", ".json", ".yml", ".yaml", ".toml", ".rs", ".sh", ".py")


def should_skip(path: Path) -> bool:
    if any(part in {"target", ".git", "node_modules"} for part in path.parts):
        return True
    rel = path.relative_to(ROOT).as_posix()
    if rel.startswith("docs/parity-site/data/"):
        return True
    return False


def iter_sources() -> list[Path]:
    out: list[Path] = []
    for pattern in ("**/*.rs", "**/*.md"):
        for p in ROOT.glob(pattern):
            if not should_skip(p):
                out.append(p)
    return sorted(set(out))


def clean_href(href: str) -> str:
    href = href.strip()
    if href.startswith("<") and ">" in href:
        href = href[1 : href.index(">")].strip()
    if " " in href:
        href = href.split(" ", 1)[0]
    href = href.split("#", 1)[0].split("?", 1)[0].strip()
    # Trailing sentence punctuation commonly glued in prose / rustdoc.
    href = href.rstrip(".,;:)")
    return href


GENERATED_PREFIXES = (
    "docs/perf/runs/",
    "parity/reports/",
    "parity/build/",
    "parity/giab/runs/",
    "parity/giab/truth/",
)


def is_generated_href(href: str) -> bool:
    h = clean_href(href).lstrip("./")
    while h.startswith("../"):
        h = h[3:]
    return any(h.startswith(prefix) for prefix in GENERATED_PREFIXES)


def looks_complete(href: str) -> bool:
    if not href or href in {".", ".."}:
        return False
    if "<" in href or ">" in href or "*" in href:
        return False
    # Template stubs like docs/perf/runs/dedicated_
    if href.endswith("_") or href.endswith("/"):
        return False
    if href.startswith("docs/"):
        return href.endswith(FILE_SUFFIXES) or bool(re.search(r"/[A-Za-z0-9_.-]+\.[A-Za-z0-9]+$", href))
    return True


def candidate_paths(source: Path, href: str) -> list[Path]:
    href = clean_href(href)
    if not looks_complete(href) or SCHEME.match(href) or href.startswith("//"):
        return []
    if href.startswith("/"):
        href = href.lstrip("/")

    out: list[Path] = []
    if href.startswith(REPO_PREFIXES) or href.startswith("docs/"):
        out.append(ROOT / href)
    if href.startswith("./") or href.startswith("../"):
        out.append((source.parent / href).resolve())
    elif href.endswith(".md") or "/" in href:
        # Prefer repo-root for workspace-style paths; also try relative.
        out.append(ROOT / href)
        out.append((source.parent / href).resolve())

    # Dedup while preserving order.
    seen: set[Path] = set()
    uniq: list[Path] = []
    for p in out:
        try:
            rp = p.resolve()
        except OSError:
            continue
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append(rp)
    return uniq


def in_repo(path: Path) -> bool:
    try:
        path.resolve().relative_to(ROOT.resolve())
        return True
    except (ValueError, OSError):
        return False


def collect_raw_refs(line: str) -> list[str]:
    refs: list[str] = []
    for m in MD_LINK.finditer(line):
        refs.append(m.group(1))
    for m in DOCS_PATH.finditer(line):
        refs.append(m.group(1))
    for m in MD_FILE.finditer(line):
        refs.append(m.group(1))
    return refs


def main() -> int:
    dead: list[str] = []
    checked = 0
    sources = iter_sources()
    for source in sources:
        text = source.read_text(encoding="utf-8", errors="replace")
        rel_src = source.relative_to(ROOT).as_posix()
        for lineno, line in enumerate(text.splitlines(), start=1):
            for raw in collect_raw_refs(line):
                if is_generated_href(raw):
                    continue
                candidates = candidate_paths(source, raw)
                if not candidates:
                    continue
                checked += 1
                if any(in_repo(c) and c.exists() for c in candidates):
                    continue
                shown = clean_href(raw)
                dead.append(f"{rel_src}:{lineno}: {shown!r} (not found from repo root or {rel_src})")

    if dead:
        # Dedup identical messages (md link + docs path double-match).
        uniq = list(dict.fromkeys(dead))
        print(f"check_doc_links: {len(uniq)} dead reference(s) (checked {checked}):", file=sys.stderr)
        for row in uniq:
            print(f"  {row}", file=sys.stderr)
        return 1

    print(f"check_doc_links: ok ({checked} references in {len(sources)} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

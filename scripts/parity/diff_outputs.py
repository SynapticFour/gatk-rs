#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path
from typing import Optional


def normalize(text: str) -> str:
    text = text.replace("\r\n", "\n")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def extract_with_regex(text: str, pattern: Optional[str]) -> str:
    if not pattern:
        return text
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if not matches:
        return ""
    if isinstance(matches[0], tuple):
        return "\n".join(" ".join(part for part in m if part) for m in matches)
    return "\n".join(matches)


EXTRACT_PREVIEW_MAX = 4096


def preview_extracted(text: str) -> str:
    if len(text) <= EXTRACT_PREVIEW_MAX:
        return text
    return text[:EXTRACT_PREVIEW_MAX] + "\n... [truncated] ..."


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--java", required=True)
    parser.add_argument("--rust", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--mode", choices=["strict", "normalized"], default="normalized")
    parser.add_argument("--json-out", required=True)
    parser.add_argument("--extract-regex", default="")
    parser.add_argument("--presence-only", action="store_true")
    args = parser.parse_args()

    java_text = Path(args.java).read_text(encoding="utf-8", errors="replace")
    rust_text = Path(args.rust).read_text(encoding="utf-8", errors="replace")
    java_text = extract_with_regex(java_text, args.extract_regex or None)
    rust_text = extract_with_regex(rust_text, args.extract_regex or None)

    if args.presence_only and args.extract_regex:
        equal = bool(java_text.strip()) and bool(rust_text.strip())
    elif args.mode == "strict":
        equal = java_text == rust_text
    else:
        equal = normalize(java_text) == normalize(rust_text)

    result = {
        "label": args.label,
        "mode": args.mode,
        "equal": equal,
        "extract_regex": args.extract_regex or None,
        "presence_only": args.presence_only,
        "java_output": args.java,
        "rust_output": args.rust,
    }
    if args.extract_regex:
        result["java_extracted"] = preview_extracted(java_text)
        result["rust_extracted"] = preview_extracted(rust_text)
    Path(args.json_out).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if equal else 1


if __name__ == "__main__":
    raise SystemExit(main())

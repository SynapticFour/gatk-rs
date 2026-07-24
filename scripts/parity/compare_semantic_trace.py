#!/usr/bin/env python3
"""Compare two semantic-trace NDJSON files and report the first divergence.

Schema: gatk_rs.hc.semantic_trace/v1

Typical usage:
  # Rust vs Rust (regression)
  python3 scripts/parity/compare_semantic_trace.py \\
      --left artifacts/rust.ndjson --right artifacts/rust2.ndjson

  # Java (projected) vs Rust
  python3 scripts/parity/project_java_to_semantic_trace.py \\
      --vcf java.vcf --regions java_regions.tsv -o artifacts/java.ndjson
  python3 scripts/parity/compare_semantic_trace.py \\
      --left artifacts/java.ndjson --right artifacts/rust.ndjson \\
      --json-out first_divergence.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

SCHEMA = "gatk_rs.hc.semantic_trace/v1"

STAGE_ORDER = {
    "activity_profile": 0,
    "active_region": 1,
    "assembly_graph": 2,
    "reference_path": 3,
    "candidate_haplotypes": 4,
    "read_likelihoods": 5,
    "genotype_likelihoods": 6,
    "inactive_rcm": 7,
    "vcf_emission": 8,
}


def load_events(path: Path) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            ev = json.loads(line)
        except json.JSONDecodeError as exc:
            raise SystemExit(f"{path}:{line_no}: invalid JSON: {exc}") from exc
        if ev.get("schema") != SCHEMA:
            raise SystemExit(
                f"{path}:{line_no}: unexpected schema {ev.get('schema')!r} (want {SCHEMA})"
            )
        events.append(ev)
    return events


def region_key(ev: dict[str, Any]) -> tuple[Any, ...]:
    r = ev.get("region") or {}
    return (
        r.get("contig"),
        r.get("start"),
        r.get("end"),
        r.get("is_active"),
        ev.get("stage"),
    )


def sort_key(ev: dict[str, Any]) -> tuple[Any, ...]:
    r = ev.get("region") or {}
    stage = ev.get("stage")
    return (
        r.get("contig") or "",
        r.get("start") if r.get("start") is not None else -1,
        r.get("end") if r.get("end") is not None else -1,
        STAGE_ORDER.get(stage, 99),
        ev.get("seq", 0),
    )


def normalize_payload(payload: Any, float_tol: float) -> Any:
    """Round floats for tolerant compare; leave structure intact."""
    if isinstance(payload, float):
        if abs(payload) == float("inf"):
            return payload
        # Quantize to near-micro precision unless caller uses looser tol via equality helper.
        return round(payload, 6)
    if isinstance(payload, list):
        return [normalize_payload(x, float_tol) for x in payload]
    if isinstance(payload, dict):
        return {k: normalize_payload(v, float_tol) for k, v in sorted(payload.items())}
    return payload


def payloads_equal(a: Any, b: Any, float_tol: float) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        if a != a and b != b:  # NaN
            return True
        return abs(a - b) <= float_tol
    if type(a) is not type(b):
        # int/float interchange
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            return abs(float(a) - float(b)) <= float_tol
        return False
    if isinstance(a, dict):
        if set(a) != set(b):
            return False
        return all(payloads_equal(a[k], b[k], float_tol) for k in a)
    if isinstance(a, list):
        if len(a) != len(b):
            return False
        return all(payloads_equal(x, y, float_tol) for x, y in zip(a, b))
    return a == b


def first_divergence(
    left: list[dict[str, Any]],
    right: list[dict[str, Any]],
    *,
    float_tol: float,
    stages: set[str] | None,
) -> dict[str, Any]:
    left_f = [e for e in left if stages is None or e.get("stage") in stages]
    right_f = [e for e in right if stages is None or e.get("stage") in stages]
    left_f.sort(key=sort_key)
    right_f.sort(key=sort_key)

    n = min(len(left_f), len(right_f))
    for i in range(n):
        le, re = left_f[i], right_f[i]
        lk, rk = region_key(le), region_key(re)
        if lk != rk:
            return {
                "status": "diverged",
                "index": i,
                "reason": "region_or_stage_mismatch",
                "left": {"key": list(lk), "seq": le.get("seq"), "impl": le.get("impl")},
                "right": {"key": list(rk), "seq": re.get("seq"), "impl": re.get("impl")},
                "left_event": le,
                "right_event": re,
            }
        lp = normalize_payload(le.get("payload"), float_tol)
        rp = normalize_payload(re.get("payload"), float_tol)
        if not payloads_equal(lp, rp, float_tol):
            return {
                "status": "diverged",
                "index": i,
                "reason": "payload_mismatch",
                "stage": le.get("stage"),
                "region": le.get("region"),
                "left": {"seq": le.get("seq"), "impl": le.get("impl"), "payload": lp},
                "right": {"seq": re.get("seq"), "impl": re.get("impl"), "payload": rp},
            }

    if len(left_f) != len(right_f):
        longer = "left" if len(left_f) > len(right_f) else "right"
        extra = left_f[n] if longer == "left" else right_f[n]
        return {
            "status": "diverged",
            "index": n,
            "reason": "length_mismatch",
            "longer": longer,
            "left_count": len(left_f),
            "right_count": len(right_f),
            "first_extra": {
                "key": list(region_key(extra)),
                "seq": extra.get("seq"),
                "impl": extra.get("impl"),
                "stage": extra.get("stage"),
            },
        }

    return {
        "status": "identical",
        "compared_events": n,
        "stages": sorted(stages) if stages else "all",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--left", required=True, type=Path, help="Reference / Java NDJSON")
    p.add_argument("--right", required=True, type=Path, help="Candidate / Rust NDJSON")
    p.add_argument("--json-out", type=Path, help="Write full report JSON")
    p.add_argument(
        "--stages",
        nargs="*",
        help="Restrict compare to these stage names (default: all)",
    )
    p.add_argument(
        "--float-tol",
        type=float,
        default=1e-5,
        help="Absolute float tolerance after round-6 normalize",
    )
    args = p.parse_args()

    left = load_events(args.left)
    right = load_events(args.right)
    stages = set(args.stages) if args.stages else None
    report = first_divergence(left, right, float_tol=args.float_tol, stages=stages)
    report["left_path"] = str(args.left)
    report["right_path"] = str(args.right)
    report["left_total"] = len(left)
    report["right_total"] = len(right)

    text = json.dumps(report, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text + "\n", encoding="utf-8")

    if report["status"] == "identical":
        print(
            f"OK: traces match ({report['compared_events']} events)",
            file=sys.stderr,
        )
        return 0

    print("FIRST DIVERGENCE", file=sys.stderr)
    print(
        f"  reason={report['reason']} index={report.get('index')}",
        file=sys.stderr,
    )
    if "stage" in report:
        print(f"  stage={report['stage']} region={report.get('region')}", file=sys.stderr)
    elif "left" in report and "key" in report["left"]:
        print(f"  left_key={report['left']['key']}", file=sys.stderr)
        print(f"  right_key={report['right']['key']}", file=sys.stderr)
    if args.json_out:
        print(f"  wrote {args.json_out}", file=sys.stderr)
    else:
        # Compact stdout for scripting
        print(text)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

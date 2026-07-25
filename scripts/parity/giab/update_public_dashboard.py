#!/usr/bin/env python3
"""Ingest nightly / genomewide JSON into docs/parity-site/data/history.json.

Also writes data/latest.json for convenience. The static HTML (index.html +
app.js) reads history.json at runtime — no bundler required.

Usage:
  # From nightly happy_summary.json
  python3 scripts/parity/giab/update_public_dashboard.py \\
    --source nightly \\
    --json parity/reports/.../happy_summary.json \\
    --site-dir docs/parity-site \\
    --commit-sha \"$GITHUB_SHA\" \\
    --workflow-run-url \"$RUN_URL\"

  # From genomewide run directory (samples.jsonl + SCOPE)
  python3 scripts/parity/giab/update_public_dashboard.py \\
    --source genomewide \\
    --run-dir parity/giab/runs/genomewide \\
    --site-dir docs/parity-site \\
    --commit-sha \"$GITHUB_SHA\"
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from datetime import datetime, timezone
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[3]
PINNED_ENV = ROOT / "docs" / "GATK_PINNED.env"


def load_pinned() -> dict[str, str]:
    out: dict[str, str] = {}
    if not PINNED_ENV.is_file():
        return {
            "GATK_PINNED_REF": "4.4.0.0",
            "GATK_PINNED_SHA": "unknown",
            "GATK_DOCKER_IMAGE": "us.gcr.io/broad-gatk/gatk:4.4.0.0",
        }
    for line in PINNED_ENV.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def load_history(path: pathlib.Path) -> dict[str, Any]:
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"meta": {}, "runs": []}


def prf_from_block(block: dict[str, Any] | None) -> dict[str, float | None]:
    if not block:
        return {"precision": None, "recall": None, "f1": None}
    return {
        "precision": _f(block.get("precision")),
        "recall": _f(block.get("recall")),
        "f1": _f(block.get("f1")),
    }


def _f(v: Any) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def pick_all_stratum(metrics: list[dict[str, Any]], label_substr: str) -> dict[str, Any] | None:
    """Prefer stratum ALL / empty for a query label (rust/java)."""
    candidates = [
        m
        for m in metrics
        if label_substr.lower() in str(m.get("query_label", "")).lower()
    ]
    if not candidates:
        candidates = list(metrics)
    for m in candidates:
        strat = str(m.get("stratum", "")).upper()
        if strat in ("ALL", "*", "", "NONE"):
            return m
    return candidates[0] if candidates else None


def ingest_nightly(summary: dict[str, Any], pinned: dict[str, str], args: argparse.Namespace) -> dict[str, Any]:
    regions = []
    metrics = []
    for row in summary.get("regions", []):
        name = row.get("region") or "unknown"
        regions.append(name)
        for engine in ("rust", "java"):
            block = row.get(engine) or {}
            for vtype in ("SNP", "INDEL"):
                prf = prf_from_block(block.get(vtype))
                if prf["f1"] is None and prf["precision"] is None:
                    continue
                metrics.append(
                    {
                        "region": name,
                        "engine": engine,
                        "variant_type": vtype,
                        **prf,
                    }
                )

    return {
        "id": f"nightly-{summary.get('generated_utc') or args.commit_sha}",
        "workflow": "nightly-equivalence",
        "generated_utc": summary.get("generated_utc")
        or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "commit_sha": args.commit_sha or summary.get("commit_sha") or "unknown",
        "workflow_run_url": args.workflow_run_url or None,
        "scope": {
            "kind": "trio_joint_e2e",
            "pipeline": "HaplotypeCaller (GVCF) → CombineGVCFs → GenotypeGVCFs → VariantFiltration",
            "samples": ["HG002", "HG003", "HG004"],
            "regions": regions,
            "assembly": "hs37d5 (GRCh37)",
            "truth": "GIAB HG002 NISTv4.2.1 (joint callset scored on HG002)",
            "eval_engine": "Illumina hap.py (Docker) vs HG002 truth",
            "java_gatk_version": pinned.get("GATK_PINNED_REF", "4.4.0.0"),
            "java_gatk_sha": pinned.get("GATK_PINNED_SHA", "unknown"),
            "java_gatk_docker": pinned.get("GATK_DOCKER_IMAGE", ""),
            "honesty": (
                "Nightly trio E2E covers only the listed chromosomes/hard slices "
                "(not full WGS). BAMs are region-sliced; metrics are hap.py vs HG002 truth "
                "after joint genotyping + SNP hard-filters."
            ),
        },
        "metrics": metrics,
    }


def ingest_genomewide(run_dir: pathlib.Path, pinned: dict[str, str], args: argparse.Namespace) -> dict[str, Any]:
    jsonl = run_dir / "samples.jsonl"
    scope_txt = ""
    scope_path = run_dir / "SCOPE.txt"
    if scope_path.is_file():
        scope_txt = scope_path.read_text(encoding="utf-8").strip()

    samples: list[str] = []
    metrics: list[dict[str, Any]] = []
    mode = "unknown"
    if jsonl.is_file():
        for line in jsonl.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            sample = row.get("sample") or "sample"
            samples.append(sample)
            mode = row.get("mode") or mode
            eq = row.get("equiv_results") or {}
            for engine, key in (("rust", "rust_vs_truth"), ("java", "java_vs_truth")):
                picked = pick_all_stratum(eq.get(key) or [], engine)
                if not picked:
                    continue
                for vtype, field in (("SNP", "snp"), ("INDEL", "indel")):
                    prf = prf_from_block(picked.get(field))
                    metrics.append(
                        {
                            "region": f"{mode}:{sample}",
                            "sample": sample,
                            "engine": engine,
                            "variant_type": vtype,
                            **prf,
                        }
                    )

    # Intervals file if present
    intervals = []
    iv_path = run_dir / "intervals.txt"
    if iv_path.is_file():
        intervals = [
            ln.strip()
            for ln in iv_path.read_text(encoding="utf-8").splitlines()
            if ln.strip()
        ]

    return {
        "id": f"genomewide-{args.commit_sha}-{mode}",
        "workflow": "genomewide-validation",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "commit_sha": args.commit_sha or "unknown",
        "workflow_run_url": args.workflow_run_url or None,
        "scope": {
            "kind": "hc_genomewide",
            "pipeline": "HaplotypeCaller (Java GATK 4.4 vs gatk-rs) → hap.py/RTG via gatk-rs-equiv",
            "samples": samples,
            "regions": intervals[:40] + (["…"] if len(intervals) > 40 else []),
            "mode_description": scope_txt or mode,
            "assembly": "hs37d5 (GRCh37)",
            "truth": "GIAB NISTv4.2.1 per sample (confident BED)",
            "eval_engine": "gatk-rs-equiv (hap.py preferred, RTG fallback)",
            "java_gatk_version": pinned.get("GATK_PINNED_REF", "4.4.0.0"),
            "java_gatk_sha": pinned.get("GATK_PINNED_SHA", "unknown"),
            "java_gatk_docker": pinned.get("GATK_DOCKER_IMAGE", ""),
            "honesty": (
                f"GIAB_MODE={mode}. Metrics are per-sample truth F1 from hap.py/RTG. "
                "Full `autosomes` is WGS-scale and runs on the paid self-hosted runner; "
                "`ci-subset` / window modes are not whole-genome claims."
            ),
        },
        "metrics": metrics,
    }


def ingest_cohort_scale(payload: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    """Accept a pre-built dashboard_run.json from run_joint_cohort_scale.sh."""
    run = dict(payload)
    run["commit_sha"] = args.commit_sha or run.get("commit_sha") or "unknown"
    if args.workflow_run_url:
        run["workflow_run_url"] = args.workflow_run_url
    run.setdefault("workflow", "joint-cohort-scale")
    run.setdefault("generated_utc", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
    return run


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--source",
        choices=("nightly", "genomewide", "cohort_scale"),
        required=True,
    )
    ap.add_argument(
        "--json",
        type=pathlib.Path,
        help="nightly happy_summary.json OR cohort_scale dashboard_run.json",
    )
    ap.add_argument("--run-dir", type=pathlib.Path, help="genomewide run directory")
    ap.add_argument(
        "--site-dir",
        type=pathlib.Path,
        default=ROOT / "docs" / "parity-site",
    )
    ap.add_argument("--commit-sha", default="unknown")
    ap.add_argument("--workflow-run-url", default="")
    ap.add_argument("--max-runs", type=int, default=60, help="Retain last N runs")
    args = ap.parse_args()

    pinned = load_pinned()
    site = args.site_dir
    data_dir = site / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    history_path = data_dir / "history.json"
    latest_path = data_dir / "latest.json"

    if args.source == "nightly":
        if not args.json or not args.json.is_file():
            print("[parity-site] missing --json happy_summary.json", file=sys.stderr)
            return 2
        summary = json.loads(args.json.read_text(encoding="utf-8"))
        run = ingest_nightly(summary, pinned, args)
    elif args.source == "cohort_scale":
        if not args.json or not args.json.is_file():
            print("[parity-site] missing --json dashboard_run.json", file=sys.stderr)
            return 2
        payload = json.loads(args.json.read_text(encoding="utf-8"))
        run = ingest_cohort_scale(payload, args)
    else:
        if not args.run_dir or not args.run_dir.is_dir():
            print("[parity-site] missing --run-dir", file=sys.stderr)
            return 2
        run = ingest_genomewide(args.run_dir, pinned, args)

    history = load_history(history_path)
    history.setdefault("runs", [])
    # De-dupe by id
    history["runs"] = [r for r in history["runs"] if r.get("id") != run["id"]]
    history["runs"].append(run)
    history["runs"] = history["runs"][-args.max_runs :]

    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    history["meta"] = {
        "title": "gatk-rs equivalence dashboard",
        "java_gatk_version": pinned.get("GATK_PINNED_REF", "4.4.0.0"),
        "java_gatk_sha": pinned.get("GATK_PINNED_SHA", "unknown"),
        "java_gatk_docker": pinned.get("GATK_DOCKER_IMAGE", ""),
        "updated_utc": now,
        "notes": (
            "Updated by nightly-equivalence / genomewide-validation / "
            "joint-cohort-scale gates."
        ),
    }

    history_path.write_text(json.dumps(history, indent=2) + "\n", encoding="utf-8")
    latest_path.write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    print(f"[parity-site] wrote {history_path} ({len(history['runs'])} runs)")
    print(f"[parity-site] wrote {latest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

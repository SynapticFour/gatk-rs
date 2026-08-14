#!/usr/bin/env python3
"""Compare Java vs Rust wall / Peak-RSS per GIAB ci-subset HC shard.

Pulls `giab-hc-*` artifacts from a GitHub Actions run (or reads a local unpack
dir) and joins `*.time.txt` (+ optional hc-mem-probe peak) by shard name.

Usage:
  # Download + summarize (needs `gh` auth):
  ./scripts/parity/giab/summarize_ci_subset_shard_perf.py --run-id 31577683299

  # From an already-downloaded artifact tree:
  ./scripts/parity/giab/summarize_ci_subset_shard_perf.py --artifact-dir /tmp/giab-31577683299

Outputs JSON (+ optional Markdown) ranking where Rust beats / lags Java.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Reuse the same GNU/macOS parsers as the harness.
_REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO / "scripts" / "perf"))
from parse_time_metrics import parse_time_log  # noqa: E402


SHARD_RE = re.compile(
    r"(?:^|/)(?P<engine>java|rust)\.(?P<shard>[^/]+)\.time\.txt$"
)
MEM_PROBE_RE = re.compile(
    r"hc-mem-probe-(?P<shard>.+)-(?P<engine>java|rust)\.log$"
)
RSS_KB_RE = re.compile(r"rss_kb=(\d+)")


def parse_time_file(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    m = parse_time_log(text)
    return {
        "wall_sec": m.get("wall_sec"),
        "max_rss_kb": m.get("peak_rss_kb"),
        "path": str(path),
    }


def probe_peak_rss_kb(path: Path) -> int | None:
    peak = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = RSS_KB_RE.search(line)
        if m:
            v = int(m.group(1))
            peak = v if peak is None else max(peak, v)
    return peak


def collect_from_dir(root: Path) -> dict[str, dict[str, dict]]:
    """shard -> engine -> metrics."""
    out: dict[str, dict[str, dict]] = {}
    for path in root.rglob("*.time.txt"):
        m = SHARD_RE.search(path.as_posix())
        if not m:
            continue
        engine, shard = m.group("engine"), m.group("shard")
        cell = out.setdefault(shard, {}).setdefault(engine, {})
        cell.update(parse_time_file(path))
    for path in root.rglob("hc-mem-probe-*.log"):
        m = MEM_PROBE_RE.search(path.name)
        if not m:
            continue
        engine, shard = m.group("engine"), m.group("shard")
        cell = out.setdefault(shard, {}).setdefault(engine, {})
        peak = probe_peak_rss_kb(path)
        if peak is not None:
            cell["probe_peak_rss_kb"] = peak
            cell["probe_path"] = str(path)
    return out


def download_run_artifacts(run_id: str, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    # Only HC shard artifacts (time + mem probe); ignore prepare/finalize bulk.
    meta = json.loads(
        subprocess.check_output(
            ["gh", "repo", "view", "--json", "nameWithOwner"],
            text=True,
        )
    )
    owner, repo = meta["nameWithOwner"].split("/", 1)
    names: list[str] = []
    page = 1
    while True:
        arts = json.loads(
            subprocess.check_output(
                [
                    "gh",
                    "api",
                    f"repos/{owner}/{repo}/actions/runs/{run_id}/artifacts"
                    f"?per_page=100&page={page}",
                ],
                text=True,
            )
        )
        batch = arts.get("artifacts") or []
        if not batch:
            break
        for a in batch:
            if a["name"].startswith("giab-hc-") and not a.get("expired"):
                names.append(a["name"])
        if len(batch) < 100:
            break
        page += 1
    if not names:
        # Fallback: download all (may be heavy once finalize lands).
        subprocess.check_call(
            ["gh", "run", "download", run_id, "-D", str(dest)],
        )
        return dest
    for name in names:
        subprocess.check_call(
            ["gh", "run", "download", run_id, "-n", name, "-D", str(dest / name)],
        )
    return dest


def fmt_sec(v: float | None) -> str:
    if v is None:
        return "n/a"
    if v >= 3600:
        h = int(v // 3600)
        m = int((v % 3600) // 60)
        s = v % 60
        return f"{h}h{m:02d}m{s:04.1f}s"
    if v >= 60:
        return f"{int(v // 60)}m{v % 60:04.1f}s"
    return f"{v:.1f}s"


def fmt_mib(kb: int | None) -> str:
    if kb is None:
        return "n/a"
    return f"{kb / 1024:.1f} MiB"


def prefer_rss_kb(cell: dict) -> int | None:
    """Prefer hc-mem-probe Peak over `/usr/bin/time` (JVM docker time RSS is ~30 MiB junk)."""
    probe = cell.get("probe_peak_rss_kb")
    timed = cell.get("max_rss_kb")
    if probe is not None and timed is not None:
        # Probe wins unless time RSS is clearly larger (rare; keep max for safety).
        return max(probe, timed)
    return probe if probe is not None else timed


def build_rows(by_shard: dict[str, dict[str, dict]]) -> list[dict]:
    rows = []
    for shard in sorted(by_shard):
        j = by_shard[shard].get("java") or {}
        r = by_shard[shard].get("rust") or {}
        jw, rw = j.get("wall_sec"), r.get("wall_sec")
        jr = prefer_rss_kb(j)
        rr = prefer_rss_kb(r)
        wall_ratio = (rw / jw) if (jw and rw and jw > 0) else None
        rss_ratio = (rr / jr) if (jr and rr and jr > 0) else None
        if wall_ratio is None:
            verdict = "incomplete"
        elif wall_ratio < 0.9:
            verdict = "rust_faster"
        elif wall_ratio > 1.1:
            verdict = "rust_slower"
        else:
            verdict = "wall_tie"
        rows.append(
            {
                "shard": shard,
                "java_wall_sec": jw,
                "rust_wall_sec": rw,
                "wall_ratio_rust_over_java": wall_ratio,
                "java_rss_kb": jr,
                "rust_rss_kb": rr,
                "rss_ratio_rust_over_java": rss_ratio,
                "java_probe_peak_rss_kb": j.get("probe_peak_rss_kb"),
                "rust_probe_peak_rss_kb": r.get("probe_peak_rss_kb"),
                "verdict": verdict,
            }
        )
    return rows


def summarize(rows: list[dict]) -> dict:
    paired = [r for r in rows if r["verdict"] != "incomplete"]
    faster = [r for r in paired if r["verdict"] == "rust_faster"]
    slower = [r for r in paired if r["verdict"] == "rust_slower"]
    tie = [r for r in paired if r["verdict"] == "wall_tie"]

    def sum_wall(key: str) -> float | None:
        vals = [r[key] for r in paired if r[key] is not None]
        return sum(vals) if vals else None

    j_sum, r_sum = sum_wall("java_wall_sec"), sum_wall("rust_wall_sec")
    return {
        "n_shards_seen": len(rows),
        "n_paired": len(paired),
        "n_rust_faster": len(faster),
        "n_rust_slower": len(slower),
        "n_wall_tie": len(tie),
        "sum_java_wall_sec": j_sum,
        "sum_rust_wall_sec": r_sum,
        "sum_wall_ratio_rust_over_java": (r_sum / j_sum)
        if (j_sum and r_sum and j_sum > 0)
        else None,
        "worst_rust_slower": sorted(
            slower,
            key=lambda r: (r["wall_ratio_rust_over_java"] or 0),
            reverse=True,
        )[:10],
        "best_rust_faster": sorted(
            faster,
            key=lambda r: (r["wall_ratio_rust_over_java"] or 99),
        )[:10],
    }


def to_markdown(run_id: str | None, rows: list[dict], summary: dict) -> str:
    lines = [
        f"# CI-subset shard perf (Java vs Rust)",
        "",
        f"Run: `{run_id or 'local'}`",
        "",
        f"- Paired shards: **{summary['n_paired']}** / {summary['n_shards_seen']}",
        f"- Rust faster (wall <0.9×): **{summary['n_rust_faster']}**",
        f"- Rust slower (wall >1.1×): **{summary['n_rust_slower']}**",
        f"- Wall tie (±10%): **{summary['n_wall_tie']}**",
    ]
    if summary.get("sum_java_wall_sec") is not None:
        lines.append(
            f"- Σ wall: Java **{fmt_sec(summary['sum_java_wall_sec'])}** · "
            f"Rust **{fmt_sec(summary['sum_rust_wall_sec'])}** · "
            f"ratio **{(summary['sum_wall_ratio_rust_over_java'] or 0):.2f}×**"
        )
    lines += [
        "",
        "| Shard | Java wall | Rust wall | Rust/Java | Java RSS | Rust RSS | Verdict |",
        "|-------|-----------|-----------|-----------|----------|----------|---------|",
    ]
    for r in rows:
        lines.append(
            "| {shard} | {jw} | {rw} | {ratio} | {jr} | {rr} | {v} |".format(
                shard=r["shard"],
                jw=fmt_sec(r["java_wall_sec"]),
                rw=fmt_sec(r["rust_wall_sec"]),
                ratio=(
                    f"{r['wall_ratio_rust_over_java']:.2f}×"
                    if r["wall_ratio_rust_over_java"] is not None
                    else "n/a"
                ),
                jr=fmt_mib(r["java_rss_kb"]),
                rr=fmt_mib(r["rust_rss_kb"]),
                v=r["verdict"],
            )
        )
    if summary["worst_rust_slower"]:
        lines += ["", "## Slowest Rust vs Java (top 10)", ""]
        for r in summary["worst_rust_slower"]:
            lines.append(
                f"- `{r['shard']}`: {fmt_sec(r['rust_wall_sec'])} vs "
                f"{fmt_sec(r['java_wall_sec'])} "
                f"({r['wall_ratio_rust_over_java']:.2f}×)"
            )
    if summary["best_rust_faster"]:
        lines += ["", "## Fastest Rust wins (top 10)", ""]
        for r in summary["best_rust_faster"]:
            lines.append(
                f"- `{r['shard']}`: {fmt_sec(r['rust_wall_sec'])} vs "
                f"{fmt_sec(r['java_wall_sec'])} "
                f"({r['wall_ratio_rust_over_java']:.2f}×)"
            )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", help="GitHub Actions run id")
    ap.add_argument("--artifact-dir", type=Path, help="Local unpacked artifacts")
    ap.add_argument("--json-out", type=Path, help="Write full JSON summary")
    ap.add_argument("--md-out", type=Path, help="Write Markdown table")
    ap.add_argument(
        "--keep-download",
        type=Path,
        help="Keep downloaded artifacts under this dir (default: temp)",
    )
    args = ap.parse_args()

    if not args.run_id and not args.artifact_dir:
        ap.error("need --run-id or --artifact-dir")

    cleanup = None
    if args.artifact_dir:
        root = args.artifact_dir
    else:
        root = args.keep_download or Path(tempfile.mkdtemp(prefix=f"giab-{args.run_id}-"))
        if not args.keep_download:
            cleanup = root
        print(f"[perf] downloading run {args.run_id} → {root}", file=sys.stderr)
        download_run_artifacts(args.run_id, root)

    try:
        by_shard = collect_from_dir(root)
        rows = build_rows(by_shard)
        summary = summarize(rows)
        payload = {
            "run_id": args.run_id,
            "summary": summary,
            "shards": rows,
        }
        if args.json_out:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(json.dumps(payload, indent=2) + "\n")
            print(f"[perf] wrote {args.json_out}", file=sys.stderr)
        md = to_markdown(args.run_id, rows, summary)
        if args.md_out:
            args.md_out.parent.mkdir(parents=True, exist_ok=True)
            args.md_out.write_text(md)
            print(f"[perf] wrote {args.md_out}", file=sys.stderr)
        else:
            print(md)
        print(json.dumps({"summary": summary}, indent=2), file=sys.stderr)
    finally:
        if cleanup and cleanup.exists():
            shutil.rmtree(cleanup, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

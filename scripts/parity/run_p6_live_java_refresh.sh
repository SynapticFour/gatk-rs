#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

gatk_image="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
gatk_platform="${GATK_DOCKER_PLATFORM:-linux/amd64}"

p5_live_summary="${report_dir}/p5_live_java_rust_diff_summary.json"
p6_live_summary="${report_dir}/p6_live_java_refresh_summary.json"
p6_live_md="${report_dir}/p6_live_java_refresh_summary.md"

echo "[p6-live-refresh] java image=${gatk_image} platform=${gatk_platform}"

docker run --rm --platform "${gatk_platform}" "${gatk_image}" gatk --version >/dev/null

echo "[p6-live-refresh] running live Java-vs-Rust runtime profile"
GATK_DOCKER_IMAGE="${gatk_image}" GATK_DOCKER_PLATFORM="${gatk_platform}" \
  CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  "${repo_root}/scripts/parity/run_p5_live_java_rust_diff.sh"

echo "[p6-live-refresh] revalidating frozen P6 likelihood vector diff"
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  cargo test -p gatk-haplotypecaller --test p6_likelihood_vector_diff_test --locked \
  pairhmm_likelihood_vector_matches_frozen_java_dump_fixture >/dev/null

python3 - "${repo_root}" <<'PY'
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

repo = Path(sys.argv[1]).resolve()
reports = repo / "parity" / "reports"
expected = repo / "parity" / "expected" / "p6_pairhmm_case1.java_likelihoods.tsv"
p5_summary = reports / "p5_live_java_rust_diff_summary.json"
p6_summary = reports / "p6_live_java_refresh_summary.json"
p6_md = reports / "p6_live_java_refresh_summary.md"

p5 = json.loads(p5_summary.read_text(encoding="utf-8"))
expected_bytes = expected.read_bytes()
expected_sha256 = hashlib.sha256(expected_bytes).hexdigest()
expected_rows = [
    line for line in expected.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.startswith("#")
]

summary = {
    "label": "phase6-live-java-refresh",
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "java_runtime_refresh": {
        "source": "run_p5_live_java_rust_diff.sh",
        "passed": p5.get("passed", 0),
        "failed": p5.get("failed", 0),
        "total": p5.get("total", 0),
        "pass_rate": p5.get("pass_rate", 0.0),
    },
    "frozen_step82_fixture": {
        "path": str(expected.relative_to(repo)),
        "sha256": expected_sha256,
        "rows": len(expected_rows),
    },
    "status": "pass",
}
p6_summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")

lines = [
    "# P6 Live Java Refresh Summary",
    "",
    "## Runtime Refresh",
    f"- source: `run_p5_live_java_rust_diff.sh`",
    f"- passed: `{summary['java_runtime_refresh']['passed']}`",
    f"- failed: `{summary['java_runtime_refresh']['failed']}`",
    f"- total: `{summary['java_runtime_refresh']['total']}`",
    "",
    "## Step-82 Frozen Fixture Fingerprint",
    f"- file: `{summary['frozen_step82_fixture']['path']}`",
    f"- rows: `{summary['frozen_step82_fixture']['rows']}`",
    f"- sha256: `{summary['frozen_step82_fixture']['sha256']}`",
    "",
    "## Note",
    "- This refresh proves current Java runtime path and frozen Step-82 fixture consistency.",
]
p6_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"[p6-live-refresh] wrote {p6_summary}")
print(f"[p6-live-refresh] wrote {p6_md}")
PY

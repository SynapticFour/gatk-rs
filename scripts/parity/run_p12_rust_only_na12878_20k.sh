#!/usr/bin/env bash
# Re-run only the Rust HaplotypeCaller for P12, reusing an existing Java VCF from a prior run.
# Speed: uses release binary when P12_CARGO_RELEASE=1 (default), higher parallelism via cargo/Rayon.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
data_dir="${P12_DATA_DIR:-${repo_root}/parity/realworld/na12878_20k_b37}"
mkdir -p "${report_dir}" "${data_dir}"
cd "${repo_root}"

java_vcf="${P12_JAVA_VCF:-${report_dir}/p12_realworld_na12878_20k.java.vcf}"
rust_vcf="${report_dir}/p12_realworld_na12878_20k.rust.vcf"
json_out="${report_dir}/p12_realworld_na12878_20k.json"
md_out="${report_dir}/p12_realworld_na12878_20k.md"
target_dir="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target-parity}"
interval="${P12_INTERVAL:-}"

# Reference required (same as full P12).
reference="${P12_REFERENCE:-}"
if [[ -z "${reference}" ]]; then
  echo "[p12-rust-only] ERROR: set P12_REFERENCE to the same FASTA used for the Java run." >&2
  exit 1
fi

if [[ ! -f "${java_vcf}" ]]; then
  echo "[p12-rust-only] ERROR: cached Java VCF missing: ${java_vcf}" >&2
  echo "  Run ./scripts/parity/run_p12_realworld_na12878_20k.sh once, or set P12_JAVA_VCF." >&2
  exit 1
fi

# Default: use all logical CPUs for Rayon + cargo build jobs unless overridden.
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"

P12_CARGO_RELEASE="${P12_CARGO_RELEASE:-1}"
cargo_run=(cargo run)
if [[ "${P12_CARGO_RELEASE}" == "1" ]]; then
  cargo_run+=(--release)
fi
cargo_run+=(--quiet --bin gatk-rs -- HaplotypeCaller -R "${reference}" -I "${data_dir}/NA12878_20k.b37.bam" -O "${rust_vcf}")
if [[ -n "${interval}" ]]; then
  cargo_run+=(-L "${interval}")
fi

mkdir -p "${target_dir}"
set +e
CARGO_TARGET_DIR="${target_dir}" "${cargo_run[@]}" >/dev/null 2>&1
rust_exit=$?
set -e
java_exit=0

notes_json="$(python3 - <<PY
import json
print(json.dumps({
  "mode": "rust_only_reuse_java_vcf",
  "java_vcf_reused": "${java_vcf}",
  "cargo_release": "${P12_CARGO_RELEASE}",
  "rayon_threads": "${RAYON_NUM_THREADS}",
}))
PY
)"

python3 "${repo_root}/scripts/parity/p12_na12878_summarize.py" \
  "${json_out}" "${md_out}" "${java_vcf}" "${rust_vcf}" "${java_exit}" "${rust_exit}" "${notes_json}"

#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
data_dir="${P12_DATA_DIR:-${repo_root}/parity/realworld/na12878_20k_b37}"
mkdir -p "${report_dir}" "${data_dir}"
cd "${repo_root}"

json_out="${report_dir}/p12_realworld_na12878_20k.json"
md_out="${report_dir}/p12_realworld_na12878_20k.md"

base_url="${P12_S3_BASE_URL:-https://gatk-test-data.s3.amazonaws.com/wgs_bam/NA12878_20k_b37}"
bam="${P12_BAM:-${data_dir}/NA12878_20k.b37.bam}"
bai="${P12_BAI:-${bam%.bam}.bai}"
if [[ "${bai}" == "${bam}" ]]; then
  bai="${bam}.bai"
fi

download_if_missing() {
  local url="$1"
  local out="$2"
  if [[ -f "${out}" ]]; then
    return 0
  fi
  curl -fsSL "${url}" -o "${out}"
}

# Default corpus download; skip when caller supplies P12_BAM (e.g. R3 dense GIAB window).
if [[ -z "${P12_BAM:-}" ]]; then
  download_if_missing "${base_url}/NA12878_20k.b37.bam" "${bam}"
  download_if_missing "${base_url}/NA12878_20k.b37.bai" "${bai}"
elif [[ ! -f "${bam}" ]]; then
  echo "[p12-realworld] P12_BAM missing: ${bam}" >&2
  exit 1
fi

reference="${P12_REFERENCE:-}"
if [[ -z "${reference}" ]]; then
  python3 - "${json_out}" "${md_out}" "${bam}" "${bai}" <<'PY'
import json
import pathlib
import sys

json_out = pathlib.Path(sys.argv[1])
md_out = pathlib.Path(sys.argv[2])
bam = pathlib.Path(sys.argv[3])
bai = pathlib.Path(sys.argv[4])
payload = {
    "label": "phase12-realworld-na12878-20k",
    "status": "data_ready_reference_missing",
    "dataset": {
        "bam": str(bam),
        "bai": str(bai),
    },
    "notes": "Set P12_REFERENCE to a b37-compatible FASTA (+ .fai/.dict) to run Java/Rust HC and strict output comparison.",
}
json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
md_out.write_text(
    "\n".join(
        [
            "# P12 Real-world NA12878 20k",
            "",
            "- status: **data_ready_reference_missing**",
            f"- bam: `{bam}`",
            f"- bai: `{bai}`",
            "- next: export `P12_REFERENCE=/path/to/human_g1k_v37.fasta` and rerun",
        ]
    )
    + "\n",
    encoding="utf-8",
)
print("[p12-realworld] dataset downloaded; reference missing, run skipped")
PY
  exit 0
fi

java_vcf="${P12_JAVA_VCF:-${report_dir}/p12_realworld_na12878_20k.java.vcf}"
rust_vcf="${P12_RUST_VCF:-${report_dir}/p12_realworld_na12878_20k.rust.vcf}"
json_out="${P12_JSON_OUT:-${json_out}}"
md_out="${P12_MD_OUT:-${md_out}}"
target_dir="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"
# NA12878_20k_b37 is sparse; chr20:413419-463418 has ~4 reads. Default to a dense chr2 window.
interval="${P12_INTERVAL:-2:92300000-92350000}"
mkdir -p "${target_dir}"

# M4 16GB defaults — do not fan out to all logical CPUs.
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
P12_CARGO_RELEASE="${P12_CARGO_RELEASE:-0}"

set +e
if [[ "${P12_SKIP_JAVA:-0}" == "1" && -s "${java_vcf}" ]]; then
  echo "[p12-realworld] skipping Java HC (P12_SKIP_JAVA=1, existing ${java_vcf})"
  java_exit=0
else
  echo "[p12-realworld] running Java HaplotypeCaller via Docker (often 5–15 min on arm64; no stdout)…"
  java_cmd=(
    docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}"
    -v "${repo_root}:${repo_root}"
    -w "${repo_root}"
    "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
    gatk HaplotypeCaller
    -R "${reference}"
    -I "${bam}"
    -O "${java_vcf}"
    --verbosity ERROR
  )
  if [[ -n "${interval}" ]]; then
    java_cmd+=(-L "${interval}")
  fi
  "${java_cmd[@]}"
  java_exit=$?
fi

profile_dir="$([[ "${P12_CARGO_RELEASE}" == "1" ]] && echo release || echo debug)"
rust_bin="${target_dir}/${profile_dir}/gatk-rs"
use_prebuilt="${P12_USE_PREBUILT:-}"
if [[ -z "${use_prebuilt}" && -x "${rust_bin}" ]]; then
  use_prebuilt=1
fi
if [[ "${use_prebuilt}" == "1" && -x "${rust_bin}" ]]; then
  echo "[p12-realworld] Rust HC via prebuilt ${rust_bin} (~2–3 min; rebuild: ./scripts/parity/build_gatk_rs_release.sh)"
  "${rust_bin}" HaplotypeCaller -R "${reference}" -I "${bam}" -O "${rust_vcf}" ${interval:+-L "${interval}"}
  rust_exit=$?
elif [[ "${P12_FORCE_CARGO_BUILD:-0}" == "1" ]]; then
  echo "[p12-realworld] building via cargo (visible progress; avoid cargo -q)…"
  build_args=(-p gatk-cli --bin gatk-rs -j "${CARGO_BUILD_JOBS}")
  if [[ "${P12_CARGO_RELEASE}" == "1" ]]; then
    build_args+=(--release)
  fi
  CARGO_TARGET_DIR="${target_dir}" cargo build "${build_args[@]}"
  echo "[p12-realworld] running ${rust_bin}…"
  "${rust_bin}" HaplotypeCaller -R "${reference}" -I "${bam}" -O "${rust_vcf}" ${interval:+-L "${interval}"}
  rust_exit=$?
else
  echo "[p12-realworld] no prebuilt ${rust_bin}; run ./scripts/parity/build_gatk_rs_release.sh or P12_FORCE_CARGO_BUILD=1" >&2
  rust_exit=127
fi
set -e

notes_json="$(python3 - <<PY
import json
print(json.dumps({
    "mode": "java_and_rust",
    "cargo_release": "${P12_CARGO_RELEASE}",
    "rayon_threads": "${RAYON_NUM_THREADS}",
}))
PY
)"
python3 "${repo_root}/scripts/parity/p12_na12878_summarize.py" \
  "${json_out}" "${md_out}" "${java_vcf}" "${rust_vcf}" "${java_exit}" "${rust_exit}" "${notes_json}"

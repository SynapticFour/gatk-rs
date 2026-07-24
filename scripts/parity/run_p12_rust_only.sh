#!/usr/bin/env bash
# P12 Rust-only path: no cargo, no docker — uses prebuilt target/release/gatk-rs only.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

reference="${P12_REFERENCE:?Set P12_REFERENCE to a b37 FASTA}"
interval="${P12_INTERVAL:-2:92300000-92350000}"
bam="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
rust_bin="${target_dir}/release/gatk-rs"
rust_vcf="${P12_RUST_VCF:-${repo_root}/parity/reports/p12_realworld_na12878_20k.rust.vcf}"
java_vcf="${P12_JAVA_VCF:-${repo_root}/parity/reports/p12_realworld_na12878_20k.java.vcf}"
json_out="${repo_root}/parity/reports/p12_realworld_na12878_20k.json"
md_out="${repo_root}/parity/reports/p12_realworld_na12878_20k.md"

if [[ ! -x "${rust_bin}" ]]; then
  echo "[p12-rust-only] missing ${rust_bin}" >&2
  echo "[p12-rust-only] run: ./scripts/parity/build_gatk_rs_release.sh" >&2
  exit 1
fi

mkdir -p "$(dirname "${rust_vcf}")"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-4}"
echo "[p12-rust-only] binary=${rust_bin}"
echo "[p12-rust-only] interval=${interval} threads=${RAYON_NUM_THREADS}"
echo "[p12-rust-only] read_augment=removed (Java-only call_region path)"
echo "[p12-rust-only] running HaplotypeCaller (~2–3 min)…"
"${rust_bin}" HaplotypeCaller \
  -R "${reference}" -I "${bam}" -O "${rust_vcf}" -L "${interval}"
rust_exit=$?
echo "[p12-rust-only] HC exit=${rust_exit} wrote ${rust_vcf}"

if [[ ! -s "${java_vcf}" ]]; then
  echo "[p12-rust-only] warn: no Java VCF at ${java_vcf}; summary will be rust-only" >&2
  java_exit=0
else
  java_exit=0
fi

notes_json='{"mode":"rust_only_prebuilt","cargo_release":"1","read_augment":"removed"}'
diff_dir="${repo_root}/parity/reports/p12_diff"
python3 "${repo_root}/scripts/parity/p12_na12878_summarize.py" \
  "${json_out}" "${md_out}" "${java_vcf}" "${rust_vcf}" "${java_exit}" "${rust_exit}" "${notes_json}" \
  "${diff_dir}"
echo "[p12-rust-only] reports: ${json_out} ${md_out}"

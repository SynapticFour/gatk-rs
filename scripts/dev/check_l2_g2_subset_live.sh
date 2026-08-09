#!/usr/bin/env bash
# Fast L2 slice: g2-subset-live dumps vs frozen Java oracles (no Docker / javac).
# Catches the Peak min-kmer RT merge regression class before CI.1 / Parity Smoke.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

CASES="${ROOT}/parity/fixtures/hc-full-parity/g2-subset-live/cases.tsv"
COMPARE="${ROOT}/scripts/parity/compare_hc_full_parity_l2.py"
JAVA_DIR="${ROOT}/parity/fixtures/hc-full-parity/java_dumps/g2-subset-live"
BAM_DIR="${ROOT}/parity/build/sam-indexed-bam"

if [[ ! -f "${CASES}" || ! -f "${COMPARE}" ]]; then
  echo "[g2-subset-live] missing fixtures/compare — skip" >&2
  exit 0
fi

if [[ ! -d "${JAVA_DIR}" ]]; then
  echo "[g2-subset-live] missing java_dumps — skip (CI will enforce)" >&2
  exit 0
fi

export OPENSSL_DIR="${OPENSSL_DIR:-/opt/homebrew/opt/openssl@3}"
unset CARGO_TARGET_DIR HTTP_PROXY HTTPS_PROXY ALL_PROXY || true

echo "[g2-subset-live] building hc_full_parity_gate_dump (release, parity_harness)"
cargo build -p gatk-haplotypecaller --features parity_harness \
  --example hc_full_parity_gate_dump --release -q
BIN="${ROOT}/target/release/examples/hc_full_parity_gate_dump"

ok=0
fail=0
while IFS=$'\t' read -r id ref sam interval _pad target maxalleles profile _expected; do
  [[ "${id}" == case_id* || "${id}" == \#* || -z "${id}" ]] && continue
  base="$(basename "${sam}" .sam)"
  bam="${BAM_DIR}/${base}.bam"
  if [[ ! -f "${bam}" ]]; then
    echo "[g2-subset-live] SKIP ${id}: missing ${bam} (stage BAMs via scripts/ci/stage_indexed_fixture_bams.sh)" >&2
    continue
  fi
  if [[ ! -f "${ref}" ]]; then
    echo "[g2-subset-live] SKIP ${id}: missing ${ref}" >&2
    continue
  fi
  rust_out="$(mktemp)"
  if ! "${BIN}" assembly-region-genotype-subset \
    "${ref}" "${bam}" "${interval}" "${target}" "${profile}" "${maxalleles}" \
    >"${rust_out}" 2>/dev/null; then
    echo "[g2-subset-live] FAIL ${id}: dump exited non-zero" >&2
    rm -f "${rust_out}"
    fail=$((fail + 1))
    continue
  fi
  java="$(ls "${JAVA_DIR}/${id}_"*.tsv 2>/dev/null | head -1 || true)"
  if [[ -z "${java}" ]]; then
    echo "[g2-subset-live] FAIL ${id}: no frozen java dump" >&2
    rm -f "${rust_out}"
    fail=$((fail + 1))
    continue
  fi
  if python3 "${COMPARE}" "${rust_out}" "${java}" >/dev/null 2>&1; then
    echo "[g2-subset-live] PASS ${id}"
    ok=$((ok + 1))
  else
    echo "[g2-subset-live] FAIL ${id} (rust vs ${java})" >&2
    fail=$((fail + 1))
  fi
  rm -f "${rust_out}"
done <"${CASES}"

if [[ "${fail}" -gt 0 ]]; then
  echo "[g2-subset-live] FAIL ok=${ok} fail=${fail}" >&2
  exit 1
fi
if [[ "${ok}" -eq 0 ]]; then
  echo "[g2-subset-live] SKIP: no cases ran (stage BAMs first)" >&2
  exit 0
fi
echo "[g2-subset-live] PASS (${ok} cases)"
exit 0

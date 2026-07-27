#!/usr/bin/env bash
# Copy tracked seed triage JSONL into parity/reports/ when missing.
# Reports are gitignored, so CI checkouts otherwise fail phase*-mismatch-triage-check.
set -euo pipefail

phase="${1:?usage: ensure_mismatch_triage.sh <p5|p6|p7|p8|p9|p11>}"
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
seed="${repo_root}/parity/fixtures/mismatch-triage/${phase}_mismatch_triage.jsonl"
dest_dir="${repo_root}/parity/reports"
dest="${dest_dir}/${phase}_mismatch_triage.jsonl"

if [[ ! -f "${seed}" ]]; then
  echo "[ensure-triage] missing seed ${seed}" >&2
  exit 2
fi
mkdir -p "${dest_dir}"
if [[ ! -f "${dest}" ]]; then
  cp "${seed}" "${dest}"
  echo "[ensure-triage] seeded ${dest} from fixtures"
fi

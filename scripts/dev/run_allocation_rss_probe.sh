#!/usr/bin/env bash
# Peak-RSS before/after for the activity-scoring allocation audit fix.
# "before" = clone-strata-every-locus mode; "after" = production borrow mode.
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"
out="${root}/scripts/dev/allocation_audit_rss.txt"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"

run_mode() {
  local mode="$1"
  local label="$2"
  echo "=== ${label} (mode=${mode}) ==="
  if /usr/bin/time -l true >/dev/null 2>&1; then
    /usr/bin/time -l cargo run -p gatk-haplotypecaller --release --example allocation_rss_probe -- "${mode}" \
      >"/tmp/alloc_${mode}.out" 2>"/tmp/alloc_${mode}.time"
  else
    /usr/bin/time -v cargo run -p gatk-haplotypecaller --release --example allocation_rss_probe -- "${mode}" \
      >"/tmp/alloc_${mode}.out" 2>"/tmp/alloc_${mode}.time"
  fi
  cat "/tmp/alloc_${mode}.out"
  # Extract peak RSS line (macOS or GNU time).
  rg -n 'maximum resident set size|Maximum resident set size|peak memory footprint' "/tmp/alloc_${mode}.time" || true
  echo
}

{
  echo "Allocation audit Peak-RSS probe — $(date -u +%Y-%m-%dT%H:%MZ)"
  echo "Fixture: synthetic 3×80 pileups × 50000 loci (activity multi-sample path)"
  echo
  run_mode clone "BEFORE (clone nonempty strata each locus)"
  run_mode borrow "AFTER (borrow strata / production API)"
  echo "Delta note: compare maximum resident set size (macOS bytes) or Maximum resident set size (Linux kbytes)."
} | tee "${out}"
echo "Wrote ${out}"

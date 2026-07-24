#!/usr/bin/env bash
# P12: build release, run with read augment ON then OFF, then cluster diagnose.
# Usage:
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_p12_both_modes.sh
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

if [[ -z "${P12_REFERENCE:-}" ]]; then
  echo "Set P12_REFERENCE to a b37 FASTA, e.g.:" >&2
  echo '  export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"' >&2
  exit 1
fi

./scripts/parity/build_gatk_rs_release.sh

echo ""
echo "=== P12 run 1/2: read augment ENABLED ==="
unset GATK_RS_HC_DISABLE_READ_AUGMENT
export P12_RUST_VCF="${repo_root}/parity/reports/p12_realworld_na12878_20k.rust.vcf"
./scripts/parity/run_p12_rust_only.sh

echo ""
echo "=== P12 run 2/2: read augment DISABLED ==="
export GATK_RS_HC_DISABLE_READ_AUGMENT=1
export P12_RUST_VCF="${repo_root}/parity/reports/p12_realworld_na12878_20k.rust.no_augment.vcf"
./scripts/parity/run_p12_rust_only.sh

unset GATK_RS_HC_DISABLE_READ_AUGMENT

echo ""
echo "=== cluster diagnose ==="
./scripts/parity/diagnose_p12_cluster_assembly.sh

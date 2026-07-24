#!/usr/bin/env bash
# VariantFiltration boundary parity: Java GATK 4.4 vs gatk-rs on synthetic sites
# that sit exactly on / just across official hard-filter thresholds.
#
# Usage:
#   ./scripts/parity/run_variant_filtration_parity.sh
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
# shellcheck source=lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"

fixture="${repo_root}/parity/variant_filtration/boundary_sites.vcf"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="${report_dir}/variant_filtration_${stamp}.log"
java_out="${report_dir}/variant_filtration.java.vcf"
rust_out="${report_dir}/variant_filtration.rust.vcf"
exec > >(tee -a "${log}") 2>&1

echo "=== VariantFiltration boundary parity ${stamp} ==="

target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
rust_bin="${target_dir}/debug/gatk-rs"
[[ -x "${rust_bin}" ]] || rust_bin="${target_dir}/release/gatk-rs"
if [[ ! -x "${rust_bin}" ]]; then
  export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
  cargo build -p gatk-cli --bin gatk-rs
  rust_bin="${target_dir}/debug/gatk-rs"
fi

# Index for Java when using intervals is not required for whole-file filtration,
# but IndexFeatureFile keeps GATK happier on some builds.
"${repo_root}/scripts/parity/run_java_gatk.sh" \
  "${report_dir}/vf.index.stdout" \
  IndexFeatureFile -I "${fixture}" || true

set +e
echo "[vf-parity] Java VariantFiltration (SNP hard-filter table)…"
"${repo_root}/scripts/parity/run_java_gatk.sh" \
  "${report_dir}/vf.java.stdout" \
  VariantFiltration \
  -V "${fixture}" \
  -O "${java_out}" \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
java_exit=$?

echo "[vf-parity] Rust VariantFiltration…"
"${rust_bin}" variant-filtration \
  -V "${fixture}" \
  -O "${rust_out}" \
  --preset snp
rust_exit=$?
set -e

echo "[vf-parity] java_exit=${java_exit} rust_exit=${rust_exit}"
if [[ "${java_exit}" -ne 0 || "${rust_exit}" -ne 0 ]]; then
  exit 1
fi

python3 - <<'PY'
from pathlib import Path
import sys

def body_filters(path: Path):
    out = {}
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        key = (f[0], int(f[1]))
        filt = f[6]
        # Normalize filter token order for comparison
        tokens = tuple(sorted(t for t in filt.split(";") if t and t != "."))
        out[key] = tokens
    return out

ja = body_filters(Path("parity/reports/variant_filtration.java.vcf"))
rb = body_filters(Path("parity/reports/variant_filtration.rust.vcf"))
keys = sorted(set(ja) | set(rb))
bad = 0
for k in keys:
    a, b = ja.get(k), rb.get(k)
    if a != b:
        print(f"[vf-parity] FILTER mismatch {k}: java={a} rust={b}")
        bad += 1
if bad:
    print(f"[vf-parity] FAIL mismatches={bad}")
    sys.exit(1)
print(f"[vf-parity] OK sites={len(keys)} identical FILTER decisions")
PY

echo "[vf-parity] log=${log}"

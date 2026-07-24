#!/usr/bin/env bash
# PR-gate: hap.py equivalence smoke on the smallest tracked region
# (parity/fixtures, chr1:1-32) via Docker for GATK 4.4 + hap.py.
#
# Observable contract:
#   Run Java GATK4 HC and gatk-rs HC on identical tiny inputs, score both
#   with Illumina hap.py against the same truth+BED, gate on |ΔF1| ≤ threshold
#   (default 0.02) via gatk-rs-equiv.
#
# Env overrides:
#   RUST_BINARY, EQUIV_BINARY, HAPPY_DOCKER_IMAGE, GATK_DOCKER_IMAGE,
#   PR_HAPPY_OUT, PR_HAPPY_F1_DELTA, PR_HAPPY_SKIP_PULL
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
fixtures="${repo_root}/parity/fixtures"
ref="${fixtures}/reference.fa"
bam="${fixtures}/sample.bam"
truth="${fixtures}/pr-happy/truth.vcf"
bed="${fixtures}/pr-happy/confident.bed"
interval="chr1:1-32"

out="${PR_HAPPY_OUT:-${repo_root}/parity/reports/pr-happy-smoke}"
f1_delta="${PR_HAPPY_F1_DELTA:-0.02}"
happy_image="${HAPPY_DOCKER_IMAGE:-jmcdani20/hap.py:v0.3.12}"
# Entrypoint inside jmcdani20/hap.py images (override if using another tag).
happy_entry="${HAPPY_DOCKER_ENTRYPOINT:-/opt/hap.py/bin/hap.py}"

rust_bin="${RUST_BINARY:-}"
equiv_bin="${EQUIV_BINARY:-}"
if [[ -z "${rust_bin}" ]]; then
  for cand in \
    "${repo_root}/target/debug/gatk-rs" \
    "${repo_root}/target/release/gatk-rs"; do
    if [[ -x "${cand}" ]]; then
      rust_bin="${cand}"
      break
    fi
  done
fi
if [[ -z "${equiv_bin}" ]]; then
  for cand in \
    "${repo_root}/target/debug/gatk-rs-equiv" \
    "${repo_root}/target/release/gatk-rs-equiv"; do
    if [[ -x "${cand}" ]]; then
      equiv_bin="${cand}"
      break
    fi
  done
fi

if [[ -z "${rust_bin}" || ! -x "${rust_bin}" ]]; then
  echo "[pr-happy] missing gatk-rs binary (build workspace or set RUST_BINARY)" >&2
  exit 2
fi
if [[ -z "${equiv_bin}" || ! -x "${equiv_bin}" ]]; then
  echo "[pr-happy] missing gatk-rs-equiv binary (build workspace or set EQUIV_BINARY)" >&2
  exit 2
fi

for f in "${ref}" "${bam}" "${truth}" "${bed}" "${ref}.fai"; do
  if [[ ! -e "${f}" ]]; then
    echo "[pr-happy] missing fixture: ${f}" >&2
    exit 2
  fi
done

if ! command -v docker >/dev/null 2>&1; then
  echo "[pr-happy] docker is required for GATK + hap.py on the PR runner" >&2
  exit 2
fi

mkdir -p "${out}/bin"
gatk_wrap="${out}/bin/gatk-docker"
happy_wrap="${out}/bin/hap.py-docker"

cat >"${gatk_wrap}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec docker run --rm \\
  --platform "${GATK_DOCKER_PLATFORM}" \\
  -v "${repo_root}:${repo_root}" \\
  -w "${repo_root}" \\
  "${GATK_DOCKER_IMAGE}" \\
  gatk "\$@"
EOF
chmod +x "${gatk_wrap}"

cat >"${happy_wrap}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec docker run --rm \\
  --platform linux/amd64 \\
  -v "${repo_root}:${repo_root}" \\
  -w "${repo_root}" \\
  "${happy_image}" \\
  "${happy_entry}" "\$@"
EOF
chmod +x "${happy_wrap}"

if [[ "${PR_HAPPY_SKIP_PULL:-0}" != "1" ]]; then
  echo "[pr-happy] ensuring Docker images…"
  docker pull --platform "${GATK_DOCKER_PLATFORM}" "${GATK_DOCKER_IMAGE}"
  docker pull --platform linux/amd64 "${happy_image}"
fi

echo "[pr-happy] GATK image:  ${GATK_DOCKER_IMAGE}"
echo "[pr-happy] hap.py image: ${happy_image}"
echo "[pr-happy] rust:         ${rust_bin}"
echo "[pr-happy] equiv:        ${equiv_bin}"
echo "[pr-happy] interval:     ${interval}"
echo "[pr-happy] out:          ${out}"

rm -rf "${out}/equiv"
mkdir -p "${out}/equiv"

set +e
"${equiv_bin}" run \
  --java-gatk-bin "${gatk_wrap}" \
  --rust-binary "${rust_bin}" \
  --reference "${ref}" \
  --bam "${bam}" \
  --truth-vcf "${truth}" \
  --confident-regions "${bed}" \
  --interval "${interval}" \
  --out "${out}/equiv" \
  --engine happy \
  --happy-bin "${happy_wrap}" \
  --f1-delta-threshold "${f1_delta}" \
  --threads 1 \
  --skip-disk-check
rc=$?
set -e

if [[ -f "${out}/equiv/REPORT.md" ]]; then
  echo "[pr-happy] --- REPORT.md (tail) ---"
  tail -n 40 "${out}/equiv/REPORT.md" || true
fi

if [[ "${rc}" -ne 0 ]]; then
  echo "[pr-happy] FAIL: gatk-rs-equiv exited ${rc}" >&2
  exit "${rc}"
fi

echo "[pr-happy] PASS (|ΔF1| ≤ ${f1_delta} on ${interval})"

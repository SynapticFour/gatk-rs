#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
image="${CI_LINUX_IMAGE:-rust:slim}"
platform="${CI_LINUX_PLATFORM:-linux/amd64}"
run_parity="${RUN_PARITY_SMOKE:-0}"

# Keep container builds laptop-survivable when this helper is used from an M4 host.
cmd='set -euo pipefail; export PATH="/usr/local/cargo/bin:/usr/local/rustup/bin:$PATH"; export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"; export CARGO_PROFILE_RELEASE_LTO="${CARGO_PROFILE_RELEASE_LTO:-false}"; export CARGO_PROFILE_RELEASE_CODEGEN_UNITS="${CARGO_PROFILE_RELEASE_CODEGEN_UNITS:-16}"; export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"; export GATK_DOCKER_IMAGE="${GATK_DOCKER_IMAGE:-broadinstitute/gatk:4.4.0.0}"; export GATK_DOCKER_PLATFORM="${GATK_DOCKER_PLATFORM:-linux/amd64}"; apt-get update; apt-get install -y pkg-config build-essential libssl-dev zlib1g-dev libbz2-dev liblzma-dev curl samtools python3 docker.io; if command -v rustup >/dev/null 2>&1; then rustup component add rustfmt clippy; fi; cargo fmt --all -- --check; cargo clippy --locked --workspace --all-targets --all-features -- -D warnings; cargo build --locked --workspace --verbose -j "${CARGO_BUILD_JOBS}"; cargo test --locked --workspace --verbose -j "${CARGO_BUILD_JOBS}"; cargo test --locked -p gatk-core --release --lib -j "${CARGO_BUILD_JOBS}"; ./scripts/parity/run_foundation_gate.sh'

if [[ "${run_parity}" == "1" ]]; then
  cmd+="; docker --version >/dev/null 2>&1 || (apt-get update && apt-get install -y docker.io); docker pull us.gcr.io/broad-gatk/gatk:4.4.0.0; GATK_DOCKER_IMAGE=us.gcr.io/broad-gatk/gatk:4.4.0.0 PARITY_REQUIRE_SAMTOOLS=1 ./scripts/parity/run_parity_smoke.sh"
fi

docker_args=(
  run --rm -t
  --platform "${platform}"
  -v "${repo_root}:/work"
  -w /work
)

if [[ -S /var/run/docker.sock ]]; then
  docker_args+=( -v /var/run/docker.sock:/var/run/docker.sock )
fi
docker_args+=( -e "PARITY_HOST_REPO_ROOT=${repo_root}" )
if [[ -n "${GATK_DOCKER_IMAGE:-}" ]]; then
  docker_args+=( -e "GATK_DOCKER_IMAGE=${GATK_DOCKER_IMAGE}" )
fi
if [[ -n "${GATK_DOCKER_PLATFORM:-}" ]]; then
  docker_args+=( -e "GATK_DOCKER_PLATFORM=${GATK_DOCKER_PLATFORM}" )
fi

echo "Running Linux CI parity in ${image} (${platform})"
docker "${docker_args[@]}" "${image}" bash -lc "${cmd}"

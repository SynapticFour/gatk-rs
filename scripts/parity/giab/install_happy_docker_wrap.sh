#!/usr/bin/env bash
# Install a Docker-backed `hap.py` wrapper for gatk-rs-equiv (preferred engine).
#
# Env:
#   HAPPY_DOCKER_IMAGE      (default jmcdani20/hap.py:v0.3.12)
#   HAPPY_DOCKER_ENTRYPOINT (default /opt/hap.py/bin/hap.py)
#   HAPPY_WRAP_DIR          (default $PWD/.giab-happy-bin)
#   HAPPY_SKIP_PULL         (1 = skip docker pull)
#   GIAB_REPO_ROOT          bind-mount root (default: repo root)
#
# Prints the wrapper path on stdout (last line); logs to stderr.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
repo_root="${GIAB_REPO_ROOT:-${repo_root}}"

happy_image="${HAPPY_DOCKER_IMAGE:-jmcdani20/hap.py:v0.3.12}"
happy_entry="${HAPPY_DOCKER_ENTRYPOINT:-/opt/hap.py/bin/hap.py}"
wrap_dir="${HAPPY_WRAP_DIR:-${repo_root}/.giab-happy-bin}"

if ! command -v docker >/dev/null 2>&1; then
  echo "[giab-happy] docker required to install hap.py wrapper" >&2
  exit 2
fi

mkdir -p "${wrap_dir}"
happy_wrap="${wrap_dir}/hap.py"

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

if [[ "${HAPPY_SKIP_PULL:-0}" != "1" ]]; then
  echo "[giab-happy] pulling ${happy_image}…" >&2
  docker pull --platform linux/amd64 "${happy_image}"
fi

echo "[giab-happy] wrapper: ${happy_wrap}" >&2
echo "${happy_wrap}"

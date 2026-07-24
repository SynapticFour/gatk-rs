#!/usr/bin/env bash
set -euo pipefail

_script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${_script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
cd "${repo_root}"

checks_json="${repo_root}/parity/checks.json"
if [[ ! -f "${checks_json}" ]]; then
  echo "Missing checks config: ${checks_json}" >&2
  exit 2
fi

run_required() {
  python3 - <<'PY2'
import json, pathlib, subprocess, sys
cfg = json.loads(pathlib.Path("parity/checks.json").read_text(encoding="utf-8"))
required = cfg.get("required", [])
if not required:
    print("No required foundation checks configured.")
    raise SystemExit(2)
for item in required:
    cid = item.get("id", "unknown")
    cmd = item.get("command")
    desc = item.get("description", "")
    klass = item.get("acceptance_class", "unspecified")
    timeout_s = int(item.get("timeout_s", 0) or 0)
    owner = item.get("owner", "unassigned")
    print(
        f"[foundation|required] {cid} "
        f"(class={klass}, owner={owner}, timeout={timeout_s or 'none'}s): {desc}"
    )
    if not cmd:
        print(f"Missing command for required check {cid}", file=sys.stderr)
        raise SystemExit(2)
    try:
        res = subprocess.run(cmd, shell=True, timeout=timeout_s if timeout_s > 0 else None)
    except subprocess.TimeoutExpired:
        print(f"Required check timed out: {cid} (timeout={timeout_s}s)", file=sys.stderr)
        raise SystemExit(124)
    if res.returncode != 0:
        print(f"Required check failed: {cid} (exit={res.returncode})", file=sys.stderr)
        raise SystemExit(res.returncode)
print("All required foundation checks passed.")
PY2
}

run_advisory() {
  python3 - <<'PY2'
import json, pathlib, subprocess
cfg = json.loads(pathlib.Path("parity/checks.json").read_text(encoding="utf-8"))
advisory = cfg.get("advisory", [])
if not advisory:
    raise SystemExit(0)
for item in advisory:
    cid = item.get("id", "unknown")
    cmd = item.get("command")
    desc = item.get("description", "")
    klass = item.get("acceptance_class", "unspecified")
    timeout_s = int(item.get("timeout_s", 0) or 0)
    owner = item.get("owner", "unassigned")
    print(
        f"[foundation|advisory] {cid} "
        f"(class={klass}, owner={owner}, timeout={timeout_s or 'none'}s): {desc}"
    )
    if not cmd:
        print("  -> skipped (no command configured)")
        continue
    try:
        res = subprocess.run(cmd, shell=True, timeout=timeout_s if timeout_s > 0 else None)
        if res.returncode != 0:
            print(f"  -> advisory check failed (non-gating): {cid} exit={res.returncode}")
    except subprocess.TimeoutExpired:
        print(f"  -> advisory check timed out (non-gating): {cid} timeout={timeout_s}s")
PY2
}

run_required
if [[ "${FOUNDATION_RUN_ADVISORY:-0}" == "1" ]]; then
  run_advisory
fi

echo "Foundation gate completed successfully."

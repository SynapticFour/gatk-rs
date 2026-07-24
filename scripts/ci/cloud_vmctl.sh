#!/usr/bin/env bash
# Start / stop / status for a paid self-hosted runner VM.
#
# Used by both:
#   - genomewide-validation.yml  (label gatk-rs-genomewide; HCLOUD_SERVER_ID / AWS_INSTANCE_ID)
#   - benchmark.yml              (label gatk-rs-benchmark; map AWS_PERF_INSTANCE_ID /
#                                 HCLOUD_PERF_SERVER_ID into AWS_INSTANCE_ID / HCLOUD_SERVER_ID)
#
# Supports:
#   CLOUD_PROVIDER=hetzner  (HCLOUD_TOKEN + HCLOUD_SERVER_ID)
#   CLOUD_PROVIDER=aws      (AWS_* + AWS_INSTANCE_ID; needs aws CLI)
#
# Important: stop = power off (preserve disk). Never terminates/deletes the VM.
# Never reuse the same instance ID for correctness and perf roles.
#
# Usage:
#   ./scripts/ci/cloud_vmctl.sh start|stop|status
#
# See docs/ci/SELF_HOSTED_RUNNER_SETUP.md and docs/ci/PERF_BENCHMARK_HOST.md
set -euo pipefail

cmd="${1:-}"
if [[ -z "${cmd}" || "${cmd}" == "-h" || "${cmd}" == "--help" ]]; then
  cat <<'EOF'
Usage: cloud_vmctl.sh start|stop|status

Environment:
  CLOUD_PROVIDER=hetzner|aws
  Hetzner: HCLOUD_TOKEN, HCLOUD_SERVER_ID
  AWS:     AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_REGION, AWS_INSTANCE_ID
  Optional: CLOUD_VM_BOOT_GRACE_SEC, GH_RUNNER_WAIT_TOKEN, GH_REPO, RUNNER_LABEL
EOF
  exit 0
fi

provider="${CLOUD_PROVIDER:-}"
if [[ -z "${provider}" ]]; then
  echo "[cloud_vmctl] ERROR: set CLOUD_PROVIDER=hetzner|aws" >&2
  exit 2
fi

log() { echo "[cloud_vmctl] $*"; }

# -----------------------------------------------------------------------------
# Hetzner
# -----------------------------------------------------------------------------
hcloud_api() {
  local method="$1" path="$2"
  curl -fsSL -X "${method}" \
    -H "Authorization: Bearer ${HCLOUD_TOKEN}" \
    -H "Content-Type: application/json" \
    "https://api.hetzner.cloud/v1${path}"
}

hetzner_require() {
  [[ -n "${HCLOUD_TOKEN:-}" ]] || { echo "[cloud_vmctl] missing HCLOUD_TOKEN" >&2; exit 2; }
  [[ -n "${HCLOUD_SERVER_ID:-}" ]] || { echo "[cloud_vmctl] missing HCLOUD_SERVER_ID" >&2; exit 2; }
}

hetzner_status() {
  hetzner_require
  local json status
  json="$(hcloud_api GET "/servers/${HCLOUD_SERVER_ID}")"
  status="$(printf '%s' "${json}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["server"]["status"])')"
  log "hetzner server ${HCLOUD_SERVER_ID} status=${status}"
  printf '%s\n' "${status}"
}

hetzner_start() {
  hetzner_require
  local status
  status="$(hetzner_status)"
  if [[ "${status}" == "running" ]]; then
    log "already running"
    return 0
  fi
  log "poweron server ${HCLOUD_SERVER_ID}"
  hcloud_api POST "/servers/${HCLOUD_SERVER_ID}/actions/poweron" >/dev/null
  # Wait until running
  for _ in $(seq 1 60); do
    status="$(hetzner_status || true)"
    [[ "${status}" == "running" ]] && break
    sleep 5
  done
  if [[ "${status}" != "running" ]]; then
    echo "[cloud_vmctl] ERROR: server did not reach running (status=${status})" >&2
    exit 1
  fi
}

hetzner_stop() {
  hetzner_require
  local status
  status="$(hetzner_status)"
  if [[ "${status}" == "off" ]]; then
    log "already off"
    return 0
  fi
  # Prefer ACPI shutdown; fall back to poweroff if needed.
  log "shutdown server ${HCLOUD_SERVER_ID} (preserve disks)"
  if ! hcloud_api POST "/servers/${HCLOUD_SERVER_ID}/actions/shutdown" >/dev/null; then
    log "shutdown failed; trying poweroff"
    hcloud_api POST "/servers/${HCLOUD_SERVER_ID}/actions/poweroff" >/dev/null
  fi
  for _ in $(seq 1 60); do
    status="$(hetzner_status || true)"
    [[ "${status}" == "off" ]] && break
    sleep 5
  done
  if [[ "${status}" != "off" ]]; then
    log "WARNING: status still ${status}; forcing poweroff"
    hcloud_api POST "/servers/${HCLOUD_SERVER_ID}/actions/poweroff" >/dev/null || true
  fi
  log "stopped (disks retained — storage still bills)"
}

# -----------------------------------------------------------------------------
# AWS EC2
# -----------------------------------------------------------------------------
aws_require() {
  [[ -n "${AWS_INSTANCE_ID:-}" ]] || { echo "[cloud_vmctl] missing AWS_INSTANCE_ID" >&2; exit 2; }
  [[ -n "${AWS_REGION:-}" ]] || { echo "[cloud_vmctl] missing AWS_REGION" >&2; exit 2; }
  if ! command -v aws >/dev/null 2>&1; then
    echo "[cloud_vmctl] aws CLI required for CLOUD_PROVIDER=aws" >&2
    exit 2
  fi
}

aws_status() {
  aws_require
  local state
  state="$(aws ec2 describe-instances \
    --region "${AWS_REGION}" \
    --instance-ids "${AWS_INSTANCE_ID}" \
    --query 'Reservations[0].Instances[0].State.Name' \
    --output text)"
  log "aws instance ${AWS_INSTANCE_ID} state=${state}"
  printf '%s\n' "${state}"
}

aws_start() {
  aws_require
  local state
  state="$(aws_status)"
  if [[ "${state}" == "running" ]]; then
    log "already running"
    return 0
  fi
  log "start-instances ${AWS_INSTANCE_ID}"
  aws ec2 start-instances --region "${AWS_REGION}" --instance-ids "${AWS_INSTANCE_ID}" >/dev/null
  aws ec2 wait instance-running --region "${AWS_REGION}" --instance-ids "${AWS_INSTANCE_ID}"
  log "running"
}

aws_stop() {
  aws_require
  local state
  state="$(aws_status)"
  if [[ "${state}" == "stopped" ]]; then
    log "already stopped"
    return 0
  fi
  log "stop-instances ${AWS_INSTANCE_ID} (preserve EBS)"
  aws ec2 stop-instances --region "${AWS_REGION}" --instance-ids "${AWS_INSTANCE_ID}" >/dev/null
  aws ec2 wait instance-stopped --region "${AWS_REGION}" --instance-ids "${AWS_INSTANCE_ID}"
  log "stopped (EBS retained — storage still bills)"
}

# -----------------------------------------------------------------------------
# Optional: wait for GitHub runner with label to become online
# -----------------------------------------------------------------------------
wait_for_runner() {
  local token="${GH_RUNNER_WAIT_TOKEN:-${GITHUB_TOKEN:-}}"
  local repo="${GH_REPO:-${GITHUB_REPOSITORY:-}}"
  local label="${RUNNER_LABEL:-gatk-rs-genomewide}"
  local timeout="${RUNNER_WAIT_TIMEOUT_SEC:-600}"
  local grace="${CLOUD_VM_BOOT_GRACE_SEC:-90}"

  if [[ -z "${token}" || -z "${repo}" ]]; then
    log "no GH_RUNNER_WAIT_TOKEN/GITHUB_TOKEN+GH_REPO — sleeping ${grace}s boot grace"
    sleep "${grace}"
    return 0
  fi

  log "waiting up to ${timeout}s for runner label=${label} on ${repo}"
  local deadline=$((SECONDS + timeout))
  while (( SECONDS < deadline )); do
    local json status
    json="$(curl -fsSL \
      -H "Authorization: Bearer ${token}" \
      -H "Accept: application/vnd.github+json" \
      "https://api.github.com/repos/${repo}/actions/runners")" || {
      sleep 10
      continue
    }
    status="$(printf '%s' "${json}" | python3 -c '
import json,sys
label=sys.argv[1]
data=json.load(sys.stdin)
for r in data.get("runners", []):
    labels={x.get("name") for x in r.get("labels", [])}
    if label in labels:
        print(r.get("status","") + " " + r.get("busy","") )
        break
else:
    print("missing")
' "${label}")"
    log "runner poll: ${status}"
    if [[ "${status}" == online* ]]; then
      log "runner online"
      return 0
    fi
    sleep 10
  done
  echo "[cloud_vmctl] ERROR: timed out waiting for runner label=${label}" >&2
  exit 1
}

# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------
case "${provider}" in
  hetzner)
    case "${cmd}" in
      status) hetzner_status ;;
      start)
        hetzner_start
        wait_for_runner
        ;;
      stop) hetzner_stop ;;
      *)
        echo "[cloud_vmctl] unknown command: ${cmd}" >&2
        exit 2
        ;;
    esac
    ;;
  aws)
    case "${cmd}" in
      status) aws_status ;;
      start)
        aws_start
        wait_for_runner
        ;;
      stop) aws_stop ;;
      *)
        echo "[cloud_vmctl] unknown command: ${cmd}" >&2
        exit 2
        ;;
    esac
    ;;
  *)
    echo "[cloud_vmctl] unsupported CLOUD_PROVIDER=${provider} (use hetzner|aws)" >&2
    exit 2
    ;;
esac

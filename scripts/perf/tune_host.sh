#!/usr/bin/env bash
# CPU tuning helpers for the dedicated gatk-rs-benchmark host.
# See docs/ci/PERF_BENCHMARK_HOST.md
set -euo pipefail

cmd="${1:-status}"

log() { echo "[tune-host] $*"; }

status() {
  echo "=== uname ==="
  uname -a
  echo "=== lscpu (summary) ==="
  if command -v rg >/dev/null 2>&1; then
    lscpu 2>/dev/null | rg -n 'Architecture|Model name|CPU\(s\)|Thread|Core|Socket|Flags|Vulnerability' || lscpu 2>/dev/null | head -40
  else
    lscpu 2>/dev/null | grep -E 'Architecture|Model name|CPU\(s\)|Thread|Core|Socket|Flags|Vulnerability' || lscpu 2>/dev/null | head -40
  fi
  echo "=== governors ==="
  if compgen -G '/sys/devices/system/cpu/cpu*/cpufreq/scaling_governor' >/dev/null; then
    for g in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
      echo "$g=$(cat "$g" 2>/dev/null || echo n/a)"
    done | sort -u
  else
    echo "(no cpufreq sysfs — common on some cloud images; document BIOS/host control)"
  fi
  echo "=== SMT / HT ==="
  if [[ -r /sys/devices/system/cpu/smt/control ]]; then
    echo "smt/control=$(cat /sys/devices/system/cpu/smt/control)"
  else
    echo "(smt control not exposed)"
  fi
  echo "=== isolcpus / cmdline ==="
  if [[ -r /proc/cmdline ]]; then
    if command -v rg >/dev/null 2>&1; then
      tr ' ' '\n' </proc/cmdline | rg 'isolcpus|nohz_full|rcu_nocbs|mitigations' || echo "(no isolcpus-style params)"
    else
      tr ' ' '\n' </proc/cmdline | grep -E 'isolcpus|nohz_full|rcu_nocbs|mitigations' || echo "(no isolcpus-style params)"
    fi
  fi
}

apply() {
  if [[ "$(id -u)" -ne 0 ]]; then
    echo "[tune-host] ERROR: 'apply' requires root (sudo ./scripts/perf/tune_host.sh apply)" >&2
    exit 1
  fi

  log "setting CPU governor → performance"
  if compgen -G '/sys/devices/system/cpu/cpu*/cpufreq/scaling_governor' >/dev/null; then
    for g in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
      echo performance >"$g" || true
    done
  elif command -v cpupower >/dev/null 2>&1; then
    cpupower frequency-set -g performance || true
  else
    log "WARN: could not set governor (install linux-tools / enable cpufreq)"
  fi

  if [[ -w /sys/devices/system/cpu/smt/control ]]; then
    cur="$(cat /sys/devices/system/cpu/smt/control)"
    if [[ "${cur}" != "off" && "${cur}" != "notsupported" && "${cur}" != "forceoff" ]]; then
      log "disabling SMT/HT (was ${cur})"
      echo off >/sys/devices/system/cpu/smt/control || log "WARN: could not disable SMT"
    else
      log "SMT already ${cur}"
    fi
  else
    log "WARN: SMT sysfs not writable — disable HT in instance BIOS/CPU options if available"
  fi

  log "done — re-run: $0 status"
  status
}

pin_help() {
  cat <<'EOF'
Pin a benchmark to dedicated cores (example: physical cores 0-3 after HT off):

  taskset -c 0-3 ./scripts/perf/run_pairhmm_speedup.sh

Or via the suite wrapper (PERF_CPU_LIST):

  PERF_CPU_LIST=0-3 ./scripts/perf/run_dedicated_benchmark_suite.sh
EOF
}

case "${cmd}" in
  status) status ;;
  apply) apply ;;
  pin-help) pin_help ;;
  -h|--help)
    echo "Usage: tune_host.sh status|apply|pin-help"
    exit 0
    ;;
  *)
    echo "unknown command: ${cmd}" >&2
    exit 2
    ;;
esac

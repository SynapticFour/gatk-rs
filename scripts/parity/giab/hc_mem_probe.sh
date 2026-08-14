#!/usr/bin/env bash
# Background memory / cgroup sampler for GIAB HC jobs (Linux CI).
# Writes one line every GIAB_MEM_PROBE_SEC (default 5) to a log file so a hard
# runner kill still leaves the last samples in the job log / artifact.
#
# Usage:
#   GIAB_MEM_PROBE_LOG=path/to.log ./scripts/parity/giab/hc_mem_probe.sh &
#   probe_pid=$!
#   ... run HC ...
#   kill "${probe_pid}" 2>/dev/null || true
set -uo pipefail

interval_sec="${GIAB_MEM_PROBE_SEC:-5}"
log="${GIAB_MEM_PROBE_LOG:-${PWD}/parity/giab/runs/ci/hc-mem-probe.log}"
mkdir -p "$(dirname "${log}")"

meminfo_kb() {
  local key="$1"
  if [[ -r /proc/meminfo ]]; then
    awk -v k="${key}" '$1 == k ":" {print $2; exit}' /proc/meminfo
  fi
}

cgroup_mem_current() {
  if [[ -r /sys/fs/cgroup/memory.current ]]; then
    cat /sys/fs/cgroup/memory.current
  elif [[ -r /sys/fs/cgroup/memory/memory.usage_in_bytes ]]; then
    cat /sys/fs/cgroup/memory/memory.usage_in_bytes
  else
    echo "na"
  fi
}

cgroup_mem_max() {
  if [[ -r /sys/fs/cgroup/memory.max ]]; then
    cat /sys/fs/cgroup/memory.max
  elif [[ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]]; then
    cat /sys/fs/cgroup/memory/memory.limit_in_bytes
  else
    echo "na"
  fi
}

cgroup_mem_peak() {
  if [[ -r /sys/fs/cgroup/memory.peak ]]; then
    cat /sys/fs/cgroup/memory.peak
  else
    echo "na"
  fi
}

# Emit to both the probe log and stdout (job log) so hard kills still show a trail.
sample() {
  local ts mem_total mem_avail swap_total swap_free cg_cur cg_max cg_peak
  local rss_kb=0 cmd="-"
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  mem_total="$(meminfo_kb MemTotal)"
  mem_avail="$(meminfo_kb MemAvailable)"
  swap_total="$(meminfo_kb SwapTotal)"
  swap_free="$(meminfo_kb SwapFree)"
  mem_total="${mem_total:-na}"
  mem_avail="${mem_avail:-na}"
  swap_total="${swap_total:-0}"
  swap_free="${swap_free:-0}"
  cg_cur="$(cgroup_mem_current)"
  cg_max="$(cgroup_mem_max)"
  cg_peak="$(cgroup_mem_peak)"

  # Prefer the true HC process Peak, not the newest match (docker wrappers often
  # match `HaplotypeCaller` with tiny host RSS while JVM / gatk-rs is the real Peak).
  local pid="" rss_kb=0 cmd="-"
  local cand best_rss=0 best_pid="" best_cmd="-"
  local any_rss=0 any_pid="" any_cmd="-"
  while IFS= read -r cand; do
    [[ -z "${cand}" || ! -r "/proc/${cand}/status" ]] && continue
    local this_rss this_cmd
    this_rss="$(awk '/^VmRSS:/ {print $2}' "/proc/${cand}/status" 2>/dev/null || echo 0)"
    this_cmd="$(tr '\0' ' ' <"/proc/${cand}/cmdline" 2>/dev/null | cut -c1-120 || echo "?")"
    if [[ "${this_rss}" =~ ^[0-9]+$ ]] && (( this_rss >= any_rss )); then
      any_rss="${this_rss}"
      any_pid="${cand}"
      any_cmd="${this_cmd}"
    fi
    # Skip docker/cli wrappers when a heavier java/gatk-rs child exists.
    if [[ "${this_cmd}" == docker\ * || "${this_cmd}" == */docker\ * ]]; then
      continue
    fi
    if [[ "${this_rss}" =~ ^[0-9]+$ ]] && (( this_rss >= best_rss )); then
      best_rss="${this_rss}"
      best_pid="${cand}"
      best_cmd="${this_cmd}"
    fi
  done < <(pgrep -f '/gatk-rs( |$)|HaplotypeCaller|/gatk[[:space:]]' 2>/dev/null || true)
  if [[ -n "${best_pid}" ]]; then
    pid="${best_pid}"
    rss_kb="${best_rss}"
    cmd="${best_cmd}"
  elif [[ -n "${any_pid}" ]]; then
    pid="${any_pid}"
    rss_kb="${any_rss}"
    cmd="${any_cmd}"
  fi
  local load
  load="$(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null || echo "?")"
  local swap_used="na"
  if [[ "${swap_total}" =~ ^[0-9]+$ && "${swap_free}" =~ ^[0-9]+$ ]]; then
    swap_used=$((swap_total - swap_free))
  fi

  local line
  line=$(printf \
    '[hc-mem] %s mem_avail_kb=%s/%s swap_used_kb=%s/%s cgroup_cur=%s cgroup_max=%s cgroup_peak=%s pid=%s rss_kb=%s load=%s cmd=%s' \
    "${ts}" \
    "${mem_avail}" "${mem_total}" \
    "${swap_used}" "${swap_total}" \
    "${cg_cur}" "${cg_max}" "${cg_peak}" \
    "${pid:-none}" "${rss_kb}" "${load}" "${cmd}")
  echo "${line}" | tee -a "${log}"
}

# Final sample on TERM/INT/EXIT so short-lived kills (and CI `trap cleanup EXIT`)
# still record Peak near process teardown instead of truncating mid-climb.
_probe_stop=0
_on_probe_signal() {
  if [[ "${_probe_stop}" -ne 0 ]]; then
    return 0
  fi
  _probe_stop=1
  trap - TERM INT EXIT
  sample || true
  echo "[hc-mem] probe stop signal=${1:-exit}" | tee -a "${log}" || true
  exit 0
}
trap '_on_probe_signal TERM' TERM
trap '_on_probe_signal INT' INT
trap '_on_probe_signal EXIT' EXIT

echo "[hc-mem] probe start interval=${interval_sec}s log=${log}" | tee -a "${log}"
sample
while [[ "${_probe_stop}" -eq 0 ]]; do
  sleep "${interval_sec}"
  sample || true
done
_on_probe_signal loop_end

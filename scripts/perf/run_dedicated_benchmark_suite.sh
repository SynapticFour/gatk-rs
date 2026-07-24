#!/usr/bin/env bash
# Orchestrate reproducible perf measurements on the dedicated benchmark host.
# Primary product: fair HC comparison (5× repeats, 3 region sizes, median±stdev).
#
# Env:
#   PERF_CPU_LIST   — optional taskset list, e.g. 0-3
#   PERF_REPEATS=5
#   PERF_SKIP_FAIR=1 — skip fair comparison (not for published numbers)
#   PERF_SKIP_MEMORY=1 / PERF_SKIP_PAIRHMM=1 — optional extras
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
run_dir="${repo_root}/docs/perf/runs/dedicated_${stamp}"
mkdir -p "${run_dir}"
export PERF_RUN_DIR="${run_dir}"
export PERF_RUNNER_LABEL="${PERF_RUNNER_LABEL:-gatk-rs-benchmark}"

echo "[dedicated-bench] stamp=${stamp} run_dir=${run_dir}"

if [[ "${PERF_SKIP_FAIR:-0}" != "1" ]]; then
  echo "[dedicated-bench] fair HC comparison (primary)…"
  # Nested under dedicated stamp so artifacts stay together; fair script also
  # writes docs/perf/FAIR_HC_COMPARISON.md + fair_hc_comparison_latest.json.
  PERF_RUN_DIR="${run_dir}/fair" \
    "${script_dir}/run_fair_hc_comparison.sh" | tee "${run_dir}/fair_hc.stdout"
  cp -f "${repo_root}/docs/perf/FAIR_HC_COMPARISON.md" "${run_dir}/FAIR_HC_COMPARISON.md" || true
  cp -f "${repo_root}/docs/perf/fair_hc_comparison_latest.json" \
    "${run_dir}/fair_hc_comparison_latest.json" || true
else
  echo "[dedicated-bench] WARN: skipping fair comparison (PERF_SKIP_FAIR=1)"
  "${script_dir}/capture_host_specs.sh"
fi

if [[ "${PERF_SKIP_MEMORY:-0}" != "1" ]]; then
  echo "[dedicated-bench] optional Peak-RSS micro-profile (tiny fixture)…"
  if [[ -n "${PERF_CPU_LIST:-}" ]] && command -v taskset >/dev/null 2>&1; then
    taskset -c "${PERF_CPU_LIST}" "${script_dir}/run_hc_memory_profile.sh" \
      | tee "${run_dir}/memory_profile.stdout"
  else
    "${script_dir}/run_hc_memory_profile.sh" | tee "${run_dir}/memory_profile.stdout"
  fi
  cp -f "${repo_root}/docs/perf/HC_MEMORY_PROFILE.md" "${run_dir}/HC_MEMORY_PROFILE.md" || true
  cp -f "${repo_root}/docs/perf/hc_memory_profile_latest.json" \
    "${run_dir}/hc_memory_profile_latest.json" || true
fi

if [[ "${PERF_SKIP_PAIRHMM:-0}" != "1" ]]; then
  echo "[dedicated-bench] optional PairHMM Criterion microbench…"
  if [[ -n "${PERF_CPU_LIST:-}" ]] && command -v taskset >/dev/null 2>&1; then
    taskset -c "${PERF_CPU_LIST}" "${script_dir}/run_pairhmm_speedup.sh" \
      | tee "${run_dir}/pairhmm_speedup.stdout"
  else
    "${script_dir}/run_pairhmm_speedup.sh" | tee "${run_dir}/pairhmm_speedup.stdout"
  fi
  cp -f "${repo_root}/docs/perf/PAIRHMM_SPEEDUP.md" "${run_dir}/PAIRHMM_SPEEDUP.md" || true
  cp -f "${repo_root}/docs/perf/pairhmm_speedup_latest.json" \
    "${run_dir}/pairhmm_speedup_latest.json" || true
fi

# Update public performance dashboard history (committed by workflow).
if [[ -f "${repo_root}/docs/perf/fair_hc_comparison_latest.json" ]]; then
  run_url=""
  if [[ -n "${GITHUB_SERVER_URL:-}" && -n "${GITHUB_REPOSITORY:-}" && -n "${GITHUB_RUN_ID:-}" ]]; then
    run_url="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${GITHUB_RUN_ID}"
  fi
  python3 "${script_dir}/update_perf_dashboard.py" \
    --summary-json "${repo_root}/docs/perf/fair_hc_comparison_latest.json" \
    --commit-sha "${GITHUB_SHA:-local}" \
    --workflow-run-url "${run_url}"
fi

cat >"${run_dir}/README.md" <<EOF
# Dedicated benchmark suite run

**UTC:** \`${stamp}\`  
**Primary report:** [FAIR_HC_COMPARISON.md](FAIR_HC_COMPARISON.md)  
**Host specs:** [host_specs.md](fair/host_specs.md) (or [HOST_SPECS.md](../../HOST_SPECS.md))

Primary Java baseline for speedups: **FASTEST_AVAILABLE**.
EOF

echo "[dedicated-bench] done → ${run_dir}"

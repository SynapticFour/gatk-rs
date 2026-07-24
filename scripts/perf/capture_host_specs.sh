#!/usr/bin/env bash
# Capture hardware / OS facts for publishing next to every benchmark number.
# Output:
#   docs/perf/HOST_SPECS.md
#   docs/perf/host_specs_latest.json
#   optional: docs/perf/runs/<stamp>/host_specs.{md,json} when PERF_RUN_DIR is set
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_root="${repo_root}/docs/perf"
mkdir -p "${out_root}"

model="$(lscpu 2>/dev/null | awk -F: '/Model name/ {gsub(/^[ \t]+/,"",$2); print $2; exit}' || true)"
cpus="$(nproc 2>/dev/null || echo unknown)"
mem_kb="$(awk '/MemTotal/ {print $2}' /proc/meminfo 2>/dev/null || echo 0)"
mem_gib="$(awk -v k="${mem_kb}" 'BEGIN { printf "%.1f", k/1024/1024 }')"
kernel="$(uname -srm)"
flags="$(awk '/^flags/{print; exit}' /proc/cpuinfo 2>/dev/null || true)"
has_avx2=0
has_avx512=0
[[ "${flags}" == *avx2* ]] && has_avx2=1
[[ "${flags}" == *avx512f* ]] && has_avx512=1
governor="n/a"
if compgen -G '/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor' >/dev/null; then
  governor="$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
fi
smt="n/a"
[[ -r /sys/devices/system/cpu/smt/control ]] && smt="$(cat /sys/devices/system/cpu/smt/control)"
hostname_s="$(hostname -s 2>/dev/null || hostname || echo unknown)"
role_label="${PERF_RUNNER_LABEL:-gatk-rs-benchmark}"

python3 - "${out_root}" "${stamp}" "${model}" "${cpus}" "${mem_gib}" "${kernel}" \
  "${has_avx2}" "${has_avx512}" "${governor}" "${smt}" "${hostname_s}" "${role_label}" \
  "${flags}" <<'PY'
import json, pathlib, sys
(
    out_root, stamp, model, cpus, mem_gib, kernel,
    has_avx2, has_avx512, governor, smt, hostname_s, role_label, flags,
) = sys.argv[1:]

summary = {
    "stamp_utc": stamp,
    "hostname": hostname_s,
    "runner_label": role_label,
    "cpu_model": model or "unknown",
    "logical_cpus": int(cpus) if str(cpus).isdigit() else cpus,
    "mem_gib": float(mem_gib),
    "kernel": kernel,
    "governor": governor,
    "smt_control": smt,
    "simd": {
        "avx2": bool(int(has_avx2)),
        "avx512f": bool(int(has_avx512)),
    },
    "cpu_flags_line": flags.strip(),
    "note": "Publish this block next to every timing / Peak-RSS number from this host.",
}

out = pathlib.Path(out_root)
json_path = out / "host_specs_latest.json"
json_path.write_text(json.dumps(summary, indent=2) + "\n")

md = f"""# Performance host specifications

**Captured (UTC):** `{stamp}`  
**Hostname:** `{hostname_s}`  
**Runner label:** `{role_label}`

| Field | Value |
|-------|--------|
| CPU model | {summary['cpu_model']} |
| Logical CPUs (`nproc`) | {summary['logical_cpus']} |
| RAM | {summary['mem_gib']} GiB |
| Kernel | `{kernel}` |
| Governor (cpu0) | `{governor}` |
| SMT / HT control | `{smt}` |
| AVX2 | {"yes" if summary["simd"]["avx2"] else "NO"} |
| AVX-512F | {"yes" if summary["simd"]["avx512f"] else "no"} |

Setup doctrine: [`docs/ci/PERF_BENCHMARK_HOST.md`](../ci/PERF_BENCHMARK_HOST.md).

> Numbers measured on GitHub-hosted runners, laptops, or the `gatk-rs-genomewide`
> correctness VM are **not** interchangeable with this host.
"""
md_path = out / "HOST_SPECS.md"
md_path.write_text(md)
print(f"Wrote {md_path}")
print(f"Wrote {json_path}")
PY

if [[ -n "${PERF_RUN_DIR:-}" ]]; then
  mkdir -p "${PERF_RUN_DIR}"
  cp -f "${out_root}/HOST_SPECS.md" "${PERF_RUN_DIR}/host_specs.md"
  cp -f "${out_root}/host_specs_latest.json" "${PERF_RUN_DIR}/host_specs.json"
  echo "[host-specs] also copied into ${PERF_RUN_DIR}"
fi

#!/usr/bin/env bash
# Quantify cargo bench profile settings vs release-equivalent builds.
#
# Does NOT modify Cargo.toml. Overrides via CARGO_PROFILE_BENCH_* / RUSTFLAGS only.
# Writes docs/perf/runs/bench_profile_<stamp>/ with per-config Criterion stdout
# and a summary TSV + Markdown.
#
# Env:
#   BENCH_MATRIX_QUICK=1     — short Criterion times (default 1)
#   BENCH_MATRIX_FILTERS     — space-separated Criterion filter regexes
#   BENCH_MATRIX_CONFIGS     — space-separated config ids
#
# Example:
#   BENCH_MATRIX_QUICK=1 ./scripts/perf/run_bench_profile_matrix.sh
#   BENCH_MATRIX_CONFIGS='current_bench release_equiv' ./scripts/perf/run_bench_profile_matrix.sh
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"
chmod +x "${script_dir}/run_bench_profile_matrix.sh" 2>/dev/null || true

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
run_dir="${repo_root}/docs/perf/runs/bench_profile_${stamp}"
mkdir -p "${run_dir}"

host_cpu="$(sysctl -n machdep.cpu.brand_string 2>/dev/null || true)"
rustc_v="$(rustc -vV | tr '\n' '; ')"
git_sha="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

cat >"${run_dir}/HOST.md" <<EOF
# Bench profile matrix host

- UTC: \`${stamp}\`
- git: \`${git_sha}\`
- uname: \`$(uname -srm)\`
- cpu: \`${host_cpu:-unknown}\`
- rustc: \`${rustc_v}\`
EOF

# Config id → env overrides applied on top of Cargo.toml [profile.bench].
# current_bench: lto=false, cgu=16, debug=true (as committed)
apply_config_env() {
  local name="$1"
  unset CARGO_PROFILE_BENCH_LTO CARGO_PROFILE_BENCH_CODEGEN_UNITS CARGO_PROFILE_BENCH_DEBUG || true
  export RUSTFLAGS="${RUSTFLAGS_BASE:-}"
  case "${name}" in
    current_bench)
      # Inherit Cargo.toml [profile.bench] as-is.
      ;;
    release_equiv)
      export CARGO_PROFILE_BENCH_LTO=fat
      export CARGO_PROFILE_BENCH_CODEGEN_UNITS=1
      export CARGO_PROFILE_BENCH_DEBUG=0
      ;;
    release_native)
      export CARGO_PROFILE_BENCH_LTO=fat
      export CARGO_PROFILE_BENCH_CODEGEN_UNITS=1
      export CARGO_PROFILE_BENCH_DEBUG=0
      export RUSTFLAGS="${RUSTFLAGS_BASE:-} -C target-cpu=native"
      ;;
    lto_thin)
      export CARGO_PROFILE_BENCH_LTO=thin
      export CARGO_PROFILE_BENCH_CODEGEN_UNITS=1
      export CARGO_PROFILE_BENCH_DEBUG=0
      ;;
    no_lto_cgu1)
      export CARGO_PROFILE_BENCH_LTO=false
      export CARGO_PROFILE_BENCH_CODEGEN_UNITS=1
      export CARGO_PROFILE_BENCH_DEBUG=0
      ;;
    no_lto_cgu16_nodebug)
      export CARGO_PROFILE_BENCH_LTO=false
      export CARGO_PROFILE_BENCH_CODEGEN_UNITS=16
      export CARGO_PROFILE_BENCH_DEBUG=0
      ;;
    *)
      echo "unknown config: ${name}" >&2
      return 1
      ;;
  esac
}

ALL_CONFIGS=(current_bench release_equiv release_native lto_thin no_lto_cgu1 no_lto_cgu16_nodebug)
if [[ -n "${BENCH_MATRIX_CONFIGS:-}" ]]; then
  # shellcheck disable=SC2206
  ALL_CONFIGS=(${BENCH_MATRIX_CONFIGS})
fi

if [[ -n "${BENCH_MATRIX_FILTERS:-}" ]]; then
  # shellcheck disable=SC2206
  FILTERS=(${BENCH_MATRIX_FILTERS})
elif [[ "${BENCH_MATRIX_QUICK:-1}" == "1" ]]; then
  FILTERS=(
    "pairhmm_logless_simd"
    "smith_waterman_align/soft_clip/128x96"
    "assembly_graph_depth/threading_build_medium_k10"
    "genotype_dense_ad/multipass_ad/R256"
  )
else
  FILTERS=(
    "pairhmm_logless_simd"
    "smith_waterman_align/soft_clip"
    "assembly_graph_depth/threading_build"
    "genotype_dense_ad/multipass_ad"
  )
fi

crit_args=(--warm-up-time 0.3 --measurement-time 1.5 --sample-size 25)
base_target="${CARGO_TARGET_DIR:-${run_dir}/target}"

summary_tsv="${run_dir}/summary.tsv"
printf 'config\tfilter\tbench_id\ttime_ns_mid\traw_line\n' >"${summary_tsv}"

parse_stdout() {
  local name="$1" filt="$2" out="$3"
  python3 - "${name}" "${filt}" "${out}" "${summary_tsv}" <<'PY'
import re, sys
name, filt, path, tsv = sys.argv[1:5]
text = open(path, encoding="utf-8", errors="replace").read()
unit = {"ns": 1.0, "µs": 1e3, "us": 1e3, "ms": 1e6, "s": 1e9}

def to_ns(num: str, u: str) -> float:
    u = u.replace("μ", "µ")
    if u == "μs":
        u = "µs"
    return float(num) * unit.get(u, float("nan"))

# Match Criterion time lines: time:   [1.2345 µs 1.2346 µs 1.2347 µs]
time_re = re.compile(
    r"time:\s+\[\s*([0-9.]+)\s*(ns|µs|us|ms|s)\s+([0-9.]+)\s*(ns|µs|us|ms|s)\s+([0-9.]+)\s*(ns|µs|us|ms|s)\s*\]"
)
# Preceding non-empty line is usually the bench id
lines = text.splitlines()
found = 0
with open(tsv, "a", encoding="utf-8") as f:
    for i, line in enumerate(lines):
        m = time_re.search(line)
        if not m:
            continue
        mid_ns = to_ns(m.group(3), m.group(4))
        bid = "(unknown)"
        for j in range(i - 1, max(-1, i - 6), -1):
            prev = lines[j].strip()
            if prev and not prev.startswith("Benchmarking") and "time:" not in prev:
                bid = prev
                break
        f.write(f"{name}\t{filt}\t{bid}\t{mid_ns}\t{line.strip()}\n")
        found += 1
if found == 0:
    with open(tsv, "a", encoding="utf-8") as f:
        f.write(f"{name}\t{filt}\tPARSE_FAIL\t\t\n")
PY
}

run_one() {
  local name="$1"
  local cfg_dir="${run_dir}/${name}"
  mkdir -p "${cfg_dir}"
  local target_dir="${base_target}/${name}"
  mkdir -p "${target_dir}"

  apply_config_env "${name}"
  {
    echo "config=${name}"
    echo "CARGO_PROFILE_BENCH_LTO=${CARGO_PROFILE_BENCH_LTO:-<toml>}"
    echo "CARGO_PROFILE_BENCH_CODEGEN_UNITS=${CARGO_PROFILE_BENCH_CODEGEN_UNITS:-<toml>}"
    echo "CARGO_PROFILE_BENCH_DEBUG=${CARGO_PROFILE_BENCH_DEBUG:-<toml>}"
    echo "RUSTFLAGS=${RUSTFLAGS:-}"
  } | tee "${cfg_dir}/env.txt" | tee -a "${run_dir}/matrix.log"

  for filt in "${FILTERS[@]}"; do
    local out="${cfg_dir}/$(echo "${filt}" | tr '/ ' '__').stdout"
    echo "[matrix] ${name} :: ${filt}" | tee -a "${run_dir}/matrix.log"
    local -a pkg_bench=()
    local -a feat_args=()
    case "${filt}" in
      pairhmm*) pkg_bench=(-p gatk-haplotypecaller --bench pairhmm) ;;
      smith_waterman*) pkg_bench=(-p gatk-haplotypecaller --bench smith_waterman) ;;
      assembly_graph*) pkg_bench=(-p gatk-haplotypecaller --bench assembly_graph) ;;
      genotype*)
        pkg_bench=(-p gatk-haplotypecaller --bench genotype_assign)
        feat_args=(--features parity_harness)
        ;;
      *)
        echo "[matrix] skip unknown filter: ${filt}" >&2
        continue
        ;;
    esac

    set +e
    CARGO_TARGET_DIR="${target_dir}" \
      cargo bench \
        "${pkg_bench[@]}" \
        ${feat_args[@]+"${feat_args[@]}"} \
        --locked \
        -- "${filt}" "${crit_args[@]}" \
      >"${out}" 2>&1
    local rc=$?
    set -e
    if [[ ${rc} -ne 0 ]]; then
      echo "[matrix] FAIL rc=${rc} → ${out}" | tee -a "${run_dir}/matrix.log"
      printf '%s\t%s\tFAIL\t\t\n' "${name}" "${filt}" >>"${summary_tsv}"
      continue
    fi
    parse_stdout "${name}" "${filt}" "${out}"
  done
}

for name in "${ALL_CONFIGS[@]}"; do
  run_one "${name}"
done

python3 - "${summary_tsv}" "${run_dir}/SUMMARY.md" <<'PY'
import csv, collections, sys
from pathlib import Path
tsv, md_path = sys.argv[1:3]
rows = list(csv.DictReader(open(tsv), delimiter="\t"))
by_filt = collections.defaultdict(list)
for r in rows:
    raw = (r.get("time_ns_mid") or "").strip()
    if not raw:
        continue
    try:
        t = float(raw)
    except ValueError:
        continue
    by_filt[r["filter"]].append((r["config"], r["bench_id"], t))

def fmt(ns: float) -> str:
    if ns >= 1e9: return f"{ns/1e9:.3f} s"
    if ns >= 1e6: return f"{ns/1e6:.3f} ms"
    if ns >= 1e3: return f"{ns/1e3:.3f} µs"
    return f"{ns:.1f} ns"

lines = ["# Bench profile matrix summary", "", f"Source: `{tsv}`", ""]
for filt, items in sorted(by_filt.items()):
    lines += [f"## `{filt}`", "",
              "| config | bench_id | time (mid) | vs current_bench |",
              "|--------|----------|-----------:|-----------------:|"]
    base = next((t for c, _, t in items if c == "current_bench"), None)
    for cfg, bid, t in items:
        delta = f"{t/base:.3f}×" if base and base > 0 else "—"
        lines.append(f"| `{cfg}` | `{bid}` | {fmt(t)} | {delta} |")
    lines.append("")
Path(md_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
print("wrote", md_path)
PY

echo "[matrix] done → ${run_dir}"
echo "[matrix] summary → ${run_dir}/SUMMARY.md"

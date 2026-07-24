#!/usr/bin/env bash
# Reproducible PairHMM Logless scalar vs SIMD microbench + HC smoke.
#
# Output:
#   docs/perf/PAIRHMM_SPEEDUP.md
#   docs/perf/pairhmm_speedup_latest.json
#   docs/perf/runs/pairhmm_<stamp>/
#
# IMPORTANT: Do not invent README speedup claims without this report.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_root="${repo_root}/docs/perf"
run_dir="${out_root}/runs/pairhmm_${stamp}"
mkdir -p "${run_dir}"

export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
export PATH="${HOME}/.cargo/bin:/usr/bin:/bin:/opt/homebrew/bin:${PATH}"

echo "[pairhmm-speedup] building release gatk-cli…"
(
  cd "${repo_root}"
  cargo build -p gatk-cli --release --locked
)

target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
RUST_BIN="${target_dir}/release/gatk-rs"
REF="${repo_root}/parity/fixtures/reference.fa"
BAM="${repo_root}/parity/fixtures/sample.bam"
INTERVAL="chr1:1-32"

echo "[pairhmm-speedup] Criterion microbench (logless_simd group)…"
bench_log="${run_dir}/criterion.txt"
(
  cd "${repo_root}"
  cargo bench -p gatk-haplotypecaller --bench pairhmm --locked -- --warm-up-time 0.5 --measurement-time 2.0 pairhmm_logless_simd
) >"${bench_log}" 2>&1 || true

echo "[pairhmm-speedup] HC smoke LOG10 vs SIMD on fixture…"
log10_vcf="${run_dir}/hc_log10.vcf"
simd_vcf="${run_dir}/hc_simd.vcf"
/usr/bin/time -l "${RUST_BIN}" HaplotypeCaller \
  -R "${REF}" -I "${BAM}" -O "${log10_vcf}" -L "${INTERVAL}" \
  --pair-hmm LOG10_PAIRHMM \
  >"${run_dir}/hc_log10.stdout" 2>"${run_dir}/hc_log10.time" || true
/usr/bin/time -l "${RUST_BIN}" HaplotypeCaller \
  -R "${REF}" -I "${BAM}" -O "${simd_vcf}" -L "${INTERVAL}" \
  --pair-hmm SIMD \
  >"${run_dir}/hc_simd.stdout" 2>"${run_dir}/hc_simd.time" || true

rustc_v="$(rustc --version 2>/dev/null || echo unknown)"
cargo_v="$(cargo --version 2>/dev/null || echo unknown)"
git_sha="$(git -C "${repo_root}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
uname_s="$(uname -srm)"
backend="$(
  cd "${repo_root}"
  cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test parse_pair_hmm --locked -- --nocapture 2>/dev/null | head -1
  python3 - <<'PY'
import platform
print(platform.machine())
PY
)"

# Parse criterion medians roughly from text (best-effort).
python3 - "${run_dir}" "${out_root}" "${stamp}" \
  "${rustc_v}" "${cargo_v}" "${uname_s}" "${git_sha}" \
  "${RUST_BIN}" "${bench_log}" <<'PY'
import json, pathlib, re, sys

run_dir, out_root, stamp, rustc_v, cargo_v, uname_s, git_sha, rust_bin, bench_log = sys.argv[1:]
text = pathlib.Path(bench_log).read_text(encoding="utf-8", errors="replace")
# Criterion lines like: simd_haps/32 … time: [.. 1.2345 ms ..]
pat = re.compile(
    r"(logless_scalar_haps|simd_haps|simd_f32_haps|log10_scalar_haps)/(\d+)\s+.*?time:\s+\[[^\]]*?([0-9.]+)\s+(ns|us|ms|s)",
    re.S,
)
rows = []
for m in pat.finditer(text):
    name, n, val, unit = m.group(1), int(m.group(2)), float(m.group(3)), m.group(4)
    scale = {"ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1.0}[unit]
    rows.append({"impl": name, "hap_count": n, "seconds": val * scale, "raw": m.group(0)[:120]})

# Speedup: scalar / simd for matching hap_count
by = {}
for r in rows:
    by.setdefault(r["hap_count"], {})[r["impl"]] = r["seconds"]
speedups = []
for n, d in sorted(by.items()):
    s = d.get("logless_scalar_haps")
    v = d.get("simd_haps")
    if s and v and v > 0:
        speedups.append({"hap_count": n, "scalar_s": s, "simd_s": v, "speedup": s / v})

host_features = []
if "arm64" in uname_s.lower() or "aarch64" in uname_s.lower():
    host_features.append("neon")
if "x86_64" in uname_s.lower():
    host_features.append("avx2_runtime_detect")

summary = {
    "stamp_utc": stamp,
    "host": uname_s,
    "git_sha": git_sha,
    "rustc": rustc_v,
    "cargo": cargo_v,
    "binary": rust_bin,
    "build": "cargo build -p gatk-cli --release --locked",
    "bench": "cargo bench -p gatk-haplotypecaller --bench pairhmm -- pairhmm_logless_simd",
    "host_features": host_features,
    "criterion_rows": rows,
    "speedups_scalar_over_simd": speedups,
    "default_pair_hmm": "LOG10_PAIRHMM",
    "default_note": "SIMD not promoted to default: GIAB/hap.py F1 gate not signed in this run; production remains Log10PairHMM.",
    "hc_smoke": {
        "interval": "chr1:1-32",
        "fixture": "parity/fixtures",
        "log10_time": "hc_log10.time",
        "simd_time": "hc_simd.time",
    },
    "unit_tests": "cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test",
}

pathlib.Path(out_root).mkdir(parents=True, exist_ok=True)
json_path = pathlib.Path(out_root) / "pairhmm_speedup_latest.json"
json_path.write_text(json.dumps(summary, indent=2) + "\n")
(pathlib.Path(run_dir) / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

speed_tbl = ""
if speedups:
    speed_tbl = "| hap_count | scalar Logless | SIMD | speedup (scalar/SIMD) |\n|---|---:|---:|---:|\n"
    for s in speedups:
        speed_tbl += f"| {s['hap_count']} | {s['scalar_s']*1e3:.3f} ms | {s['simd_s']*1e3:.3f} ms | **{s['speedup']:.2f}×** |\n"
else:
    speed_tbl = "_No Criterion medians parsed — see `criterion.txt` raw log._\n"

md = f"""# PairHMM speedup profile (reproducible)

**Generated (UTC):** `{stamp}`  
**Host:** `{uname_s}`  
**Runner:** [`scripts/perf/run_pairhmm_speedup.sh`](../../scripts/perf/run_pairhmm_speedup.sh)  
**Raw run:** `docs/perf/runs/pairhmm_{stamp}/`

## Build / versions

- rustc: `{rustc_v}`
- cargo: `{cargo_v}`
- git: `{git_sha}`
- build: `cargo build -p gatk-cli --release --locked`
- features detected (coarse): `{host_features}`

## Microbench (Criterion `pairhmm_logless_simd`)

Read length 200 bp, parity indel/GCP quals 45/45/10, baseQ 30.

{speed_tbl}

Exact command:

```bash
cargo bench -p gatk-haplotypecaller --bench pairhmm --locked -- pairhmm_logless_simd
```

## Production default

**Still `LOG10_PAIRHMM`.** SIMD is available via `--pair-hmm SIMD` / `FASTEST_AVAILABLE`,
but is **not** the HC default until a signed GIAB/hap.py run shows no F1 regression
([`docs/CLAIM_MATRIX.md`](../CLAIM_MATRIX.md) — GIAB ci-subset not yet signed).

Unit gate (SIMD vs scalar Logless):  
`cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test`

## HC smoke (fixture)

```bash
./target/release/gatk-rs HaplotypeCaller -R parity/fixtures/reference.fa \\
  -I parity/fixtures/sample.bam -O /tmp/hc.vcf -L chr1:1-32 --pair-hmm LOG10_PAIRHMM
./target/release/gatk-rs HaplotypeCaller ... --pair-hmm SIMD
```

See `{run_dir}/hc_*.time` for wall/RSS from `/usr/bin/time -l`.
"""
md_path = pathlib.Path(out_root) / "PAIRHMM_SPEEDUP.md"
md_path.write_text(md)
print(f"Wrote {md_path}")
print(f"Wrote {json_path}")
if speedups:
    best = max(speedups, key=lambda s: s["speedup"])
    print(f"Best measured speedup: {best['speedup']:.2f}× at hap_count={best['hap_count']}")
else:
    print("No speedup parsed from Criterion output")
PY

echo "[pairhmm-speedup] done → docs/perf/PAIRHMM_SPEEDUP.md"

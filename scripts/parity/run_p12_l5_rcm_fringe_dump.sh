#!/usr/bin/env bash
# P12 L5.2 — RCM locus dump on read fringe (92305500–92305640) Rust vs Java.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

ref="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
bam="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
interval="${P12_RCM_FRINGE_INTERVAL:-2:92305500-92305640}"
padding="${P12_RCM_PADDING:-100}"
out_dir="${repo_root}/parity/reports"
mkdir -p "${out_dir}"

rust_out="${out_dir}/p12_rcm_fringe.rust.tsv"
java_out="${out_dir}/p12_rcm_fringe.java.tsv"

cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
"${cargo_run[@]}" reference-confidence-locus "${ref}" "${bam}" "${interval}" "${padding}" \
  >"${rust_out}"

./scripts/parity/run_hc_full_parity_java_dump.sh reference-confidence-locus \
  "${ref}" "${bam}" "${interval}" "${padding}" \
  2>/dev/null | rg -v '^(INFO|WARN|DEBUG|SLF4J|\tat )' >"${java_out}"

echo "[p12-rcm-fringe] rust=${rust_out} java=${java_out}"
python3 - "${rust_out}" "${java_out}" <<'PY'
import csv, sys
from pathlib import Path

def load(path):
    rows = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("contig") is False and line.split("\t", 1)[0] not in ("2", "chr2"):
            if line.startswith("contig\t"):
                continue
        if line.startswith("contig\t"):
            hdr = line
            continue
        parts = line.split("\t")
        if len(parts) < 6 or not parts[1].isdigit():
            continue
        rows[int(parts[1])] = {"gq": parts[5], "dp": parts[6]}
    return rows

rust, java = load(sys.argv[1]), load(sys.argv[2])
m = [(p, java[p], rust[p]) for p in sorted(set(rust) & set(java)) if rust[p] != java[p]]
print(f"positions={len(rust)} mismatches={len(m)}")
for p, j, r in m[:10]:
    print(f"  {p}: java gq={j['gq']} dp={j['dp']} | rust gq={r['gq']} dp={r['dp']}")
if len(m) > 10:
    print(f"  ... {len(m)-10} more")
PY

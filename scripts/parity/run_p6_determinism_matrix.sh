#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

summary_json="${report_dir}/p6_determinism_matrix_summary.json"
rows_json="${report_dir}/p6_determinism_matrix_rows.json"
: > "${rows_json}"

threads=(1 2 4 8)
repeats=(1 2)

for t in "${threads[@]}"; do
  for r in "${repeats[@]}"; do
    echo "[p6-determinism] threads=${t} repeat=${r}"
    set +e
    out=$(cd "${repo_root}" && RAYON_NUM_THREADS="${t}" cargo test -p gatk-haplotypecaller --test p6_scalar_vector_equivalence_test --locked scalar_and_vectorized_pairhmm_are_equivalent 2>&1)
    code=$?
    set -e
    ok_int=0
    if [[ "${code}" -eq 0 ]]; then
      ok_int=1
    fi
    python3 - <<PY >> "${rows_json}"
import json
print(json.dumps({
  "threads": ${t},
  "repeat": ${r},
  "exit_code": ${code},
  "ok": bool(${ok_int}),
}))
PY
  done
done

python3 - "${repo_root}" <<'PY'
import json
import pathlib
import sys
repo = pathlib.Path(sys.argv[1])
rows_path = repo / "parity" / "reports" / "p6_determinism_matrix_rows.json"
summary_path = repo / "parity" / "reports" / "p6_determinism_matrix_summary.json"
rows = [json.loads(line) for line in rows_path.read_text(encoding="utf-8").splitlines() if line.strip()]
all_ok = all(r["ok"] for r in rows)
summary = {
  "label": "phase6-determinism-matrix",
  "matrix": {"threads": [1,2,4,8], "repeats": 2},
  "rows": rows,
  "pass": all_ok,
}
summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(f"[p6-determinism] wrote {summary_path}")
print(f"[p6-determinism] pass={all_ok}")
if not all_ok:
  raise SystemExit(1)
PY

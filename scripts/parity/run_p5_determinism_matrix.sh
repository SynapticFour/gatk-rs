#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

summary_json="${report_dir}/p5_determinism_matrix_summary.json"
tmp_json="${report_dir}/p5_determinism_matrix_rows.json"
: > "${tmp_json}"

threads=(1 2 4 8)
repeats=(1 2)

rows=()
for t in "${threads[@]}"; do
  for r in "${repeats[@]}"; do
    echo "[p5-determinism] threads=${t} repeat=${r}"
    log="${report_dir}/p5_determinism_t${t}_r${r}.log"
    set +e
    (cd "${repo_root}" && RAYON_NUM_THREADS="${t}" cargo test -p gatk-haplotypecaller --test p5_assembly_regression_test --locked outputs_are_stable_across_repeated_runs_and_input_order) >"${log}" 2>&1
    code=$?
    set -e
    ok_int=0
    if [[ "${code}" -eq 0 ]]; then
      ok_int=1
    else
      echo "[p5-determinism] FAIL threads=${t} repeat=${r} exit=${code} (see ${log})" >&2
      tail -n 40 "${log}" >&2 || true
    fi
    python3 - <<PY >> "${tmp_json}"
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
rows_path = repo / "parity" / "reports" / "p5_determinism_matrix_rows.json"
summary_path = repo / "parity" / "reports" / "p5_determinism_matrix_summary.json"
rows = [json.loads(line) for line in rows_path.read_text(encoding="utf-8").splitlines() if line.strip()]
all_ok = all(r["ok"] for r in rows)
summary = {
  "label": "phase5-determinism-matrix",
  "matrix": {"threads": [1,2,4,8], "repeats": 2},
  "rows": rows,
  "pass": all_ok,
}
summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(f"[p5-determinism] wrote {summary_path}")
print(f"[p5-determinism] pass={all_ok}")
if not all_ok:
  raise SystemExit(1)
PY

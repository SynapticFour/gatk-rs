#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
tmp_dir="${report_dir}/p3-region-records-diff-tmp"
mkdir -p "${tmp_dir}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
run_rust="${repo_root}/scripts/parity/run_rust_gatk.sh"

source_sam="${repo_root}/parity/fixtures/p3_region_reads.sam"
source_bam="${tmp_dir}/p3_region_reads.bam"
samtools view -bS "${source_sam}" > "${source_bam}"
samtools index "${source_bam}"

# Write one case result as valid JSON (qnames may contain noise from cargo 2>&1 capture).
write_check_json() {
  local out_path="$1"
  python3 - "$@" <<'PY'
import json, pathlib, re, sys

out_path = pathlib.Path(sys.argv[1])
label = sys.argv[2]
equal_s = sys.argv[3]  # "true" | "false" | "null"
region = sys.argv[4]
java_path = pathlib.Path(sys.argv[5])
rust_path = pathlib.Path(sys.argv[6])
reason = sys.argv[7] if len(sys.argv) > 7 else ""
java_exit = int(sys.argv[8]) if len(sys.argv) > 8 and sys.argv[8] != "" else None
rust_exit = int(sys.argv[9]) if len(sys.argv) > 9 and sys.argv[9] != "" else None

# SAM QNAME charset; drops cargo/ANSI/log lines mixed into rust stdout via 2>&1.
_QNAME = re.compile(r"^[!-?A-~]{1,254}$")

def java_qnames(path: pathlib.Path) -> list[str]:
    out = []
    if not path.is_file():
        return out
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not ln or ln.startswith("@"):
            continue
        out.append(ln.split("\t", 1)[0])
    return out

def rust_qnames(path: pathlib.Path) -> list[str]:
    out = []
    if not path.is_file():
        return out
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or not _QNAME.match(s):
            continue
        out.append(s)
    return out

payload = {
    "label": label,
    "mode": "phase3-region-records-diff",
    "region": region,
    "java_output": str(java_path),
    "rust_output": str(rust_path),
}
if equal_s == "null":
    payload["equal"] = None
    payload["skipped"] = True
    payload["reason"] = reason or "java_gatk_missing"
elif equal_s == "tool_fail":
    payload["equal"] = False
    payload["reason"] = "tool_exit_nonzero"
    payload["java_exit"] = java_exit
    payload["rust_exit"] = rust_exit
else:
    jq = java_qnames(java_path)
    rq = rust_qnames(rust_path)
    payload["equal"] = jq == rq
    payload["java_qnames"] = jq
    payload["rust_qnames"] = rq

out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print("1" if payload.get("equal") is True else ("skip" if payload.get("skipped") else "0"))
PY
}

passed=0
failed=0
skipped=0
checks_json=()

run_case() {
  local label="$1"
  local region="$2"
  local java_stdout="${report_dir}/${label}.java.stdout"
  local rust_stdout="${report_dir}/${label}.rust.stdout"
  local java_sam="${report_dir}/${label}.java.sam"
  local check_json="${report_dir}/${label}.json"

  set +e
  "${run_java}" "${java_stdout}" PrintReads -I "${source_bam}" -L "${region}" -O "${java_sam}"
  java_exit=$?
  "${run_rust}" "${rust_stdout}" ListReadsInRegion -I "${source_bam}" -L "${region}"
  rust_exit=$?
  set -e

  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    write_check_json "${check_json}" "${label}" "null" "${region}" "${java_stdout}" "${rust_stdout}" "java_gatk_missing" >/dev/null
    skipped=$((skipped + 1))
    checks_json+=("${check_json}")
    return
  fi

  if [[ "${java_exit}" -ne 0 || "${rust_exit}" -ne 0 ]]; then
    write_check_json "${check_json}" "${label}" "tool_fail" "${region}" "${java_stdout}" "${rust_stdout}" "" "${java_exit}" "${rust_exit}" >/dev/null
    failed=$((failed + 1))
    checks_json+=("${check_json}")
    return
  fi

  # Compare PrintReads SAM QNAMEs vs ListReadsInRegion lines (ignore cargo/ANSI noise in 2>&1 capture).
  result="$(write_check_json "${check_json}" "${label}" "compare" "${region}" "${java_sam}" "${rust_stdout}")"
  if [[ "${result}" == "1" ]]; then
    passed=$((passed + 1))
  else
    failed=$((failed + 1))
  fi
  checks_json+=("${check_json}")
}

run_case "p3-region-records-chr1-1-16" "chr1:1-16"
run_case "p3-region-records-chr1-17-32" "chr1:17-32"

summary_json="${report_dir}/p3-region-records-diff.json"
python3 - "${passed}" "${failed}" "${skipped}" "${summary_json}" "${checks_json[@]}" <<'PY'
import json, pathlib, sys

passed, failed, skipped = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
summary_path = pathlib.Path(sys.argv[4])
checks = [json.loads(pathlib.Path(p).read_text(encoding="utf-8")) for p in sys.argv[5:]]
summary_path.write_text(
    json.dumps({"passed": passed, "failed": failed, "skipped": skipped, "checks": checks}, indent=2)
    + "\n",
    encoding="utf-8",
)
PY

if [[ "${failed}" -gt 0 ]]; then
  echo "P3 region records differential checks failed: ${failed}"
  exit 1
fi

echo "P3 region records differential checks passed: ${passed}, skipped: ${skipped}"

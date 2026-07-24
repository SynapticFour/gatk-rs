#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
cd "${repo_root}"

manifest="${repo_root}/parity/fixtures/p11_field_diff_cases.tsv"
if [[ ! -f "${manifest}" ]]; then
  echo "Missing corpus manifest: ${manifest}" >&2
  exit 2
fi

target_dir="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target-parity}"
mkdir -p "${target_dir}"
json_out="${report_dir}/p11_hc_output_field_diff_corpus.json"
tmp_dir="${report_dir}/p11-field-diff-corpus-tmp"
mkdir -p "${tmp_dir}"

python3 - "${manifest}" "${repo_root}" "${tmp_dir}" "${json_out}" "${target_dir}" <<'PY'
import csv
import json
import pathlib
import subprocess
import sys

manifest = pathlib.Path(sys.argv[1])
repo = pathlib.Path(sys.argv[2])
tmp_dir = pathlib.Path(sys.argv[3])
json_out = pathlib.Path(sys.argv[4])
target_dir = pathlib.Path(sys.argv[5])
docker_image = "us.gcr.io/broad-gatk/gatk:4.4.0.0"
docker_platform = "linux/amd64"

def run(cmd, env=None, check=True):
    res = subprocess.run(cmd, cwd=repo, env=env, capture_output=True, text=True)
    if check and res.returncode != 0:
        raise RuntimeError(f"command failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr}")
    return res

def count_variants(path: pathlib.Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for l in path.read_text(encoding="utf-8", errors="replace").splitlines() if l and not l.startswith("#"))

def first_variant(path: pathlib.Path):
    if not path.exists():
        return None
    for l in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not l or l.startswith("#"):
            continue
        cols = l.split("\t")
        if len(cols) < 8:
            return None
        out = {
            "CHROM": cols[0],
            "POS": cols[1],
            "REF": cols[3],
            "ALT": cols[4],
            "QUAL": cols[5],
            "FILTER": cols[6],
            "INFO": cols[7],
        }
        if len(cols) >= 10:
            out["FORMAT"] = cols[8]
            out["SAMPLE"] = cols[9]
        return out
    return None

def parse_info(s):
    d = {}
    if not s or s == ".":
        return d
    for tok in s.split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            d[k] = v
        else:
            d[tok] = "true"
    return d

def parse_sample(variant):
    if not variant:
        return {}
    fmt = variant.get("FORMAT")
    sample = variant.get("SAMPLE")
    if not fmt or not sample:
        return {}
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: vals[i] if i < len(vals) else "." for i, k in enumerate(keys)}

rows = list(csv.DictReader(manifest.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
if not rows:
    raise SystemExit("p11 corpus manifest is empty")

results = []
failed = 0
for row in rows:
    cid = row["case_id"]
    reference = repo / row["reference"]
    java_input = repo / row["java_input"]
    rust_input = row["rust_input"]
    interval = row["interval"]
    activate = row["activate_output"] == "1"
    expect_variant = row["expect_variant"] == "1"
    java_vcf = tmp_dir / f"{cid}.java.vcf"
    rust_vcf = tmp_dir / f"{cid}.rust.vcf"

    java_input_for_hc = java_input
    java_tmp_bam = None
    if java_input.suffix.lower() == ".sam":
        java_tmp_bam = tmp_dir / f"{cid}.java.bam"
        run([
            "docker", "run", "--rm", "--platform", docker_platform,
            "-v", f"{repo}:{repo}",
            "-w", str(repo),
            docker_image,
            "gatk", "SortSam",
            "-I", str(java_input),
            "-O", str(java_tmp_bam),
            "-SO", "coordinate",
            "--CREATE_INDEX", "true",
            "--QUIET", "true",
        ])
        java_input_for_hc = java_tmp_bam

    java_cmd = [
        "docker", "run", "--rm", "--platform", docker_platform,
        "-v", f"{repo}:{repo}",
        "-w", str(repo),
        docker_image,
        "gatk", "HaplotypeCaller",
        "-R", str(reference),
        "-I", str(java_input_for_hc),
        "-O", str(java_vcf),
        "--standard-min-confidence-threshold-for-calling", "0.0",
        "--verbosity", "ERROR",
    ]
    if interval and interval != "-":
        java_cmd.extend(["-L", interval])
    java_res = run(java_cmd, check=False)

    env = dict(**__import__("os").environ)
    env["CARGO_TARGET_DIR"] = str(target_dir)
    if activate:
        env.pop("GATK_RS_HC_SCAFFOLD_OUTPUT", None)
    else:
        env["GATK_RS_HC_SCAFFOLD_OUTPUT"] = "1"
    rust_cmd = [
        "cargo", "run", "--quiet", "--bin", "gatk-rs", "--",
        "HaplotypeCaller",
        "-R", row["reference"],
        "-I", rust_input,
        "-O", str(rust_vcf),
    ]
    if interval and interval != "-":
        rust_cmd.extend(["-L", interval])
    rust_res = run(rust_cmd, env=env, check=False)

    java_n = count_variants(java_vcf) if java_res.returncode == 0 else 0
    rust_n = count_variants(rust_vcf) if rust_res.returncode == 0 else 0
    java_first = first_variant(java_vcf)
    rust_first = first_variant(rust_vcf)
    status = "pass"
    mismatches = []

    if java_res.returncode != 0:
        status = "java_fail"
        mismatches.append("java_exit")
    elif rust_res.returncode != 0:
        status = "rust_fail"
        mismatches.append("rust_exit")
    elif expect_variant:
        if java_n == 0 or rust_n == 0:
            status = "fail"
            mismatches.append("variant_presence")
        else:
            for k in ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"]:
                if (java_first or {}).get(k) != (rust_first or {}).get(k):
                    mismatches.append(k)
            ji, ri = parse_info((java_first or {}).get("INFO", ".")), parse_info((rust_first or {}).get("INFO", "."))
            for k in ["AC", "AF", "AN", "DP"]:
                if k == "AF":
                    try:
                        if abs(float(ji.get(k, "nan")) - float(ri.get(k, "nan"))) > 0.01:
                            mismatches.append(f"INFO.{k}")
                    except Exception:
                        mismatches.append(f"INFO.{k}")
                elif ji.get(k) != ri.get(k):
                    mismatches.append(f"INFO.{k}")
            js, rs = parse_sample(java_first), parse_sample(rust_first)
            for k in ["GT", "AD", "DP", "GQ", "PL"]:
                if js.get(k) != rs.get(k):
                    mismatches.append(f"SAMPLE.{k}")
            if mismatches:
                status = "fail"
    else:
        if java_n != 0 or rust_n != 0:
            status = "fail"
            mismatches.append("expected_no_variant")

    if status != "pass":
        failed += 1
    results.append({
        "case_id": cid,
        "status": status,
        "expect_variant": expect_variant,
        "activate_output": activate,
        "java_exit": java_res.returncode,
        "rust_exit": rust_res.returncode,
        "java_variant_record_count": java_n,
        "rust_variant_record_count": rust_n,
        "mismatches": mismatches,
        "java_first_variant": java_first,
        "rust_first_variant": rust_first,
    })

summary = {
    "label": "phase11-hc-output-field-diff-corpus",
    "status": "pass" if failed == 0 else "fail",
    "failed_cases": failed,
    "total_cases": len(results),
    "cases": results,
}
json_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(f"[p11-corpus] status={summary['status']} failed={failed}/{len(results)}")
if failed:
    raise SystemExit(1)
PY

rm -f "${tmp_dir}"/*.java.vcf "${tmp_dir}"/*.rust.vcf "${tmp_dir}"/*.java.bam "${tmp_dir}"/*.java.bam.bai 2>/dev/null || true

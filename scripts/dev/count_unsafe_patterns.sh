#!/usr/bin/env bash
# Ratchet: production unwrap()/expect()/panic!() and .clone() counts must not grow.
#
# Scopes: gatk-haplotypecaller/src, gatk-core/src, gatk-cli/src
# Excludes:
#   - any path segment named tests/
#   - #[cfg(test)] modules inside .rs files (same policy as check_io_unwrap_policy.py)
#   - full-line // comments
#
# Usage:
#   ./scripts/dev/count_unsafe_patterns.sh           # check both ratchets
#   ./scripts/dev/count_unsafe_patterns.sh unwrap     # unwrap only
#   ./scripts/dev/count_unsafe_patterns.sh clone      # clone only
#   ./scripts/dev/count_unsafe_patterns.sh print      # print counts, never fail
#
# Baselines: .quality-gates/unwrap_baseline.txt, .quality-gates/clone_baseline.txt
# Raising a baseline requires an intentional commit message containing:
#   baseline-bump: <reason>
# (enforced in CI when comparing to the PR base branch).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

mode="${1:-check}"
gates_dir="${repo_root}/.quality-gates"
unwrap_baseline_file="${gates_dir}/unwrap_baseline.txt"
clone_baseline_file="${gates_dir}/clone_baseline.txt"

count_json="$(python3 - <<'PY'
from __future__ import annotations

import json
import re
from pathlib import Path

ROOTS = [
    Path("gatk-haplotypecaller/src"),
    Path("gatk-core/src"),
    Path("gatk-cli/src"),
]
HIT_UNWRAP = re.compile(r"\.unwrap\(\)|\.expect\(|panic!\(")
HIT_CLONE = re.compile(r"\.clone\(\)")


def test_module_ranges(lines: list[str]) -> list[tuple[int, int]]:
    ranges: list[tuple[int, int]] = []
    i = 0
    while i < len(lines):
        if re.search(r"#\[cfg\(test\)\]", lines[i]):
            j = i + 1
            while j < len(lines) and (
                lines[j].strip() == ""
                or lines[j].strip().startswith("#[")
                or lines[j].strip().startswith("//")
            ):
                j += 1
            if j < len(lines) and re.match(r"\s*(pub\s+)?mod\s+\w+", lines[j]):
                k = j
                while k < len(lines) and "{" not in lines[k]:
                    k += 1
                if k < len(lines):
                    bal = 0
                    for t in range(k, len(lines)):
                        bal += lines[t].count("{") - lines[t].count("}")
                        if bal == 0:
                            ranges.append((i + 1, t + 1))
                            break
        i += 1
    return ranges


def in_test(n: int, ranges: list[tuple[int, int]]) -> bool:
    return any(a <= n <= b for a, b in ranges)


unwrap = 0
clone = 0
unwrap_by_file: dict[str, int] = {}
clone_by_file: dict[str, int] = {}

for root in ROOTS:
    if not root.is_dir():
        continue
    for path in sorted(root.rglob("*.rs")):
        if "tests" in path.parts:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        ranges = test_module_ranges(lines)
        u = c = 0
        for n, line in enumerate(lines, 1):
            if in_test(n, ranges):
                continue
            stripped = line.strip()
            if stripped.startswith("//"):
                continue
            if HIT_UNWRAP.search(line):
                u += 1
            if HIT_CLONE.search(line):
                c += 1
        rel = str(path)
        unwrap += u
        clone += c
        if u:
            unwrap_by_file[rel] = u
        if c:
            clone_by_file[rel] = c

print(
    json.dumps(
        {
            "unwrap": unwrap,
            "clone": clone,
            "unwrap_by_file": unwrap_by_file,
            "clone_by_file": clone_by_file,
        }
    )
)
PY
)"

unwrap_count="$(python3 -c 'import json,sys; print(json.loads(sys.argv[1])["unwrap"])' "${count_json}")"
clone_count="$(python3 -c 'import json,sys; print(json.loads(sys.argv[1])["clone"])' "${count_json}")"

read_baseline() {
  local file="$1"
  local name="$2"
  if [[ ! -f "${file}" ]]; then
    echo "[unsafe-patterns] missing baseline ${file}" >&2
    exit 2
  fi
  local raw
  raw="$(tr -d '[:space:]' <"${file}")"
  if [[ ! "${raw}" =~ ^[0-9]+$ ]]; then
    echo "[unsafe-patterns] ${name} baseline must be a non-negative integer (got: ${raw})" >&2
    exit 2
  fi
  printf '%s\n' "${raw}"
}

print_top() {
  local kind="$1"
  python3 -c '
import json,sys
d=json.loads(sys.argv[1])
key=sys.argv[2]
items=sorted(d[key].items(), key=lambda kv: (-kv[1], kv[0]))[:12]
for path,n in items:
    print(f"  {n:4d}  {path}")
' "${count_json}" "${kind}"
}

echo "[unsafe-patterns] production unwrap/expect/panic!: ${unwrap_count}"
echo "[unsafe-patterns] production .clone(): ${clone_count}"
echo "[unsafe-patterns] top unwrap files:"
print_top unwrap_by_file
echo "[unsafe-patterns] top clone files:"
print_top clone_by_file

if [[ "${mode}" == "print" ]]; then
  exit 0
fi

fail=0
check_one() {
  local kind="$1"
  local count="$2"
  local baseline_file="$3"
  if [[ "${mode}" != "check" && "${mode}" != "${kind}" ]]; then
    return 0
  fi
  local baseline
  baseline="$(read_baseline "${baseline_file}" "${kind}")"
  echo "[unsafe-patterns] ${kind}: count=${count} baseline=${baseline}"
  if (( count > baseline )); then
    echo "[unsafe-patterns] FAIL: ${kind} count ${count} exceeds baseline ${baseline}" >&2
    echo "[unsafe-patterns] Fix new sites, or raise ${baseline_file} in the same PR with commit message:" >&2
    echo "[unsafe-patterns]   baseline-bump: <why the ratchet must move>" >&2
    fail=1
  elif (( count < baseline )); then
    echo "[unsafe-patterns] note: ${kind} count ${count} is below baseline ${baseline} — consider lowering the baseline in this PR"
  fi
}

check_one unwrap "${unwrap_count}" "${unwrap_baseline_file}"
check_one clone "${clone_count}" "${clone_baseline_file}"

# When a baseline file itself was raised vs the merge base, demand an explicit commit rationale.
if [[ -n "${GITHUB_BASE_REF:-}" ]] || [[ -n "${UNSAFE_PATTERNS_BASE_REF:-}" ]]; then
  base_ref="${UNSAFE_PATTERNS_BASE_REF:-origin/${GITHUB_BASE_REF}}"
  check_baseline_bump() {
    local file="$1"
    local label="$2"
    if ! git cat-file -e "${base_ref}:${file}" 2>/dev/null; then
      return 0
    fi
    local old new
    old="$(git show "${base_ref}:${file}" | tr -d '[:space:]')"
    new="$(tr -d '[:space:]' <"${file}")"
    if [[ "${old}" =~ ^[0-9]+$ && "${new}" =~ ^[0-9]+$ ]] && (( new > old )); then
      if ! git log "${base_ref}..HEAD" --format=%B | grep -qiE 'baseline-bump:'; then
        echo "[unsafe-patterns] FAIL: ${label} baseline raised ${old} → ${new} without 'baseline-bump:' in a commit message" >&2
        fail=1
      else
        echo "[unsafe-patterns] ok: ${label} baseline bump ${old} → ${new} documented via baseline-bump:"
      fi
    fi
  }
  if [[ "${mode}" == "check" || "${mode}" == "unwrap" ]]; then
    check_baseline_bump ".quality-gates/unwrap_baseline.txt" "unwrap"
  fi
  if [[ "${mode}" == "check" || "${mode}" == "clone" ]]; then
    check_baseline_bump ".quality-gates/clone_baseline.txt" "clone"
  fi
fi

exit "${fail}"

# Quality-gate ratchets

Numeric ceilings checked by `scripts/dev/count_unsafe_patterns.sh` on every PR.

| File | Metric | Scope |
|------|--------|--------|
| `unwrap_baseline.txt` | `.unwrap()` / `.expect(` / `panic!(` | Production code under `gatk-haplotypecaller/src`, `gatk-core/src`, `gatk-cli/src` |
| `clone_baseline.txt` | `.clone()` | Same scope |

**Excluded:** `tests/` path segments and `#[cfg(test)]` modules (test harness may panic/assert freely).

## Policy

- Counts **must not exceed** the baseline (CI fails the PR).
- Prefer fixing new sites. If a raise is unavoidable, bump the baseline file **in the same PR** and include a commit message line:

  ```
  baseline-bump: <short reason>
  ```

  Silent baseline increases are rejected when the PR base branch is available.
- When a cleanup lowers the count, **lower the baseline in the same PR** so the ratchet stays tight.

## Local

```bash
./scripts/dev/count_unsafe_patterns.sh          # check both
./scripts/dev/count_unsafe_patterns.sh print    # counts only
./scripts/dev/count_unsafe_patterns.sh unwrap
./scripts/dev/count_unsafe_patterns.sh clone
```

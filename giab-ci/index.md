# GIAB genome-wide equivalence dashboard

- **Run dir:** `/home/runner/work/gatk-rs/gatk-rs/parity/giab/runs/ci`
- **Generated (UTC):** 2026-08-15T10:12:00Z
- **Scope:** CI-SUBSET (default “genome-wide” in this repo): FULL chr20 + FULL chr21 + one 50kb probe on each other autosome. Not all bases of chr1–19/22. Omits Peak-hang shards 00_chr20_w47 + 01_chr21_w10 unless GIAB_INCLUDE_HANG_SHARDS=1 (docs/perf/CI_SUBSET_HANG_W47_W10.md).

## Samples

| Sample | Mode | Gate | Max\|ΔF1\| | Java wall | Java RSS | Rust wall | Rust RSS |
|--------|------|:----:|----------:|----------:|---------:|----------:|---------:|
| HG001 | ci-subset | FAIL | 0.0246 | — | — | — | — |

## Notes

- Primary equivalence metric is **Rust−Java F1 delta** via `gatk-rs-equiv` (hap.py / RTG), not absolute F1.
- Wall / RSS from `/usr/bin/time -v` (Linux) or `/usr/bin/time -l` (macOS).
- “Genome-wide” follows `GIAB_MODE` — see `SCOPE.txt` in the run directory.


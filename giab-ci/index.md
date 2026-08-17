# GIAB genome-wide equivalence dashboard

- **Run dir:** `/home/runner/work/gatk-rs/gatk-rs/parity/giab/runs/ci`
- **Generated (UTC):** 2026-08-17T07:47:23Z
- **Scope:** WALL-LOSERS: eight 1 Mb dense campaign windows (chr20/21 w09/w11/w26/w29). Product wall (no GATK_RS_HC_SEQUENTIAL). Peak abort retained.

## Samples

| Sample | Mode | Gate | Max\|ΔF1\| | Java wall | Java RSS | Rust wall | Rust RSS |
|--------|------|:----:|----------:|----------:|---------:|----------:|---------:|
| HG001 | wall-losers | PASS | 0.0004 | — | — | — | — |

## Notes

- Primary equivalence metric is **Rust−Java F1 delta** via `gatk-rs-equiv` (hap.py / RTG), not absolute F1.
- Wall / RSS from `/usr/bin/time -v` (Linux) or `/usr/bin/time -l` (macOS).
- “Genome-wide” follows `GIAB_MODE` — see `SCOPE.txt` in the run directory.


# GIAB genome-wide equivalence dashboard

- **Run dir:** `/home/runner/work/gatk-rs/gatk-rs/parity/giab/runs/ci`
- **Generated (UTC):** 2026-08-01T10:08:49Z
- **Scope:** SMOKE: three ~50kb windows (chr20/chr21/P12). P12 reads from NA12878_20k evidence class; chr20/21 from HG001 30×. Full-30× P12 is benchmark-host only. Not genome-wide.

## Samples

| Sample | Mode | Gate | Max\|ΔF1\| | Java wall | Java RSS | Rust wall | Rust RSS |
|--------|------|:----:|----------:|----------:|---------:|----------:|---------:|
| HG001 | smoke | FAIL | — | — | — | — | — |

## Notes

- Primary equivalence metric is **Rust−Java F1 delta** via `gatk-rs-equiv` (hap.py / RTG), not absolute F1.
- Wall / RSS from `/usr/bin/time -v` (Linux) or `/usr/bin/time -l` (macOS).
- “Genome-wide” follows `GIAB_MODE` — see `SCOPE.txt` in the run directory.


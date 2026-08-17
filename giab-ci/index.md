# GIAB genome-wide equivalence dashboard

- **Run dir:** `/home/runner/work/gatk-rs/gatk-rs/parity/giab/runs/ci`
- **Generated (UTC):** 2026-08-17T21:39:03Z
- **Scope:** CI-SUBSET (default “genome-wide” in this repo): FULL chr20 + FULL chr21 + one 50kb probe on each other autosome. Not all bases of chr1–19/22.

## Samples

| Sample | Mode | Gate | Max\|ΔF1\| | Java wall | Java RSS | Rust wall | Rust RSS |
|--------|------|:----:|----------:|----------:|---------:|----------:|---------:|
| HG001 | ci-subset | PASS | 0.0099 | — | — | — | — |

## Notes

- Primary equivalence metric is **Rust−Java F1 delta** via `gatk-rs-equiv` (hap.py / RTG), not absolute F1.
- Wall / RSS from `/usr/bin/time -v` (Linux) or `/usr/bin/time -l` (macOS).
- “Genome-wide” follows `GIAB_MODE` — see `SCOPE.txt` in the run directory.


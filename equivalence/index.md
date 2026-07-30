# Equivalence dashboard (nightly trio E2E)

Full spine: **HaplotypeCaller (GVCF) → CombineGVCFs → GenotypeGVCFs → VariantFiltration**
for GIAB Ashkenazi trio HG002/HG003/HG004, scored with **hap.py** vs HG002 truth.

- **Generated (UTC):** 2026-07-30T11:54:17Z
- **Commit:** `23d524b064ff6101f2bfe50e7a32e98e90926414`
- **Run dir:** `/home/runner/work/gatk-rs/gatk-rs/parity/giab/runs/nightly_trio_20260730T060719Z`
- **Baseline:** `docs/equivalence/baseline.json`
- **Regression threshold (|ΔF1| drop):** 0.02
- **Regressions:** 0

## Regions

| Region | Kind | Status | Rust SNP F1 | Rust INDEL F1 | Java SNP F1 | Java INDEL F1 |
|--------|------|--------|------------:|--------------:|------------:|--------------:|
| chr20 | chromosome | running | — | — | — | — |

## Regressions vs last green

_None above threshold._

## Notes

- BAMs are **region-sliced** with `samtools view -L` (no full WGS download).
- Hard regions are capped slices of GIAB stratification BEDs (segdups / TR / alldifficult / MHC).
- Soft gate: regressions open a GitHub issue (`equivalence-regression`); the workflow does not hard-fail.


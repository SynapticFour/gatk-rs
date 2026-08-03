# Equivalence dashboard (nightly trio E2E)

Full spine: **HaplotypeCaller (GVCF) → CombineGVCFs → GenotypeGVCFs → VariantFiltration**
for GIAB Ashkenazi trio HG002/HG003/HG004, scored with **hap.py** vs HG002 truth.

- **Generated (UTC):** 2026-08-03T12:03:13Z
- **Commit:** `50e95d679cb8fd6c139b99ffbe373e6814efb0cb`
- **Run dir:** `/home/runner/work/gatk-rs/gatk-rs/parity/giab/runs/nightly_trio_20260803T061551Z`
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


# Equivalence dashboard (nightly trio E2E)

Full spine: **HaplotypeCaller (GVCF) → CombineGVCFs → GenotypeGVCFs → VariantFiltration**
for GIAB Ashkenazi trio HG002/HG003/HG004, scored with **hap.py** vs HG002 truth.

- **Generated (UTC):** 2026-08-04T11:57:44Z
- **Commit:** `0a450d3059bc81a543d6b32e66676a79c4b8766c`
- **Run dir:** `/home/runner/work/gatk-rs/gatk-rs/parity/giab/runs/nightly_trio_20260804T061124Z`
- **Baseline:** `docs/equivalence/baseline.json`
- **Regression threshold (|ΔF1| drop):** 0.02
- **Regressions:** 0

## Regions

| Region | Kind | Status | Rust SNP F1 | Rust INDEL F1 | Java SNP F1 | Java INDEL F1 |
|--------|------|--------|------------:|--------------:|------------:|--------------:|
| chr20 | chromosome | happy_failed | — | — | — | — |
| chr21 | chromosome | running | — | — | — | — |

## Regressions vs last green

_None above threshold._

## Notes

- BAMs are **region-sliced** with `samtools view -L` (no full WGS download).
- Hard regions are capped slices of GIAB stratification BEDs (segdups / TR / alldifficult / MHC).
- Soft gate: regressions open a GitHub issue (`equivalence-regression`); the workflow does not hard-fail.


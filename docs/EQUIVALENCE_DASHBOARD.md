# Equivalence dashboard (nightly trio E2E)

**Not a signed product claim.** Until the table below leaves “pending”, no
hap.py number from this page is evidence. Authority remains
[`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). The nightly publisher
(`scripts/parity/giab/publish_nightly_equivalence.py`, invoked by
`.github/workflows/nightly-equivalence.yml`) replaces this file when it runs.

This file is **updated automatically** by `.github/workflows/nightly-equivalence.yml`.

Full spine: **HaplotypeCaller (GVCF) → CombineGVCFs → GenotypeGVCFs → VariantFiltration**
for GIAB Ashkenazi trio HG002/HG003/HG004, scored with **hap.py** vs HG002 truth.

- **Generated (UTC):** _(pending first nightly run)_
- **Commit:** `—`
- **Baseline:** `docs/equivalence/baseline.json`
- **Regression threshold (|ΔF1| drop):** `0.02` (configurable via workflow input)

## Regions

| Region | Kind | Status | Rust SNP F1 | Rust INDEL F1 | Java SNP F1 | Java INDEL F1 |
|--------|------|--------|------------:|--------------:|------------:|--------------:|
| _(pending)_ | | | | | | |

## Notes

- BAMs are region-sliced with `samtools view -L` (no full WGS download).
- Hard regions are capped slices of GIAB stratification BEDs.
- Soft gate: F1 regressions open a GitHub issue labeled `equivalence-regression`.
- HTML mirror (when Pages is enabled): `/equivalence/` on the project GitHub Pages site.

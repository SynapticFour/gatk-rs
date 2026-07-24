# J6 — NA12878 scale + GIAB truth (L6)

**Layer:** L6 (scale + external truth)  
**Gate script:** `scripts/parity/run_hc_full_parity_j6_truth.sh`  
**Sign-off:** [`docs/CLAIM_MATRIX.md`](../../../../docs/CLAIM_MATRIX.md)

## What this gate proves

1. **Scale:** Java and Rust `HaplotypeCaller` both complete on NA12878_20k_b37 over the **L6 GIAB window** (`2:92000000-92400000` by default) with **variant-set parity** (`parity_status=variant_parity`). This window **includes** the P12 spine (`92300000-92350000`); the spine alone has **zero** GIAB high-confidence truth sites.
2. **Truth (P13):** Both callsets are compared to GIAB HG001 v4.2.1 (GRCh37) inside high-confidence regions, with **overall and stratified (SNP / INDEL) F1** metrics.
3. **Gate:** Rust F1 must track Java F1 within documented bounds (`thresholds.json`).

## Prerequisites

| Asset | Default path |
|-------|----------------|
| Reference (hs37d5 simple) | `parity/realworld/assets/hs37d5.simple.fa` |
| GIAB truth VCF | `parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz` |
| GIAB high-confidence BED | `parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.bed` |
| NA12878 20k BAM | `parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam` (auto-download) |

Stage reference + truth:

```bash
./scripts/parity/realworld/03_stage_reference_and_truth.sh
# REALWORLD_STOP_AFTER_ASSETS=1 exits after download only
```

## Run

```bash
export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
./scripts/parity/run_hc_full_parity_j6_truth.sh
```

Reports: `parity/reports/hc-full-parity-j6/j6_truth_summary.{json,md}`

## Thresholds

See [`thresholds.json`](./thresholds.json):

- Overall: `rust_f1 >= 0.95 * java_f1` and `rust_f1 >= java_f1 - 0.05`
- SNP stratum: same as overall
- INDEL stratum: `0.90 * java_f1`, max delta `0.08` (sparse indels on 20k slice)

Adjust thresholds only with sign-off doc update + canonical log refresh.

## Dense GIAB window (R3 / L7)

```bash
export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
J6_DENSE=1 ./scripts/parity/run_hc_full_parity_j6_truth.sh
```

Reports: `parity/reports/hc-full-parity-j6-dense/`  
Thresholds: [`thresholds_dense.json`](./thresholds_dense.json)

### L7 second non-chr2 slice (chr21)

Uses a separate BAM OUT_DIR + report dir (`J6_REPORT_DIR`) so chr20 artifacts are not clobbered:

```bash
./scripts/parity/run_hc_full_parity_j6_dense_chr21.sh
```

Thresholds: [`thresholds_dense_chr21.json`](./thresholds_dense_chr21.json)  
Default interval: `21:41200001-41250000` (~200 HC-BED truth sites)

Optional overrides: `J6_DENSE_INTERVAL`, `J6_DENSE_OUT_DIR`, `J6_REPORT_DIR`, `J6_THRESHOLDS`.

## CI

- Weekly spine: `.github/workflows/p12-l6-scale.yml` (Sunday 04:00 UTC, `workflow_dispatch`)
- L7 dense + spine: `.github/workflows/p12-l7-dense-spine.yml` (Sunday 05:00 UTC, `workflow_dispatch`, path-filtered PR)

Manual full chr2 exploration (not CI default):

```bash
J6_INTERVAL=2:92000000-92400000 ./scripts/parity/run_hc_full_parity_j6_truth.sh
```

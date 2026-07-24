# Equivalence proof surface

Runnable evidence for claims in [`docs/CLAIM_MATRIX.md`](../../docs/CLAIM_MATRIX.md).

| Location | Role |
|----------|------|
| [`gatk-rs-equiv/`](../../gatk-rs-equiv/) | GIAB / hap.py / vcfeval CLI + differential fuzz driver |
| [`fuzz/`](../../fuzz/) | LibFuzzer target (`hc_differential`) sharing scenarios with equiv |
| [`scripts/parity/`](../../scripts/parity/) | L2 synthetic gates, P12 NA12878 windows, GIAB staging (`giab/`) |
| [`scripts/parity/run_combine_gvcfs_parity.sh`](../../scripts/parity/run_combine_gvcfs_parity.sh) | Mini-cohort CombineGVCFs Java vs Rust (REF/ALT/PL) |
| [`scripts/parity/run_genotype_gvcfs_parity.sh`](../../scripts/parity/run_genotype_gvcfs_parity.sh) | Mini-cohort Combine→GenotypeGVCFs Java vs Rust |
| [`scripts/parity/giab/run_trio_joint_genotyping_e2e.sh`](../../scripts/parity/giab/run_trio_joint_genotyping_e2e.sh) | E2E HC(gVCF)→Combine→Genotype (smoke or GIAB HG002/3/4) |
| [`scripts/parity/run_variant_filtration_parity.sh`](../../scripts/parity/run_variant_filtration_parity.sh) | Hard-filter boundary FILTER decisions vs Java VariantFiltration |
| [`parity/combine_gvcfs/mini/`](../../parity/combine_gvcfs/mini/) | Synthetic two-sample gVCF + tiny FASTA for Combine/Genotype |
| [`parity/fixtures/`](../../parity/fixtures/) | Tracked fixtures consumed by those scripts |
| [`docs/GATK_PINNED.env`](../../docs/GATK_PINNED.env) | Pinned GATK 4.4.0.0 oracle coordinates |

Historical L2–L14 sign-off writeups live on git branch `pre-cleanup-archive`.

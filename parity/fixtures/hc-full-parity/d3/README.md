# Phase D.3 — HQ soft-clip mean

Gate: `scripts/parity/run_hc_full_parity_d3_soft_clip.sh`

Producer: `soft-clip-mean <ref.fa> <sam|bam> <interval_cli>`

TSV column `hq_soft_clip_mean` is the **RCM** `MathUtils.RunningAverage` mean from `calcGenotypeLikelihoodsOfRefVsAny` (same update rule as [`hq_soft_clip_running_mean_rcm_path`](../../../../gatk-haplotypecaller/src/activity_scoring.rs) / GATK `AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD`).

Fixtures used for Java L2 (`HcFullParityGateDump soft-clip-mean`) should include an **`@RG` with `SM:`** so `HaplotypeCallerEngine` can initialize samples; Rust does not require this.

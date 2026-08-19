# PairHMM speedup profile (reproducible)

**Generated (UTC):** `20260724T055240Z`  
**Host:** `Darwin 25.5.0 arm64`  
**Runner:** [`scripts/perf/run_pairhmm_speedup.sh`](../../scripts/perf/run_pairhmm_speedup.sh)  
**Raw run:** `docs/perf/runs/pairhmm_20260724T055240Z/`

## Build / versions

- rustc: `rustc 1.88.0 (6b00bc388 2025-06-23)`
- cargo: `cargo 1.88.0 (873a06493 2025-05-10)`
- git: `ecb97e5`
- build: `cargo build -p gatk-cli --release --locked`
- features detected (coarse): `['neon']`

## Microbench (Criterion `pairhmm_logless_simd`)

Read length 200 bp, parity indel/GCP quals 45/45/10, baseQ 30.

| hap_count | scalar Logless | SIMD | speedup (scalar/SIMD) |
|---|---:|---:|---:|
| 8 | 1.254 ms | 2.334 ms | **0.54×** |
| 32 | 5.069 ms | 3.625 ms | **1.40×** |
| 64 | 10.165 ms | 7.219 ms | **1.41×** |


Exact command:

```bash
cargo bench -p gatk-haplotypecaller --bench pairhmm --locked -- pairhmm_logless_simd
```

## Production default

**Code default is `FASTEST_AVAILABLE`** (`HcLikelihoodEngineConfig` /
`pair_hmm_impl_from_config` when unset → host SIMD when present, else Logless).

**Signed promotion** (claim matrix / GIAB ci-subset F1) is still required before
marketing SIMD-as-oracle equivalence. Unit gate remains mandatory on every PairHMM
change:

`cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test`

Architecture C (opt-in wavefront, not production default):
[`PAIRHMM_WAVEFRONT.md`](PAIRHMM_WAVEFRONT.md).

Criterion phenotype matrix (read 100/200/300 × hap 8/32/64/128):  
`cargo bench -p gatk-haplotypecaller --bench pairhmm --locked -- pairhmm_logless_simd`

Smith-Waterman TLS scratch is **bounded** and contig-scale matrices are **refused**
(`smith_waterman::oversized_matrix_is_refused`) — Peak-RSS safety, not a downsampler.

## Phase B note (2026-08-05)

Packed Logless reuses DP scratch across haplotypes; k-best borrows the graph when
cycle-stripping is off. Default PairHMM remains Log10.

## HC smoke (fixture)

```bash
./target/release/gatk-rs HaplotypeCaller -R parity/fixtures/reference.fa \
  -I parity/fixtures/sample.bam -O /tmp/hc.vcf -L chr1:1-32 --pair-hmm LOG10_PAIRHMM
./target/release/gatk-rs HaplotypeCaller ... --pair-hmm SIMD
```

See `docs/perf/runs/pairhmm_20260724T055240Z/hc_*.time` for wall/RSS from `/usr/bin/time -l`.

## GIAB / default promotion

GIAB `ci-subset` equivalence is **not signed** in [`CLAIM_MATRIX.md`](../CLAIM_MATRIX.md).
This run validated SIMD vs scalar Logless unit tests + HC fixture smoke only.
Default is already `FASTEST_AVAILABLE` in tree; keep the **signed hap.py / L9 F1**
gate green before claiming GIAB equivalence. Do not silently change resolve order
(SIMD → Logless → Log10) without rematch + F1.

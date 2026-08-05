# HaplotypeCaller memory profile (reproducible)

**Generated (UTC):** `20260804T120757Z`  
**Host:** `Darwin 25.5.0 arm64`  
**Git:** `7347f88`  
**Runner script:** [`scripts/perf/run_hc_memory_profile.sh`](../../scripts/perf/run_hc_memory_profile.sh)  
**Raw run directory:** `docs/perf/runs/20260804T120757Z/`

**Public memory claim status:** **not allowed** from this run (need realistic profile on dedicated `gatk-rs-benchmark` host).

Profiles measured: `smoke`.

## Memory footprint pass (2026-08-04)

Engineering on branch `perf/giab-memory-footprint` (Air M4 16 GiB gold standard):

| Change | Where |
|--------|--------|
| `SharedBamRecord` (`Arc`) — no deep clone on region fill; progressive shard release | [`shared_bam.rs`](../../gatk-haplotypecaller/src/shared_bam.rs), [`assembly_region_iterator.rs`](../../gatk-haplotypecaller/src/assembly_region_iterator.rs) |
| Sequential regions / stream-merge VCF batches | [`run.rs`](../../gatk-haplotypecaller/src/run.rs) (`GATK_RS_HC_SEQUENTIAL`, lower large-region threshold) |
| PairHMM / SW TLS scratch shrink | `pairhmm_log10` / `pairhmm_logless` / `smith_waterman` |
| Arc reference windows | [`reference_context.rs`](../../gatk-haplotypecaller/src/reference_context.rs) |
| Finalize disk hygiene (hardlink / delete shards+BAM) | [`run_genomewide_equivalence.sh`](../../scripts/parity/giab/run_genomewide_equivalence.sh) |

Local recipe: [`AIR_M4_GIAB_RECIPE.md`](AIR_M4_GIAB_RECIPE.md). Pre-change anchors: [`HC_MEMORY_BASELINE_20260804.md`](HC_MEMORY_BASELINE_20260804.md).

**Do not** re-run 2 Mb realistic / full 30× P12 on this 16 GiB host until dense-window RSS is re-measured with `GATK_RS_HC_SEQUENTIAL=1`. **Do not** dispatch signed `ci-subset` until that proof exists.

## Rust showcase Phase A (ownership) — 2026-08-05

Roadmap: [`RUST_SHOWCASE_ROADMAP.md`](RUST_SHOWCASE_ROADMAP.md).

| Change | Where |
|--------|--------|
| Byte-native assembly / kmers (`Vec<u8>`) | `assembly.rs`, `read_threading_graph.rs`, … |
| Single finalize buffer for assemble + PairHMM | `assembly_based_caller.rs`, `engine.rs` |
| COW-aware unique ownership before realign | `shared_bam::into_unique_records`, `engine.rs` |
| Sequential hap scoring under `GATK_RS_HC_SEQUENTIAL=1` | `likelihood_engine.rs` |
| Optional jemalloc | `gatk-cli` feature `jemalloc` + `MALLOC_CONF` notes in Air recipe |

**Measured this land (no realistic BAM staged):** `check_hc_rss_regression.sh` OK — dict Peak-RSS ≈ 6.0 MiB; bomb/50 kb/100 kb HC Peak-RSS **not** re-recorded (assets missing). Excellence N-7 band freeze PASS. No public Peak-RSS claim from smoke.

## A. Trivial smoke — reproducibility reference only

**Label:** Trivial smoke (reproducibility only)  
**Interval:** `chr1:1-32`  
**BAM / ref:** checked-in `parity/fixtures/`  

> **Not for marketing.** Peak-RSS here is dominated by JVM/runtime fixed cost
> on a 32 bp window. Do **not** derive “X% less memory” from this table.

| Engine | Peak RSS | Wall time |
|--------|----------|-----------|
| **gatk-rs** (Rust release) | **9.48 MiB (9712 KiB)** | 0.63 s |
| **Java GATK 4.4.0.0** | **430.66 MiB (441000 KiB)** | 4.908 s |

| Java / Rust Peak-RSS | 45.41× |
| Rust as fraction of Java Peak-RSS | 2.2% |
| Absolute delta (Java − Rust) | 421.18 MiB |


## Exact commands

### Rust (build once)

```bash
cargo build -p gatk-cli --release --locked
# rustc: rustc 1.95.0 (59807616e 2026-04-14) (Homebrew)
# cargo: cargo 1.95.0 (f2d3ce0bd 2026-03-21) (Homebrew)
# git: 7347f88
```

Per-profile command lines are in each `docs/perf/runs/20260804T120757Z/<profile>/` summary
and the JSON under `rust.cmdline` / `java.cmdline.txt`.

### Java GATK 4.4

- Pin: `GATK_PINNED_SHA=2dbc025821bc5f686c423ff332a41e6cef892a77` (`docs/GATK_PINNED.env`)
- Image: `us.gcr.io/broad-gatk/gatk:4.4.0.0`
- JVM options (pipeline-realistic): `-Xms1g -Xmx4g`

```bash
./scripts/perf/run_hc_memory_profile.sh
# smoke only:
#   HC_MEM_PROFILES=smoke ./scripts/perf/run_hc_memory_profile.sh
# realistic only (stages 2 Mb GIAB window if needed):
#   HC_MEM_PROFILES=realistic ./scripts/perf/run_hc_memory_profile.sh
```

When Docker is used, Peak-RSS is sampled from `/proc/*/status` **VmHWM**
for `java`/`gatk` **inside** the Linux container (the Broad 4.4 image has no
GNU `/usr/bin/time`). Host `time docker …` is never used for RSS.

## Re-run

```bash
./scripts/perf/run_hc_memory_profile.sh
# overrides:
#   HC_MEM_PROFILES=smoke,realistic
#   HC_MEM_REALISTIC_INTERVAL=20:10000000-12000000
#   JAVA_XMX=4g JAVA_XMS=1g
# Dedicated host (published claims):
#   see docs/ci/PERF_BENCHMARK_HOST.md + Actions workflow benchmark.yml
```

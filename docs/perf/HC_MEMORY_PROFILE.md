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

**Measured 2026-08-05** (Air M4 Darwin arm64, `GATK_RS_HC_SEQUENTIAL=1`, `RAYON_NUM_THREADS=1`, release `gatk-rs` without jemalloc; BAM `parity/realworld/na12878_giab_window_mem_500kb_b37`, ref `hs37d5.simple.fa`):

| Window | Interval | Peak-RSS | Exit |
|--------|----------|----------|------|
| bomb | `20:10098500-10099500` | **26.86 MiB** | 0 |
| 50 kb | `20:10000000-10050000` | **31.50 MiB** | 0 |
| 100 kb | `20:10000000-10100000` | **2632 MiB** then abort | 1 |
| 500 kb | `20:10000000-10500000` | skipped after 100 kb failure | — |

Pre-pass anchors (`HC_MEMORY_BASELINE_20260804.md`): bomb ~38 MiB, 50 kb ~114 MiB. **Dense-window Peak-RSS drop is real on bomb/50 kb**; 100 kb still fails on this 16 GiB host — **not** a public claim and **not** a GIAB `ci-subset` sign. Excellence N-7 band freeze PASS.

### 100 kb follow-up (`perf/fix-100kb-rss`, 2026-08-05 evening)

Ownership + DP fail-closed land (clip-in-place, SharedBam realign, 8 M PairHMM/SW cell cap, stream SeqGraphs, sequential BAM header share). Resource-limited re-measure (`nice -n 15`, `CARGO_BUILD_JOBS=1`):

| Window | Peak-RSS | Exit | Notes |
|--------|----------|------|-------|
| bomb | 38.97 MiB | 0 | OK |
| 50 kb | 48.58 MiB | 0 | OK (higher than morning run; still ≪ pre-pass 114 MiB) |
| 100 kb | **~136 MiB live RSS** after 10 min | stopped | No longer climbs to ~2.6 GiB; CPU-bound on 1001 tiles — do **not** finish on a multitasking 16 GiB Air; re-run overnight / dedicated host before signing GIAB |

**Overnight 2026-08-06:** sequential 100 kb reached Peak-RSS **~825 MiB** then jetsam (`terminated abnormally`, load spike ~282) under ~8 GiB compressor — not a return of the 2.6 GiB climb. Soft TLS keep left PairHMM/SW arenas sticky across phases.

**Follow-up (TLS):** hard-clear PairHMM/SW TLS after regions; release SW before PairHMM and PairHMM before realign; Rayon `broadcast` TLS clear after parallel batches. Bomb/50 kb: **24.6 / 40.4 MiB** exit 0. **jemalloc on Darwin raised Peak** (attempt Peak ~5.4 GiB) — prefer system allocator for Peak gates.

**Later overnight (same day):** three sequential attempts still peaked **2.2–4.2 GiB** with a **sawtooth** live RSS (~200↔1600 MiB), then jetsam. 10 k region-read refuse did not fire → spikes are not “>10 k reads alone.”

### Theory (2026-08-06)

Industry context: Broad HC Peak targets ≤3 GB typical / ≤6 GB worst ([gatk#2591](https://github.com/broadinstitute/gatk/issues/2591)); peaks from high-depth/repetitive pileups and bushy haplotype graphs. Spark notes lazy region materialization ([gatk#4376](https://github.com/broadinstitute/gatk/issues/4376)). Rust DP TLS high-water is a known Peak pin ([rammap#9](https://github.com/jwanglab/rammap/pull/9)).

**Our Peak formula (sequential):**  
`Peak ≈ live_shard_BAM + region_Arcs + finalize_owned + assembly_working_set + max(PairHMM_TLS, SW_TLS) + allocator_dirt`

Layers: (A) eager whole-`-L` shard load + dual residency with `region.reads`; (B) finalize owned copies; (C) assembly/kbest; (D) capped TLS (~256 MiB) — not the multi-GiB source alone; (E) allocator dirt (jemalloc worsened Darwin Peak).

**Ownership reign-in (in tree):** detach-on-fill (`mem::replace` out of `all_records` into `region.reads` + owned `previous_region_reads`); finalize via `into_unique_records` first; drop RT/SeqGraph after haplotype extract; `GATK_RS_HC_RSS_TRACE=1` per-region / phase RSS log. Fail-closed 8 M DP / 10 k region-read ceilings remain safety nets.

**Remeasure 2026-08-06 (post detach, no jemalloc, sequential):** bomb **30.9 MiB** / 50 kb **48.3 MiB** exit 0. 100 kb live sawtooth **~55 → 633 → 2537 → 866 → 539 MiB** over ~3 min (within-region spike still GiB-class; inter-region baseline recovers — confirms sawtooth theory). Earlier `time -l` Peak **~1817 MiB** then jetsam. Full exit-0 still owed.

### Within-region Peak cut (2026-08-06 follow-up)

| Change | Where |
|--------|--------|
| Mid-phase RSS sampler + shared locus (`GATK_RS_HC_RSS_TRACE=1`) | `runtime_config.rs` — samples every 100 ms when RSS ≥800 MiB |
| Clear previous-region Arc pin during `callRegion`; commit after | `for_each_assembly_region` + `AssemblyRegionIterator` |
| Unique finalize take + re-share untrimmed (`call_region_mut`) | `assembly_based_caller.rs`, `engine.rs`, `run.rs` |
| Owned k-best (no `graph.clone` on cycle strip) | `kbest_haplotype.rs` |

**Named spike (TRACE):** first climb past 800 MiB was **`20:10098169-10098441`** (`reads=97`) — after finalize ~22 MiB, live RSS rose **~150→870 MiB/s** before `after_assemble`.

**Root cause (2026-08-09):** `supplement_p12_cluster_coupled_haplotypes` walked **all** configured+expanded k-mers whenever `scoring` was set. Under `strict_java`, that is every region → Peak on bushy non-P12 loci. Fix: non-P12 **early-stops** after alts (then 2 empty extracts) or 4 consecutive empty RT extracts; P12 keeps the full coupled-bridge walk. SeqGraph path still merges min variation k-mer separately.

**Remeasure (workspace `target/release`, unset `CARGO_TARGET_DIR`, no jemalloc, sequential, 800 MiB fuse):**

| Window | Peak-RSS | nvars | vs prior |
|--------|----------|-------|----------|
| bomb | **28.5 MiB** | 12/12 | was ~28 MiB / 12 |
| spike `20:10098000-10098600` | **23.5 MiB** | 7 | was **~884 MiB** SIGKILL |
| 50 kb | **47.2 MiB** | 23/31 | 8 sites still need later RT k-mers (Peak tradeoff; empty-streak-only reopens spike Peak) |
| 100 kb | **64.1 MiB** exit 0 | 47 | was **~2 GiB** / jetsam |
| **500 kb** `20:10000000-10500000` | **177.0 MiB** exit 0 | 104 | first clean exit-0 on this host post-Peak-cut (screen+caffeinate+fuse) |

**Unset `CARGO_TARGET_DIR` when measuring** — Cursor sandbox cache can leave a stale `target/release/gatk-rs`. 100 kb + 500 kb exit 0 with recorded Peak are in hand; holdout HC sanity re-checked below. Still **do not** sign `ci-subset` until FORMAT/F1 holdout gates + doctrine checklist are green. Optional: recover the 8 missing 50 kb sites without reopening the spike Peak.

Overnight recipe unchanged: **no jemalloc**, `GATK_RS_HC_SEQUENTIAL=1`, `RAYON_NUM_THREADS=1`.

**Related (not the HC BAM spike):** CombineGVCFs / GenotypeGVCFs still `read_all_records` whole inputs (cohort/chr-scale multi-GiB risk). ReferenceWindowCache warns if a long contig is loaded whole without `.fai`.

**Holdout Rust HC sanity (2026-08-09 post-Peak-cut, same env, 50 kb windows, 800 MiB fuse):**

| Window | Interval | Peak-RSS | Variants | Exit | Prior anchor |
|--------|----------|----------|----------|------|--------------|
| chr21 | `21:41200001-41250000` | **38.55 MiB** | 24 | 0 | 28.27 MiB / 26 |
| chr20 holdout | `20:15000000-15050000` | **37.08 MiB** | 3 | 0 | 32.53 MiB / 4 |

Logs: `/tmp/hc-rss-holdouts/`. nvars −2/−1 vs prior anchors tracks the RT-supplement early-stop Peak tradeoff (same class as the 8 missing 50 kb sites). Peak remains ≪ fuse; no jetsam.

**FORMAT + F1 holdout gates (2026-08-09, Peak-cut `strict_java` rust VCFs regenerated; Java cached):** **FAIL — do not sign `ci-subset`.**

| Slice | Rust sites (Java) | F1 rust (need) | FORMAT |
|-------|-------------------|----------------|--------|
| holdout `20:15Mb` | 3 (55) | **0.074** (≥0.90×Java) | fail — soft AD/DP 2/2 |
| chr21 `21:41.2Mb` | 24 (230) | **0.197** (≥0.95×Java) | fail — soft AD/DP 11/22 |
| chr20 dense `20:10Mb` | 23 (144) | **0.237** (≥0.95×Java) | fail — hard GT 3/16 + soft |

Precision is fine (P≈1.0); **recall collapsed** vs Jul-22 rust (~221 chr21 sites, F1≈0.99). Root cause (2026-08-09): Peak-cut `MAX_KBEST_PATH_EDGES=256` was shorter than padded HC paths (~1–2 kb of k-mers) → RT/Seq k-best `paths=0` everywhere, dangling fragments pruned, sparse read-variation emits only. **Fix:** raise path-edge cap to **4096** (heap/expansion Peak guards retained; non-P12 RT-supplement early-stop retained). Remeasure (`/tmp/hc-callrate-verify/`, 800 MiB fuse): spike **36.3 MiB** / 14 sites; chr21 **92.9 MiB** / **221** sites; 50 kb **103 MiB** / 139 sites.

**FORMAT + F1 after path-edge fix (2026-08-09T051310Z, `run_l9_signoff_gates.sh`):** **PASS** (still do not auto-sign `ci-subset` without doctrine checklist).

| Slice | Rust sites (Java) | F1 rust | FORMAT |
|-------|-------------------|---------|--------|
| holdout `20:15Mb` | 55 (55) | **1.000** | pass |
| chr21 `21:41.2Mb` | 221 (230) | **0.993** | pass |
| chr20 dense `20:10Mb` | 139 (144) | **1.000** | pass |

Log: `parity/reports/l9_holdout_format_f1_20260809T051310Z.log`. Prior fail log: `parity/reports/l9_holdout_format_f1_20260809T041657Z.log`.

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

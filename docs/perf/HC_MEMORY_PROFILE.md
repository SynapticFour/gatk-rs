# HaplotypeCaller memory profile (reproducible)

**Generated (UTC):** `20260724T181512Z` (failure) + engineering follow-up `20260724T213000Z`  
**Host:** `Darwin 25.5.0 arm64` (16 GiB MacBook Air — do **not** re-run the 2 Mb window here)  
**Runner script:** [`scripts/perf/run_hc_memory_profile.sh`](../../scripts/perf/run_hc_memory_profile.sh)  
**Raw run directory (original failure):** `docs/perf/runs/20260724T181512Z/`

**Public memory claim status:** **not allowed** until realistic profile is re-run on dedicated `gatk-rs-benchmark` host after the k-best fix below.

Profiles measured: `smoke, realistic` (realistic originally **failed**).

## Root cause (engineering, 2026-07-24)

The ~60 GiB “peak memory footprint” on `20:10000000-12000000` was **not** a measurement artifact.

| Finding | Evidence |
|--------|----------|
| Abrupt / locus-triggered | 10–50 kb OK (~90–114 MiB); ≥100 kb exploded |
| Minimal bomb window | `20:10098500-10099500` and especially active region `20:10098169-10098441` |
| Mechanism | Unbounded **k-best haplotype frontier** on bushy/cyclic assembly graphs: `BinaryHeap<PathState>` + per-push `path_bases()` / edge-list growth; cyclic graphs previously skipped cycle stripping |
| Secondary | Dead `span_records` BAM clone in the assembly iterator; SW/PairHMM lacked contig-scale refusal (chr20 × ~read → ~60 GiB matrices) |

### Fixes landed

1. **k-best** (`kbest_haplotype.rs` / `seq_kbest_haplotype.rs`): heap cap, path-edge cap, expansion cap; prefer cycle stripping before preserving topology.
2. **SW / PairHMM**: refuse oversized DP before allocate/scan.
3. **Iterator**: remove unused `span_records` full-span BAM clone.
4. **PR gate**: [`scripts/ci/check_hc_rss_regression.sh`](../../scripts/ci/check_hc_rss_regression.sh) — always unit bounds + optional HC on the 1 kb bomb window (≤256 MiB).

## Deep-pile / P12 hosted OOM (engineering, 2026-07-31)

GIAB smoke window `2:92300000-92350000` (P12 spine on full HG001 30×) killed hosted Rust HC (`runner has received a shutdown signal`) after RSS climbed ~2.3 GiB → ~15 GiB.

| Finding | Evidence |
|--------|----------|
| Positional DS **is** on | Production path uses `WalkerTraversalConfig::gatk_haplotype_caller_production` → cap 50/start |
| Residual after DS | ~2.1M raw → **~537k kept** (~20k unique starts); mid-window mean depth ~41k× |
| Not missing Java DS | GATK 4.4 has no extra per-region 1000-read cap; staggered starts defeat per-start caps |
| Amplification | Shard materializes all post-DS records; Rayon held **N** regions each cloning overlapping BAM reads |

### Fixes landed (parity-preserving)

1. **Sequential large regions** (`run.rs`): regions with ≥8192 reads flush alone (no Rayon siblings) — Java-like one-heavy-region peak shape; **same evidence**.
2. **Hard refuse** (`MAX_READS_PER_ASSEMBLY_REGION` = 100 000): fail closed like PairHMM oversized DP when a single region still exceeds a safe ceiling — **not** a genotype-contract downsampler.
3. **RSS unit**: `oversized_assembly_region_read_count_is_refused` in `check_hc_rss_regression.sh`.

### Local proof (2026-08-01, 16 GiB Mac)

| Setup | Result |
|-------|--------|
| Full HG001 30× P12, `--threads 4` | Peak-RSS ~2.7 GiB then abnormal exit (macOS peak footprint ≫ RSS) |
| Full HG001 30× P12, `--threads 1` | **exit 137** (SIGKILL) after ~6 min; max post-DS overlap ~32k reads/500 bp (under 100k refuse) |
| Post-DS residual | ~537k reads in 50 kb; refuse ceiling does not fire |

**Smoke hygiene:** `GIAB_MODE=smoke` stages P12 from **NA12878_20k** (`giab_stage_smoke_bam_hybrid`); chr20/21 remain HG001 30×. Full-30× P12 stays a benchmark-host / dedicated-RAM gate — no silent force-DS.

### Post-fix spot checks (this host, RSS watchdog)

| Window | Mode | Peak-RSS | Result |
|--------|------|----------|--------|
| `20:10098500-10099500` (1 kb bomb) | full HC | **~38 MiB** | OK, 12 variants |
| `20:10000000-10050000` (50 kb) | full HC | ~114 MiB (pre-kbest ladder) | OK |
| `20:10000000-12000000` (2 Mb) | full HC | — | **Not re-run on 16 GiB laptop**; use benchmark host |

## A. Trivial smoke — reproducibility reference only

**Label:** Trivial smoke (reproducibility only)  
**Interval:** `chr1:1-32`  
**BAM / ref:** checked-in `parity/fixtures/`  

> **Not for marketing.** Peak-RSS here is dominated by JVM/runtime fixed cost
> on a 32 bp window. Do **not** derive “X% less memory” from this table.

| Engine | Peak RSS | Wall time |
|--------|----------|-----------|
| **gatk-rs** (Rust release) | **9.52 MiB (9744 KiB)** | n/a |
| **Java GATK 4.4.0.0** | **437.49 MiB (447988 KiB)** | 4.4 s |

| Java / Rust Peak-RSS | 45.98× |
| Rust as fraction of Java Peak-RSS | 2.2% |
| Absolute delta (Java − Rust) | 427.97 MiB |

## B. Realistic GIAB-dense window — public-claim basis

**Label:** Realistic GIAB-dense multi-Mb window  
**Interval:** `20:10000000-12000000` (multi-Mb; default 2 Mb on chr20 dense locus)  
**BAM:** staged NA12878 NIST 30× slice (`parity/realworld/na12878_giab_window_mem_2mb_b37/`)  

> This realistic profile was **not** measured on the dedicated `gatk-rs-benchmark` host (see [`HOST_SPECS.md`](HOST_SPECS.md) / [`docs/ci/PERF_BENCHMARK_HOST.md`](../ci/PERF_BENCHMARK_HOST.md)). Numbers below are engineering evidence only — **do not** use them for a public “X% less memory” claim until re-run on that host.

**Status (original `20260724T181512Z`):** measurement **failed** — Rust Peak-RSS ~2.9 GiB, peak memory footprint ~60 GiB, process died (`20001 traversal tile(s)`).

**Status after k-best fix:** 1 kb bomb window recovers (~38 MiB). **Re-run the full 2 Mb realistic profile on `gatk-rs-benchmark`** (or any host with ≥32 GiB RAM) before updating public claims:

```bash
export CARGO_TARGET_DIR="$PWD/target"   # avoid stale sandbox binaries
HC_MEM_PROFILES=realistic ./scripts/perf/run_hc_memory_profile.sh
```

## Exact commands

### Rust (build once)

```bash
cargo build -p gatk-cli --release --locked
```

### Memory regression (PR Check)

```bash
./scripts/ci/check_hc_rss_regression.sh
# optional override when assets staged:
#   HC_RSS_INTERVAL=20:10098500-10099500 HC_RSS_MAX_MIB=256 ./scripts/ci/check_hc_rss_regression.sh
```

### Java GATK 4.4

- Pin: `GATK_PINNED_SHA=2dbc025821bc5f686c423ff332a41e6cef892a77` (`docs/GATK_PINNED.env`)
- Image: `us.gcr.io/broadinstitute/gatk:4.4.0.0` (or Broad `us.gcr.io/broad-gatk/gatk:4.4.0.0`)
- JVM options (pipeline-realistic): `-Xms1g -Xmx4g`

```bash
./scripts/perf/run_hc_memory_profile.sh
# smoke only:
#   HC_MEM_PROFILES=smoke ./scripts/perf/run_hc_memory_profile.sh
# realistic only (stages 2 Mb GIAB window if needed):
#   HC_MEM_PROFILES=realistic ./scripts/perf/run_hc_memory_profile.sh
```

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

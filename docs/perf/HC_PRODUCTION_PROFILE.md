# Production HaplotypeCaller profiling

Observe-only instrumentation for **real** HC workloads (not Criterion microbenches).
Does **not** change genotype / emit results.

## Enable

```bash
export GATK_RS_HC_PROFILE=/path/to/out/hc_profile.json   # or a directory
# Recommended: also enable TRACE so assemble sub-phases dual-write into the profile
export GATK_RS_HC_RSS_TRACE=1
export RAYON_NUM_THREADS=2
# unset GATK_RS_HC_SEQUENTIAL   # product wall
```

Or use the wrapper:

```bash
./scripts/perf/run_hc_profile.sh \
  --out-dir docs/perf/runs/hc_profile_w09 \
  -- \
  -R parity/realworld/assets/hs37d5.simple.fa \
  -I parity/realworld/na12878_ci_loser_windows/01_chr21_w09.bam \
  -O /tmp/w09_profile.vcf \
  -L 21:9500000-9700000
```

Outputs:

| File | Contents |
|------|----------|
| `hc_profile.json` | Machine-readable `gatk-rs.hc_profile.v1` |
| `hc_profile.md` | Human-readable stage / PairHMM / genotype tables |
| `hc_rss.trace` | Existing `HC_RSS_TRACE` lines (when TRACE on) |

## What is measured

### Stages (wall + process CPU when available)

`input_bam_decode`, `read_preprocessing`, `active_region_construction`,
`event_discovery`, `assembly_graph_construction`, `graph_pruning`,
`haplotype_generation`, `smith_waterman`, `pairhmm`, `likelihood_processing`,
`genotype_assignment`, `ad_annotation`, `vcf_emission`,
`synchronization_waiting`, `allocations`, `other`

Per stage: **calls**, **wall**, **CPU** (when sampled), **avg wall/call**,
**alloc bytes/events** (best-effort).

Sources:

- Explicit RAII guards: BAM fill, active-region next, finalizeRegion, PairHMM,
  likelihood normalize/filter, realign SW, genotype assign, AD pileups, VCF emit,
  post-wave TLS barrier
- Dual-write from `rss_trace_checkpoint` phase names → coarse stages (assemble
  graph / prune / k-best / etc.)
- Run-level CPU via `getrusage(RUSAGE_SELF)`; parallel efficiency ≈ `cpu / (wall × threads)`

### PairHMM extras

- reads / haplotypes / read×hap pairs, mean haps/read
- read & hap length histograms (25 bp buckets)
- SIMD pack units, pack occupancy %, prefix-reuse %, leftover %
- DP cells evaluated / avoided via `hapStartIndex` (when prefix path runs)

### Genotyping extras

- sites, candidate alleles, genotype states, PL vector sizes
- time/site and time/state (region wall amortized per call)
- AD wall (pileup AD calls) and event-rebuild wall (per-hap EventMap build)

## Honesty

- Stage walls from TRACE dual-write can **overlap** nested guards; use PairHMM /
  genotype nested counters for leaf detail.
- `synchronization_waiting` currently covers post-wave TLS barriers, not full
  rayon idle time — use `parallel_efficiency` for imbalance.
- Allocation bytes are best-effort (`note_alloc_bytes`); not a full allocator trace.
- Prefer product thr=2 + `FASTEST_AVAILABLE` for wall decisions; Peak
  (`GATK_RS_HC_SEQUENTIAL=1`) is a different shape.

## Code

- `gatk-haplotypecaller/src/hc_profile/`
- Env inventory: `GATK_RS_HC_PROFILE` in `runtime_config::DebugConfig`
- Wrapper: `scripts/perf/run_hc_profile.sh`

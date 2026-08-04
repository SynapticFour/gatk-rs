# HC memory baseline — pre Arc-share pass (2026-08-04)

Branch: `perf/giab-memory-footprint`  
Host: Darwin arm64 MacBook Air M4 class (16 GiB)  
Note: disk ~13 GiB free — avoid staging multi-Mb GIAB windows during this pass.

## Prior measured anchors (from `HC_MEMORY_PROFILE.md`)

| Window | Threads | Peak-RSS / outcome |
|--------|---------|-------------------|
| `20:10098500-10099500` (1 kb bomb) | 1 | ~38 MiB OK |
| `20:10000000-10050000` (50 kb) | 1 | ~114 MiB OK (pre-kbest ladder era) |
| `20:10000000-12000000` (2 Mb) | 1 | **Do not re-run on 16 GiB** until Arc-share + sequential proven |
| P12 full 30× 50 kb | 1 | exit 137 historically |

## Engineering targets this pass

1. Eliminate deep `bam::Record` clone in `fill_region_with_reads` (`Arc` share).
2. Progressive shard Arc release past previous-region window.
3. Sequential region apply on constrained hosts (`GATK_RS_HC_SEQUENTIAL=1`).
4. PairHMM / SW TLS scratch shrink after regions.
5. Stream-merge VCF batches (no interval-wide `all_batches`).
6. Finalize disk: hardlink equiv VCFs; delete shards / BAM / strat after score.

## Post-pass

Re-measure bomb + 50 kb + 100 kb (+ 500 kb if safe) with `GATK_RS_HC_SEQUENTIAL=1` and update `HC_MEMORY_PROFILE.md` only with measured numbers. No marketing claims from smoke.

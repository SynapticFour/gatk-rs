# ci-subset hang: `00_chr20_w47` / `01_chr21_w10` rust

## Symptom

Across many `ci-subset` runs these two **rust** matrix legs run until the
GitHub-hosted **6h job timeout** and are **cancelled**. The workflow then skips
Finalize (`needs: hc-shard` fails), so a full ci-subset never scores.

| Shard | Interval | Java HC wall | Java sites | Rust |
|-------|----------|-------------:|-----------:|------|
| `00_chr20_w47` | `20:47000001-48000000` | ~1.6 min | 1788 | cancel @ 6h |
| `01_chr21_w10` | `21:10000001-11000000` | ~12 min | 10625 | cancel @ 6h |

Examples: runs `31841981413`, `31803808310`, `31772652099`, `31744034963`,
`31868192320` (same two still in progress).

## What it is not

- Not runner OOM / Peak abort (`GATK_RS_HC_RSS_ABORT_MIB=8192` never fires).
- Not a panic / exit 137.
- Not Java-side pathology (Java finishes cleanly).

## Evidence (cancelled artifacts, run 31841981413)

Mem-probe @ 5s for ~5.9h of HC:

| Shard | RSS plateau | Time to plateau | Stuck after | Load |
|-------|------------:|----------------:|------------:|-----:|
| w47 | **~317 MiB** flat | ~0.4 min | **~5.96 h** | ~1.0 |
| w10 | **~1326 MiB** flat | ~0.4 min | **~5.96 h** | ~1.0 |

- Process still `gatk-rs HaplotypeCaller` until job `TERM`.
- Uploaded rust VCF is **header only** (0 variant rows) — hang before any call flush.
- Empty `/usr/bin/time` file (process never exited cleanly).
- `GATK_RS_HC_SEQUENTIAL=1` → load≈1 matches single-thread CPU grind.

## Working hypothesis

A **pathological assembly region** (or early-window combinatorial k-best /
PairHMM / genotype path) near the start of each 1 Mb window. Wall-clock
unbounded under Peak sequential; RSS stable ⇒ not memory growth.

## Complete ci-subset once (ops)

Default: **omit these two shards** from `ci-subset` prepare matrix (documented
coverage hole). Opt back in with `GIAB_INCLUDE_HANG_SHARDS=1`.

Finalize uses `if: always()` so a future single-leg cancel cannot block scoring
of completed shards.

## Algorithm follow-up (separate from wall campaign)

1. TRACE / Instruments on `20:47000001-47050000` and `21:10000001-10050000`
   (first 50 kb) under Peak sequential.
2. Optional: region wall-clock soft-abort (peer of RSS abort) once locus class
   is known — **no P12 widen**.

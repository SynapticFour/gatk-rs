# ci-subset hang: `00_chr20_w47` / `01_chr21_w10` rust — **fixed**

## Symptom (historical)

Across many `ci-subset` runs these two **rust** matrix legs ran until the
GitHub-hosted **6h job timeout** and were **cancelled**. Finalize then skipped
(`needs: hc-shard` fails), so a full ci-subset never scored.

| Shard | Interval | Java HC wall | Java sites | Rust (pre-fix) |
|-------|----------|-------------:|-----------:|-----------------:|
| `00_chr20_w47` | `20:47000001-48000000` | ~1.6 min | 1788 | cancel @ 6h |
| `01_chr21_w10` | `21:10000001-11000000` | ~12 min | 10625 | cancel @ 6h |

Examples: runs `31841981413`, `31803808310`, `31772652099`, `31744034963`,
`31868192320`.

## What it was not

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
- Uploaded rust VCF was **header only** (0 variant rows).
- `GATK_RS_HC_SEQUENTIAL=1` → load≈1 matches single-thread CPU grind.

## Root cause

Local quarter-slice bisect (`parity/giab/runs/hang-repro/`): only **`w47_q1`**
hung; hot locus **`20:47131737-47131971`**.

TRACE: `rt_graph_built` (~4461 nodes) then stall **before** `kbest_begin`.

In `graph_for_kbest` → `find_cycle_guilty` (`kbest_haplotype.rs`): recursive DFS
over a **bushy diamond DAG** that also fails `has_cycle()`, **without per-vertex
memoization of `reaches_sink`**. Path count is exponential in diamond depth →
CPU hang at flat RSS until the 6h cancel.

## Fix

Memoize `reaches_sink: Vec<Option<bool>>` in `find_cycle_guilty` (same observable
cycle-edge / non-reaching-vertex removal contract; linear in nodes+edges on this
class). TRACE checkpoints: `kbest_cycle_strip_begin` / `kbest_cycle_strip_done`.

Unit: `graph_for_kbest_memoizes_bushy_cyclic_dag`.

Local Peak-sequential proof (post-fix):

| Window | Wall | Sites | Notes |
|--------|-----:|------:|-------|
| `w47_hot` (`20:47130000-47140000`) | ~10 s | 56 | hot locus cycle-strip **~1.5 ms** → `kbest_begin` |
| `w47_q1` (250 kb) | ~20 s | 412 | cycle-strip **~1.9 ms** |
| full `w47` 1 Mb | **~69 s** | 1208 | no hang |
| full `w10` 1 Mb | ~20–40 min (Peak seq.) | (flush at end) | progresses; was 6h flat hang pre-fix |

## Ops note

Temporary `GIAB_INCLUDE_HANG_SHARDS` omit of w47/w10 was **removed** once the
algorithm fix landed. Finalize retains `if: always()` so a future single-leg
cancel cannot block scoring of completed shards.

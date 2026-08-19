# Read-threading k-mer representation (packed keys)

## Problem

Production HC defaults try `k ∈ {10, 25}` then expands by +10. The previous builder
allocated `Arc<[u8]>` for **every sliding window** in preprocess / find-start /
extend, then keyed maps with those Arcs. On deep regions that is O(reads ×
usable_len) heap traffic plus RC bumps.

## Solution (in tree)

[`kmer_key.rs`](../../gatk-haplotypecaller/src/kmer_key.rs) + wiring in
[`read_threading_graph.rs`](../../gatk-haplotypecaller/src/read_threading_graph.rs):

| Window | Key |
|--------|-----|
| Pure ACGT, `k ≤ 32` | `u64` (2 bits/base, MSB-first) |
| Pure ACGT, `33 ≤ k ≤ 64` | `u128` |
| Contains `N` / non-ACGT, or `k > 64` | `Arc<[u8]>` (exact bytes) |

Preprocess / find-start use [`RollingKmer`](../../gatk-haplotypecaller/src/kmer_key.rs)
(O(1) advance after the first window when bases stay ACGT).

Node payloads remain `Arc<[u8]>` (path bases / dumps). Unique / non-unique maps
hash `KmerKey` — **zero alloc** for the common ACGT k=10/25 preprocess path.

`AssemblyGraph::kmer_to_id` is a `HashMap` (lookup only). Dump/sort still orders
by k-mer bytes where needed.

## Comparison of approaches

| Approach | Allocations (preprocess) | Hash / compare | Status |
|----------|--------------------------|----------------|--------|
| **A. `Arc<[u8]>` per window** | 1 Arc + RC per window | hash/eq over `k` bytes | previous default; kept as fallback for N / k>64 |
| **B. Integer-packed `KmerKey`** | 0 for ACGT k≤64 | word hash/eq | **shipped** (+ rolling) |
| **C. Interned vertex IDs** | dictionary + id maps | u32 compare | deferred — B already removes dominant Arc churn for production k |

## Deterministic order (must keep)

`outgoing` / `incoming` stay **`BTreeSet`**. Iteration order drives
`extend_chain_by_one` first-suffix-match when multiple outs share a base —
that is an **observable** topology choice. Documented in the builder module docs.

Where order is **not** observable:

- `unique_kmers` / `non_unique_kmers` → `HashMap` / `HashSet`
- `AssemblyGraph::kmer_to_id` → `HashMap`
- edge storage → `HashMap` (sorted only when dumping)

`read_threading_assembler` still uses `BTreeSet` for `VariationEvent` collections
where sorted event identity is part of deterministic haplotype discovery.

## Prove

```bash
cargo test -p gatk-haplotypecaller --lib kmer_key --locked
cargo test -p gatk-haplotypecaller --lib read_threading_graph --locked
cargo test -p gatk-haplotypecaller --test indel_threading_graph --locked

cargo bench -p gatk-haplotypecaller --bench assembly_graph --locked -- \
  'kmer_key_representation|threading_build'
```

## Measured (local Criterion, darwin, 2026-08-19)

Fixtures: ~120 bp ACGT reads with small indel alleles; depths 16 / 64 / 512.

### A vs B — sliding-window map fill (`kmer_key_representation`)

| k | A `Arc<[u8]>` | B packed+rolling | Speedup |
|---|---------------|------------------|---------|
| 10 | 244 µs | 106 µs | **~2.3×** |
| 25 | 278 µs | 110 µs | **~2.5×** |

(Median of Criterion `time` interval midpoints; re-run on your host for CI noise.)

HashMap insert/lookup still dominates residual time. Memory: A stores `k` bytes +
Arc header per distinct key during fill; B stores 8 or 16 bytes per packed key.

### Graph construction (`assembly_graph_depth/threading_build_*`)

Packed builder absolute times (post-change, realistic fixtures):

| Depth | k=10 | k=25 |
|-------|------|------|
| low (16) | 193 µs | 219 µs |
| medium (64) | 593 µs | 721 µs |
| high (512) | 4.31 ms | 4.86 ms |

Earlier Criterion “regressions” vs prior baselines were from **short fixtures**
that could not form k=25 windows (near-empty early exit). Do not compare those.

Cache misses / heap profiles: not collected here; use Instruments /
`dhat` if needed. Dominant win is eliminating per-window `Arc` in preprocess.

## Non-goals

- No change to pruning multiplicity / candidate haplotype algorithm
- Junction-tree / nearby-kmer corrector still use `Vec<u8>` maps (secondary paths)
- Interned integer vertex IDs for *all* k-mers (option C) deferred

# Smith-Waterman for production HaplotypeCaller

## Role in wall time

Phase 7 Instruments (`docs/perf/PHASE7_WALL_PROFILE.md`) ranked
`smith_waterman::align_uppercase_ready` **#3** on a mid-run GIAB sample (behind
PairHMM scoring and RT assembly). SW matters, but typical matrices are **short**
(read×hap ≈ 100–250; padded hap-to-ref ≈ 70–320), so a Farrar/striped SIMD
rewrite is deferred until a rematch shows SW still dominates after this leaf.

## Investigation summary

| # | Question | Finding | Action |
|---|----------|---------|--------|
| 1 | Full score + backtrack both needed? | End-cell pick needs only last row + last column of scores; CIGAR follows `btrack` only | **Rolling score rows** + `last_col[]` |
| 2 | Rolling scores when only final score? | HC always needs CIGAR, never score-only | N/A for main path; rolling still applies |
| 3 | Compact backtrack? | Track stores ±gap length up to dim; under `MAX_SW_DIM=100k` needs `i32`. `i16` unsafe without lowering the cap | Keep `i32` btrack |
| 4 | Narrower score ints? | Hap-to-ref `match=200` overflows `i16` by ~len 164 | Keep `i32` scores |
| 5 | SIMD scoring kernel? | Affine gap has left-to-right dependency (`best_gap_h`); striped SW is non-trivial | Deferred |
| 6 | Striped/blocked SW? | Setup tax vs short HC lengths; parity risk | Deferred |
| 7 | Exact / prefix / small-indel paths? | Exact substring via `lastIndexOf` already; equal-len SNP/MNP skipped in `haplotype_cigar` | Keep; no new heuristic indel path |
| 8 | Uppercase copies? | Skipped when already ASCII-upper (HC common) | Keep |
| 9 | Scratch / TLS? | Dropped `mem::take` restore; work in place under one `RefCell` borrow; `clear()` now drops Peak | Done |
| 10 | CIGAR backtrack bottleneck? | O(path) ≪ O(nrows·ncols) for HC sizes | No change |

## Layout (current)

- `btrack`: flat `nrow × ncol` `i32`
- scores: `sw_prev` / `sw_cur` (`ncol` each) + `last_col` (`nrow`)
- SoftClip / Ignore end cell from `last_col` + last row; Indel corner unchanged

Memory for scores: \(O(n+m)\) instead of \(O(nm)\); backtrack unchanged.

## Prove

```bash
cargo test -p gatk-haplotypecaller --lib smith_waterman --locked
cargo test -p gatk-haplotypecaller --lib haplotype_cigar --locked

cargo bench -p gatk-haplotypecaller --bench smith_waterman --locked
```

Oracle: `rolling_matches_full_matrix_oracle` compares production CIGAR/offset to an
independent full-matrix fill.

## Measured (local Criterion, 2026-08-19)

SoftClip / Indel vs prior full-score-matrix baseline (Criterion `change`):

| Size | SoftClip before → after | Δ |
|------|-------------------------|---|
| 64×48 | 9.0 → 7.3 µs | **≈−19%** |
| 128×96 | 38 → 30 µs | **≈−20%** |
| 256×192 | 150 → 122 µs | **≈−19%** |

HC-focused (`smith_waterman_hc`), absolute post-change:

| Workload | Time |
|----------|------|
| padded Indel core 80 / 150 / 250 | ~26 / 73 / 183 µs |
| read→hap SoftClip 120×100 / 200×151 / 280×151 | ~30 / 74 / 104 µs |
| exact substring fast path | ~138 ns |

Exact-substring SoftClip stays on the `lastIndexOf` path (sub-µs class; see
`exact_substring_fast_path` bench).

Score-plane TLS footprint drops from \(O(nm)\) to \(O(n+m)\); backtrack plane
unchanged (still the Peak driver for large cells).

## Gate to start striped SIMD (post #128 rematch)

Local 200 kb heads on remaining wall-losers still rank **PairHMM #1, SW #2**
(w11 / w26 / w29). SW share remains large enough to justify a striped experiment,
but only with:

1. Existing oracle `rolling_matches_full_matrix_oracle` extended to striped path
2. SoftClip / Indel Criterion benches showing wall↓ on HC lengths (80–280)
3. No CIGAR/pos drift vs production rolling path

Until those land, keep rolling + `last_index_of` — do **not** merge striped without
oracle proof.

## Non-goals

- No striped SIMD SW without the gate above
- No score-only API (HC always needs CIGAR)
- No P12 band widening / heuristic indel shortcuts that alter CIGAR

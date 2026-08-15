# Phase 9 — PairHMM pack occupancy (product wall)

Goal: cut `after_pairhmm` share (~42% on w11 product TRACE) toward &lt;1.0× Java
median on `wall-losers` (baseline ~1.79×).

## Landed in this phase

| Change | Contract |
|--------|----------|
| NEON equal-length **`by_len`** packs (AVX2-style HashMap group) | Same Logless numerics; replaces look-ahead walk |
| TRACE `neon_pack2=` / `neon_leftover=` on `after_pairhmm` | Observe-only; no genotype effect |

## Deferred

- Equal-length **read** packs (transpose SIMD axis) — highest remaining leverage on Darwin.
- GKL-style hap flank windows — only if ≡ scalar Logless under SIMD unit test.
- Shared-prefix prior/column reuse across haps.
- `f32` SIMD promotion (Criterion slower than f64).

## Prove

```bash
cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test --locked
# Product TRACE rematch (thr=2, no sequential):
# after_pairhmm Σδ vs callrate-era; inspect neon_pack2 / neon_leftover rates
```

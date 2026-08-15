# Phase 9 — PairHMM pack occupancy + hap prefix reuse (product wall)

Goal: cut `after_pairhmm` share (~42% on w11 product TRACE) toward &lt;1.0× Java
median on `wall-losers` (baseline ~1.79×).

## Landed

| Change | Contract |
|--------|----------|
| NEON equal-length **`by_len`** packs (AVX2-style) | Same Logless numerics (#118) |
| TRACE `neon_pack2=` / `neon_leftover=` on `after_pairhmm` | Observe-only |
| Occupancy dig: w11 50 kb leftover **~6.7%** | Hap packing already dense |
| Java **`hapStartIndex` prefix reuse** in packed Logless | Same scores as full recompute; fewer DP cells on shared assembly prefixes |
| NEON: same-length groups of **≥3** use prefix-reuse chain; pairs stay pack2 | SIMD kept where packs are short |

## Deferred

- Equal-length **read** packs (transpose SIMD axis).
- GKL flank windows — only if ≡ scalar Logless under SIMD unit test.
- `f32` SIMD promotion (Criterion slower than f64).

## Prove

```bash
cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test --locked
# Product TRACE rematch vs callrate-era / phase9 pack baseline
```

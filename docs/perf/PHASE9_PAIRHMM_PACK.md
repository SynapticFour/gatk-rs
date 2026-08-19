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
| TRACE `neon_prefix_reuse=` separate from `neon_leftover=` | Observe-only (prefix ≠ scalar leftover) |
| Order same-length haps by sequence before reuse/pack | Longer shared prefixes; score-invariant |
| AVX2: same-length groups of **≥5** use prefix-reuse (pack4 for shorter) | Dense GIAB shared prefixes |
| NEON: TLS `by_len` + in-place sort (no idxs clone); leftover via `score_one_hap` | Same scores; less alloc churn |
| AVX2: TLS `by_len` + leftover via `score_one_hap` (mirror NEON) | Same scores; less alloc churn |

## Prefix-vs-pack A/B (2026-08-19, post wall-losers 1.15× rematch)

| Knob | Result | Keep? |
|------|--------|------:|
| `PREFIX_REUSE_OVER_SIMD_FRAC` 0.35→0.50 | Occupancy unchanged (~3% packs / ~95% prefix) on w11 200 kb — mean consecutive prefix already &gt;0.50 | **REVERT** |
| `PREFIX_REUSE_MIN_HAPS_NEON` 3→6 / AVX2 5→8 | Occupancy **identical** (same pack/prefix counts) — groups already large | **REVERT** |

Named constants kept at original 3/5. Next PairHMM bet: read-axis packs / wavefront
with TRACE Σδ — hap-axis prefix path is saturated on dense losers.

## Deferred

- Equal-length **read** packs (transpose SIMD axis) — only after product TRACE shows
  hap-axis still leaves wall-losers median &gt;1.0×.
- GKL flank windows — only if ≡ scalar Logless under SIMD unit test.
- `f32` SIMD promotion (Criterion slower than f64).
- **Architecture C row-wavefront** — implemented opt-in as `--pair-hmm WAVEFRONT`
  ([`PAIRHMM_WAVEFRONT.md`](PAIRHMM_WAVEFRONT.md)); not in `FASTEST_AVAILABLE` yet.
  True GKL anti-diagonal remains the next measurement target after wavefront TRACE.

## Prove

```bash
cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test --locked
# Product TRACE rematch vs callrate-era / phase9 pack baseline
```

## TRACE note (phase10 tip, local product thr=2)

On w11 densest-ish 50 kb, `after_pairhmm` Σ ≈ 31–36 s vs assign ≈ 28–32 s —
hap-axis still large but not sole leaf. Prefer wall-losers CI for signed ratios
before investing in read-axis packs.
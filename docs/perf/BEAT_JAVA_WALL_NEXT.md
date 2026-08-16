# Beat Java wall — hard look (phase6→10)

Goal: product-shaped wall (thr=2, no sequential) **&lt;1.0×** Java on dense GIAB, with
algorithm parity (no P12 band widening). Peak RSS already wins on probe.

## What the evidence says now (phase10 / #120)

Signed product lane: `GIAB_MODE=wall-losers` on tip `7991f28`
([run 31940745824](https://github.com/SynapticFour/gatk-rs/actions/runs/31940745824)):

| Metric | Value |
|--------|------:|
| Median wall Rust/Java | **~1.62×** (was ~1.79× baseline) |
| Σ wall | **1.74×** (46m12 / 26m30) |
| Worst shard | `01_chr21_w09` **2.58×** |
| Best dense | `00_chr20_w09` **1.38×** |
| Peak RSS | ~0.14–0.61× Java (still wins) |
| Equivalence gate | **PASS** `max_\|ΔF1\|=0.0004` |

Peak-mode ci-subset dense median (~1.57×) is **not** the product-wall claim —
`GATK_RS_HC_SEQUENTIAL=1` there. Use wall-losers for beat-Java wall.

## Highest-leverage wall bets (evidence-class)

| Priority | Bet | Status |
|----------|-----|--------|
| 1 | **Collapse multi-pass AD + CIGAR/seq cache** | Landed (#120): TLS `AdDecodeCache`, pad/slice reuse, single-pass softclip counts |
| 2 | **Softclip-aware TLS base lookup** | Landed: `AdDecodeCache::softclip_base_at_ref_1based` (SoftClip-as-ref ≠ AD `query_index`) |
| 3 | **PairHMM leaf / packs** | Landed: NEON + AVX2 TLS `by_len` / leftover `score_one_hap`; **read-axis packs deferred** until TRACE proves hap-axis still leaves &gt;1.0× |
| 4 | **Gate parity_spine when redundant** | Indel CIGAR-complete skip kept; SNP skip only when no alt haps |
| 5 | **Realign SW** | Landed: `last_index_of` first-byte reject (phase8) — no further cheap reject |

Prove: `pairhmm_simd_vs_scalar_test`, softclip engine tests, mega TRACE + wall-losers rematch; no P12 band widening.

## Reject / defer

- Widening P12 bands or deleting pins to “win” wall.
- Nested Rayon on PairHMM (already collapsed to one axis).
- Full EventMap precache in allele filtering (regressed prep wall badly).
- Marketing genome-wide equivalence from Peak-mode CI Σ.
- Read-axis PairHMM packs before TRACE Σδ proof on current hap-axis.

## Immediate next experiment

1. ~~Phase6–9 / wall-losers baseline~~ — see table above (#112–#120).
2. **chr21 w09 genotype** — softclip TLS + profile `try_genotype` / AD / htslib (2.58× shard).
3. **AVX2 PairHMM hygiene** (TLS `by_len` + leftover `score_one_hap`); TRACE `after_pairhmm` before read-axis.
4. **Call-rate / L8 (separate track)** — #116 emit pileup AD landed; wall-losers still ~7k strong SNP hets with **no Rust POS** → discovery/EventMap undercall, not emit-depth. See [`CALLRATE_EMIT_AD_DIG.md`](CALLRATE_EMIT_AD_DIG.md).

### Local TRACE rematch (phase10 tip, product thr=2)

| Window | assign_genotype_s | pairhmm_s | realign_s |
|--------|------------------:|----------:|----------:|
| mega `21:9825–9828k` | ~147–172 | ~15–20 | ~4–7 |
| w10 hot 50 kb | ~14 | ~12–14 | ~7–8 |
| w11 densest-ish 50 kb | ~28–32 | ~31–36 | ~14–18 |

Laptop assign noise is high; prefer CI wall-losers for signed wall ratios.

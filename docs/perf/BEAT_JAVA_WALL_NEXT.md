# Beat Java wall — hard look (phase6→10)

Goal: product-shaped wall (thr=2, no sequential) **&lt;1.0×** Java on dense GIAB, with
algorithm parity (no P12 band widening). Peak RSS already wins on probe.

## What the evidence says now (phase10 / #120 rematch)

Signed product lane: `GIAB_MODE=wall-losers` on `main` @ `6a478ca`
([run 32003572266](https://github.com/SynapticFour/gatk-rs/actions/runs/32003572266)):

| Metric | Value |
|--------|------:|
| Median wall Rust/Java | **~1.61×** (was ~1.62× @ 7991f28 / ~1.79× baseline) |
| Σ wall | **1.65×** (44m14 / 26m51; was 1.74×) |
| Worst shard | `01_chr21_w09` **2.01×** (was **2.58×**) |
| Best dense | `00_chr20_w09` **1.39×** |
| Peak RSS | ~0.16–0.62× Java (still wins) |
| Equivalence gate | **PASS** `max_\|ΔF1\|=0.0004` |
| Call-rate rust/java sites | **0.857** (26544 / 30984) |

Workflow red only on **Publish** (protected-main GH006) — Finalize PASS.
Peak-mode ci-subset dense median (~1.55×) is **not** the product-wall claim —
`GATK_RS_HC_SEQUENTIAL=1` there. Use wall-losers for beat-Java wall.

## Highest-leverage wall bets (evidence-class)

| Priority | Bet | Status |
|----------|-----|--------|
| 1 | **Collapse multi-pass AD + CIGAR/seq cache** | Landed (#120): TLS `AdDecodeCache`, pad/slice reuse, single-pass softclip counts |
| 2 | **Softclip-aware TLS base lookup** | Landed: `AdDecodeCache::softclip_base_at_ref_1based` (SoftClip-as-ref ≠ AD `query_index`) |
| 3 | **PairHMM leaf / packs** | Landed: NEON + AVX2 TLS `by_len` / leftover `score_one_hap`; **read-axis packs deferred** until TRACE proves hap-axis still leaves &gt;1.0× |
| 4 | **Gate parity_spine when redundant** | Indel CIGAR-complete skip kept; SNP no-alt skip **reverted** (blocked materialize / p5 j2 fixture) |
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
2. ~~Softclip TLS / AVX2 hygiene~~ — landed on #120 tip.
3. ~~Product wall rematch~~ — [32003572266](https://github.com/SynapticFour/gatk-rs/actions/runs/32003572266):
   w09 **2.58→2.01×**, Σ **1.74→1.65×**; still 8/8 rust_slower.
4. **chr21 genotype leaf** — local TRACE (product thr=2):
   - mega `21:9825–9828k`: **assign 130 s** vs pairhmm 9 s / realign 3 s; hot region
     `21:9826233` alone ~37 s assign.
   - w11 densest-ish: assign **32 s** / pairhmm **24 s** / realign **9 s**.
   Next wall cut = structural genotype/AD on mega loci, not PairHMM packs.
5. **Call-rate / L8** — spine strong-het gates + sync retain (see
   [`CALLRATE_EMIT_AD_DIG.md`](CALLRATE_EMIT_AD_DIG.md)); rematch 5/12 miss sites emitted.
6. **Publish hygiene** — stop `git push` to protected `main` in GIAB Publish jobs.

### Local TRACE rematch (phase10 tip, product thr=2)

| Window | assign_genotype_s | pairhmm_s | realign_s |
|--------|------------------:|----------:|----------:|
| mega `21:9825–9828k` | ~147–172 | ~15–20 | ~4–7 |
| w10 hot 50 kb | ~14 | ~12–14 | ~7–8 |
| w11 densest-ish 50 kb | ~28–32 | ~31–36 | ~14–18 |

Laptop assign noise is high; prefer CI wall-losers for signed wall ratios.

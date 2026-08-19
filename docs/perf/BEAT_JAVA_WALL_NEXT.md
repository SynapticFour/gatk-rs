# Beat Java wall — hard look (phase6→10 + ledger tip)

Goal: product-shaped wall (thr=2, no sequential) **&lt;1.0×** Java on dense GIAB, with
algorithm parity (no P12 band widening). Peak RSS already wins on probe.

## What the evidence says now (ledger tip / #128 rematch)

Signed product lane: `GIAB_MODE=wall-losers` on `perf/hc-wall-ledger-g1-p1-p2`
([run 32244936304](https://github.com/SynapticFour/gatk-rs/actions/runs/32244936304)):

| Metric | Value |
|--------|------:|
| Median wall Rust/Java | **~1.15×** (was **~1.61×** @ 6a478ca / 32003572266) |
| Σ wall | **1.27×** (32m54 / 25m53; was 1.65×) |
| Worst shard | `01_chr21_w11` **1.65×** (was w09 **2.01×**) |
| Best dense | `01_chr21_w29` **0.82×** (Rust faster) |
| Peak RSS | ~0.14–0.66× Java (still wins) |
| Equivalence gate | **PASS** `max_\|ΔF1\|=0.0004` |

Living method + experiment table: [`PERFORMANCE_LEDGER.md`](PERFORMANCE_LEDGER.md).

Peak-mode ci-subset dense median is **not** the product-wall claim —
`GATK_RS_HC_SEQUENTIAL=1` there. Use wall-losers for beat-Java wall.

## Highest-leverage wall bets (evidence-class)

| Priority | Bet | Status |
|----------|-----|--------|
| 1 | **Collapse multi-pass AD + CIGAR/seq cache** | Landed (#120): TLS `AdDecodeCache`, pad/slice reuse, single-pass softclip counts |
| 2 | **Softclip-aware TLS base lookup** | Landed: `AdDecodeCache::softclip_base_at_ref_1based` |
| 3 | **Genotype reshape / AD memo (G1)** | Landed (#128): TLS borrow likelihood rows + indel AD memo — mega assign no longer dominant |
| 4 | **PairHMM leaf / packs + prefix threshold** | Landed packs + hapStartIndex; **retune** `PREFIX_REUSE_OVER_SIMD_FRAC` with TRACE A/B on remaining losers |
| 5 | **Gate parity_spine when redundant** | Indel CIGAR-complete skip kept; SNP no-alt skip **reverted** |
| 6 | **Realign SW** | Rolling scores + oracle landed; striped/SIMD **deferred** until loser profiles show SW ≥ PairHMM |

Prove: `pairhmm_simd_vs_scalar_test`, softclip engine tests, mega TRACE + wall-losers rematch; no P12 band widening.

## Reject / defer

- Widening P12 bands or deleting pins to “win” wall.
- Nested Rayon on PairHMM (already collapsed to one axis).
- Full EventMap precache in allele filtering (regressed prep wall badly).
- Marketing genome-wide equivalence from Peak-mode CI Σ.
- Read-axis PairHMM packs before TRACE Σδ proof on current hap-axis.
- Striped SIMD SW before rematch shows SW still dominates after rolling leaf.

## Immediate next experiment

1. ~~Phase6–10 / wall-losers baselines~~ — see table above (#112–#128).
2. ~~G1/P1/P2 tip rematch~~ — [32244936304](https://github.com/SynapticFour/gatk-rs/actions/runs/32244936304):
   median **1.61→1.15×**, Σ **1.65→1.27×**; first Rust wall wins on dense shards.
3. **Profile remaining slower shards** — `01_chr21_w11`, `00_chr20_w26/w29`
   (200 kb heads locally; full 1 Mb is CI).
4. ~~**PairHMM prefix-vs-pack A/B**~~ — frac + min-haps A/B **reverted** (occupancy
   unchanged); next = read-axis / wavefront with TRACE.
5. **ci-subset** — F1 safety only (Peak sequential); run [32263139535](https://github.com/SynapticFour/gatk-rs/actions/runs/32263139535).
6. **Call-rate / L8** — spine strong-het gates + sync retain (see
   [`CALLRATE_EMIT_AD_DIG.md`](CALLRATE_EMIT_AD_DIG.md)).
7. **Publish hygiene** — stop `git push` to protected `main` in GIAB Publish jobs.

### Local TRACE rematch (phase10 tip → ledger tip)

| Window | assign_genotype_s | pairhmm_s | realign_s |
|--------|------------------:|----------:|----------:|
| mega `21:9825–9828k` (pre-G1) | ~147–172 | ~15–20 | ~4–7 |
| mega (post-G1 tip) | **~1.2** | ~13 | ~4.7 |
| w09 200 kb (tip profile) | ~3.2 | ~22 | ~22 |
| w11 densest-ish 50 kb (phase10) | ~28–32 | ~31–36 | ~14–18 |

Laptop assign noise is high; prefer CI wall-losers for signed wall ratios.

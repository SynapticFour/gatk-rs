# Phase 7 — profile-driven wall cuts

## Profile (macOS `sample`, w09 200kb, thr=2, FastestAvailable)

Hottest frames (40s mid-run):

1. `likelihood_engine::score_read_against_haplotypes`
2. `read_threading_assembler::build_threading_graph_core` (+ hash / RT graph)
3. `smith_waterman::align_uppercase_ready` (realign)
4. `pairhmm_simd::pack::score_haps_logless_packed_f64`
5. `event_map::variation_events_for_haplotype` / `EventMap::from_haplotype_and_reference`
6. `try_genotype_variation_event`

## Cuts in this phase

| Area | Change | Contract |
|------|--------|----------|
| Allele filter | `PerHaplotypeVariationEvents` only on HC-inverse PL ranking (not mark path) | Same support truth; avoid mark-path precache tax |
| PairHMM NEON/AVX2 | Leftover hap scores reuse prebuilt transitions (`*_with_transitions`) | Same Logless numerics; SIMD unit gate |
| Fair wall CI | `GIAB_MODE=wall-losers` — 8 dense 1 Mb windows, `GATK_RS_HC_SEQUENTIAL=0` | Product wall lane; Peak abort kept |
| L8 | [`L8_HOLDOUT_F1_TRACK.md`](L8_HOLDOUT_F1_TRACK.md) | Separate parity track |

## Deferred (next Instruments leaf)

- ~~SW / realign (`align_uppercase_ready`)~~ → [`PHASE8_SW_REALIGN.md`](PHASE8_SW_REALIGN.md)
- RT graph hash cost — assemble pie; keep from regressing k-best win on w11

## Prove

- `cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test`
- L2 / hap.py smoke on PR
- Local TRACE rematch w09 shares vs phase6
- Optional: `gh workflow run giab-genomewide.yml -f mode=wall-losers`

# Beat Java wall — hard look (phase6)

Goal: product-shaped wall (thr=2, no sequential) **&lt;1.0×** Java on dense GIAB, with
algorithm parity (no P12 band widening). Peak RSS already wins on probe.

## What the evidence says now

1. **CI dense median ~1.57×** (cancelled phase5 run) — flat vs phase4. Ownership /
   TLS / genotype-cache work did not move the CI median.
2. **Local TRACE:** w09 shares flat; **w11 k-best share dropped ~19%→8%** (real
   win), but PairHMM + genotype + realign still dominate (~38%+30%+14%). w26 flat.
3. **CI wall ≠ product wall** — `giab-genomewide.yml` sets `GATK_RS_HC_SEQUENTIAL=1`
   for Peak. Treat CI shard times as Peak-constrained; use local TRACE / a future
   non-sequential wall lane for beat-Java claims.
4. **PairHMM default is already `FASTEST_AVAILABLE`** in code. Remaining gate is
   signed dense/holdout F1, not another default flip.
5. **Java probe Peak was truncated** on 02_probes; probe + summarizer fixes land in
   phase6 so the next ci-subset can quote honest RSS.

## Highest-leverage wall bets (evidence-class, not locus pins)

| Priority | Bet | Why | Gate |
|----------|-----|-----|------|
| 1 | **PairHMM region wall** | Still #1 share on all losers; NEON_F64 already on path | TRACE rematch; SIMD unit; F1 holdout |
| 2 | **Genotype / EventMap** | Rose to ~30% share on w11; prior ownership cuts insufficient | TRACE share drop; w09 call-set identity |
| 3 | **Realign** | Stable ~12–15% | TRACE; call-set identity |
| 4 | **k-best / RT extract** | w11 already improved; keep from regressing; extract-store still visible | TRACE |
| 5 | **Fair wall CI lane** | Separate job **without** sequential, thr=2 | Compare to Java same lane |
| 6 | **L8 holdout F1** | Separate parity track (rust F1≈0.84) | `run_l9_signoff_gates.sh` |

## Reject / defer

- Widening P12 bands or deleting pins to “win” wall.
- Nested Rayon on PairHMM (already collapsed to one axis).
- Full EventMap precache in allele filtering (regressed prep wall badly).
- Marketing genome-wide equivalence from Peak-mode CI Σ.

## Immediate next experiment

1. ~~Finish phase6 TRACE / probe fix / ci-subset~~ (landed #112).
2. ~~Profile + PairHMM/EventMap + wall-losers lane~~ (landed #113).
3. ~~wall-losers Finalize dig~~ — concat bug, not F1 ([`WALL_LOSERS_F1_DIG.md`](WALL_LOSERS_F1_DIG.md)).
4. ~~SW / realign leaf~~ — [`PHASE8_SW_REALIGN.md`](PHASE8_SW_REALIGN.md).
5. Dispatch `GIAB_MODE=ci-subset` on phase8 PR (Peak rematch).
6. Re-run `wall-losers` after concat fix (product wall + honest F1 delta).
7. L8 holdout F1 — [`L8_HOLDOUT_F1_TRACK.md`](L8_HOLDOUT_F1_TRACK.md).
8. Next leaf: RT graph hash / deeper PairHMM after phase8 rematch.

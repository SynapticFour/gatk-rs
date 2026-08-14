# Phase 5 wall TRACE baseline (before Layers A–E cuts)

Frozen from phase4 post-cut / prior loser TRACE. Product shape: thr=2, no `GATK_RS_HC_SEQUENTIAL`.
PairHMM default remains as configured; campaign wall path uses `--pair-hmm FASTEST_AVAILABLE`.

## Windows

| Label | Notes |
|-------|-------|
| w09 200kb | `21:9500000-9700000`; post phase4 cuts |
| w11 head | ~200kb head of chr21_w11 BAM |
| w26 head | ~200kb head of chr20_w26 BAM |

## Top phases (sum delta_ms)

### w09 200kb postcuts (phase4)
- wall ~72s; `kbest_begin` 14.6s; `after_pairhmm` 17.6s; `prep_realign` 8.3s; `after_genotype` 5.8s; Σ~67s

### w11 head
- `after_pairhmm` 269.6s; `after_genotype` 200.4s; `kbest_begin` 156.3s; `prep_realign` 92.7s; Σ~810s

### w26 head
- `after_pairhmm` 50.5s; `after_genotype` 25.2s; `prep_realign` 19.1s; Σ~127s

## CI prelim (phase4 run 31744034963, honest dense)
- median wall Rust/Java ~1.57×; probe Peak RSS Rust/Java ~0.37×
- losers: w09 ~2.45×, w11 ~2.31×, w26 ~2.10×, w29 ~2.57×

# Phase 6 wall / mem baseline (post #110+#111)

Reassessment after phase5 merge and hot-path alloc PR (#111). Product wall shape:
`RAYON_NUM_THREADS=2`, **no** `GATK_RS_HC_SEQUENTIAL`, `--pair-hmm FASTEST_AVAILABLE`.

## CI-subset run `31772652099` (phase5 branch; **cancelled** before finalize)

Artifacts retained (~112/114 paired HC shards). Honest dense filter (exclude
near-empty Rust &lt;2s walls):

| Metric | Value |
|--------|------:|
| Dense median wall Rust/Java | **~1.57×** (unchanged vs phase4 prelim) |
| Σ wall (all paired, includes empties) | ~1.46× |
| Probe Peak RSS Rust | **~60.6 MiB** (full-ish coverage of 2m23s) |
| Probe Peak RSS Java | **≥951.6 MiB** (**truncated** — log stops ~15s after JVM appears; wall ~1m42s) |

Campaign losers (wall Rust/Java from `/usr/bin/time`; RSS columns in summarizer tables
are still often docker-time junk for Java — use probe):

| Shard | Wall R/J |
|-------|---------:|
| 01_chr21_w09 | 2.42× |
| 01_chr21_w11 | 2.39× |
| 00_chr20_w26 | 2.15× |
| 00_chr20_w11 | 1.64× |
| 00_chr20_w29 | 1.60× |
| 00_chr20_w09 | 1.41× |
| 01_chr21_w26 | 1.65× |
| 01_chr21_w29 | 1.28× |

**Do not** treat CI Σ or CI shard wall as product wall: workflow sets
`GATK_RS_HC_SEQUENTIAL=1` (Peak-RSS mode). Local TRACE is the wall pie.

## Probe lifetime fix

`hc_mem_probe.sh` previously used `pgrep -n` (newest match), which often pinned the
`docker run` wrapper (~30 MiB) and could miss JVM Peak; cancelled jobs also truncated
mid-climb. Phase6 fix: pick **max VmRSS** among HC-like PIDs (skip docker when a
heavier child exists) and emit a **final sample on TERM/INT/EXIT**.

Summarizer now prefers `probe_peak_rss_kb` over `/usr/bin/time` RSS.

## Local TRACE rematch (`phase6-rematch/`, binary #111)

Phenotype **counts** match phase4/phase5 on w09/w11 (same work shape). Absolute
seconds drift with host load — compare **shares**.

### w09 200kb shares

| Phase | phase4 Σ share | phase6 Σ share |
|-------|---------------:|---------------:|
| `after_pairhmm` | ~26% | ~28% |
| `kbest_begin` | ~22% | ~18% |
| `prep_realign` | ~12% | ~13% |
| `after_genotype` | ~9% | ~8% |
| `rt_extract_cache_store` | ~8% | ~9% |

### w11 head shares (vs phase5 baseline head)

| Phase | phase5 Σ share | phase6 Σ share |
|-------|---------------:|---------------:|
| `after_pairhmm` | ~33% | ~38% |
| `after_genotype` | ~25% | ~30% |
| `kbest_begin` | ~19% | **~8%** |
| `prep_realign` | ~12% | ~14% |
| `rt_extract_cache_store` | ~3% | ~4% |

### w26 head shares

Flat vs phase5 baseline (PairHMM ~37%, genotype ~19%, realign ~15%, RT extract
~11%). Counts identical.

**Verdict:** k-best share improved on w11; wall still dominated by PairHMM +
genotype + realign. No P12 widening. w29 not rematched locally (no loser BAM pin).

## PairHMM default / F1

Code default is already [`PairHmmImpl::FastestAvailable`](../../gatk-haplotypecaller/src/likelihood_engine.rs)
(host SIMD when present). Docs that still say “production remains Log10” are stale —
promotion is **signed F1**, not a further default flip.

- Unit: `pairhmm_simd_vs_scalar_test` (7/7) green on this tree.
- Holdout regen (`20:15000000-15050000`, current release + default FastestAvailable):
  java=55 rust=41 shared=40; **rust F1≈0.844** vs java F1=1.0 — **gate fail**
  (`parity/reports/hc-full-parity-j6-dense-holdout/phase6_fastest_f1.*`).
- Dense holdout F1 remains a **separate L8 parity track** (do not conflate with wall;
  do not widen P12).

## L8 holdout

Keep L8 F1 / FORMAT holdout on its own track. Wall PRs must not “fix” holdout by
widening bands or deleting P12 pins.

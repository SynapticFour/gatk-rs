# gatk-rs documentation

Every document under `docs/`, grouped as a tutorial, how-to guides, reference, and explanation. One line each.

Read the claim matrix before citing a result. Everything else is evidence, history, or a specialist note. Product claims are only the **Yes** rows in [CLAIM_MATRIX.md](CLAIM_MATRIX.md).

Notes under `parity/` are 6R/7R investigation records and k-best policy. They are not acceptance by themselves. `perf/` holds performance notes. A public speed or memory claim still requires the dedicated host in [ci/PERF_BENCHMARK_HOST.md](ci/PERF_BENCHMARK_HOST.md). [EQUIVALENCE_DASHBOARD.md](EQUIVALENCE_DASHBOARD.md) has an empty F1 table, so that file is not a measured result. `parity-site/data/history.json` contains historical runs. Those rows are not a current-main **Yes**.

## Tutorial

| Document | What it is |
|----------|------------|
| [GETTING-STARTED.md](GETTING-STARTED.md) | Getting started |
| [USER_GUIDE.md](USER_GUIDE.md) | User guide |

## How-to guides

| Document | What it is |
|----------|------------|
| [PILOT_GUIDE.md](PILOT_GUIDE.md) | Pilot guide — verify gatk-rs on your own data |
| [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) | Release checklist |
| [ci/PERF_BENCHMARK_HOST.md](ci/PERF_BENCHMARK_HOST.md) | Dedicated performance-benchmark host |
| [perf/AIR_M4_GIAB_RECIPE.md](perf/AIR_M4_GIAB_RECIPE.md) | Air M4 16 GB GIAB recipe (memory + disk) |
| [perf/BEAT_JAVA_WALL_NEXT.md](perf/BEAT_JAVA_WALL_NEXT.md) | Beat Java wall — hard look (phase6→10 + ledger tip) |
| [perf/BENCH_PROFILE_MATRIX.md](perf/BENCH_PROFILE_MATRIX.md) | Cargo bench profile audit |
| [perf/CALLRATE_EMIT_AD_DIG.md](perf/CALLRATE_EMIT_AD_DIG.md) | Call-rate dig — ci-subset ΔF1 + L8 holdout (post hang-fix) |
| [perf/CI_SUBSET_HANG_W47_W10.md](perf/CI_SUBSET_HANG_W47_W10.md) | ci-subset hang: `00_chr20_w47` / `01_chr21_w10` rust — **fixed** |
| [perf/DEDICATED_BENCH_LATEST.md](perf/DEDICATED_BENCH_LATEST.md) | Latest dedicated-host benchmark pointer |
| [perf/DENSE_LOSER_TRACE_RECIPE.md](perf/DENSE_LOSER_TRACE_RECIPE.md) | Dense loser TRACE recipe (phase5 wall campaign) |
| [perf/FAIR_HC_COMPARISON.md](perf/FAIR_HC_COMPARISON.md) | Fair HaplotypeCaller comparison (dedicated host) |
| [perf/GENOTYPE_ASSIGN_COMPLEXITY.md](perf/GENOTYPE_ASSIGN_COMPLEXITY.md) | Genotype assignment complexity & dense-region wall |
| [perf/HC_MEMORY_BASELINE_20260804.md](perf/HC_MEMORY_BASELINE_20260804.md) | HC memory baseline — pre Arc-share pass (2026-08-04) |
| [perf/HC_MEMORY_PROFILE.md](perf/HC_MEMORY_PROFILE.md) | HaplotypeCaller memory profile (reproducible) |
| [perf/HC_PRODUCTION_PROFILE.md](perf/HC_PRODUCTION_PROFILE.md) | Production HaplotypeCaller profiling |
| [perf/HOST_SPECS.md](perf/HOST_SPECS.md) | Performance host specifications |
| [perf/L8_HOLDOUT_F1_TRACK.md](perf/L8_HOLDOUT_F1_TRACK.md) | L8 holdout F1 track (separate from wall) |
| [perf/PAIRHMM_SPEEDUP.md](perf/PAIRHMM_SPEEDUP.md) | PairHMM speedup profile (reproducible) |
| [perf/PAIRHMM_WAVEFRONT.md](perf/PAIRHMM_WAVEFRONT.md) | PairHMM wavefront backend (Architecture C, phase-1) |
| [perf/PERFORMANCE_LEDGER.md](perf/PERFORMANCE_LEDGER.md) | HaplotypeCaller performance ledger |
| [perf/PHASE5_WALL_BASELINE.md](perf/PHASE5_WALL_BASELINE.md) | Phase 5 wall TRACE baseline (before Layers A–E cuts) |
| [perf/PHASE6_WALL_BASELINE.md](perf/PHASE6_WALL_BASELINE.md) | Phase 6 wall / mem baseline (post #110+#111) |
| [perf/PHASE7_WALL_PROFILE.md](perf/PHASE7_WALL_PROFILE.md) | Phase 7 — profile-driven wall cuts |
| [perf/PHASE8_SW_REALIGN.md](perf/PHASE8_SW_REALIGN.md) | Phase 8 — SW / realign leaf |
| [perf/PHASE9_PAIRHMM_PACK.md](perf/PHASE9_PAIRHMM_PACK.md) | Phase 9 — PairHMM pack occupancy + hap prefix reuse (product wall) |
| [perf/RT_KMER_PACKED_KEYS.md](perf/RT_KMER_PACKED_KEYS.md) | Read-threading k-mer representation (packed keys) |
| [perf/RUST_MEMORY_CACHE_AUDIT.md](perf/RUST_MEMORY_CACHE_AUDIT.md) | Rust memory / cache audit — HaplotypeCaller |
| [perf/RUST_SHOWCASE_ROADMAP.md](perf/RUST_SHOWCASE_ROADMAP.md) | Rust-native HC performance showcase (parity-safe) |
| [perf/SW_HC_OPTIMIZE.md](perf/SW_HC_OPTIMIZE.md) | Smith-Waterman for production HaplotypeCaller |
| [perf/WALL_LOSERS_F1_DIG.md](perf/WALL_LOSERS_F1_DIG.md) | wall-losers Finalize dig (run 31857571738) |
| [perf/hc_memory_profile_latest.json](perf/hc_memory_profile_latest.json) | Recorded data. Not a current-main Yes row by itself. |
| [perf/pairhmm_speedup_latest.json](perf/pairhmm_speedup_latest.json) | Recorded data. Not a current-main Yes row by itself. |

## Reference

| Document | What it is |
|----------|------------|
| [CLAIM_MATRIX.md](CLAIM_MATRIX.md) | Claim matrix — what gatk-rs HaplotypeCaller asserts |
| [EQUIVALENCE_DASHBOARD.md](EQUIVALENCE_DASHBOARD.md) | Equivalence dashboard (nightly trio E2E) |
| [FOR-EVALUATORS.md](FOR-EVALUATORS.md) | For evaluators |
| [GATK_PINNED.env](GATK_PINNED.env) | Recorded data. Not a current-main Yes row by itself. |
| [PARITY.md](PARITY.md) | Parity and validation |
| [PARITY_MILESTONE_6R.md](PARITY_MILESTONE_6R.md) | 6R milestone — canonical mid-B HaplotypeCaller path |
| [SUPPORTED_FEATURE_MATRIX.md](SUPPORTED_FEATURE_MATRIX.md) | Supported feature matrix — HaplotypeCaller |
| [equivalence/baseline.json](equivalence/baseline.json) | Recorded data. Not a current-main Yes row by itself. |
| [parity-site/data/history.json](parity-site/data/history.json) | Recorded data. Not a current-main Yes row by itself. |
| [parity-site/data/latest.json](parity-site/data/latest.json) | Recorded data. Not a current-main Yes row by itself. |
| [parity-site/data/perf_history.json](parity-site/data/perf_history.json) | Recorded data. Not a current-main Yes row by itself. |
| [parity/6R.107_FINAL_CLOSURE.md](parity/6R.107_FINAL_CLOSURE.md) | 6R.107 — Final closure / verification |
| [parity/6R.107_TEST_MAINTENANCE.md](parity/6R.107_TEST_MAINTENANCE.md) | 6R.107 test maintenance — HOLDOUT_6R104 / 6R.105 |
| [parity/6R.312_FINAL_CLOSURE.md](parity/6R.312_FINAL_CLOSURE.md) | 6R.312 — Frozen canonical Java-parity checkpoint |
| [parity/6R.314_HOLDOUT.md](parity/6R.314_HOLDOUT.md) | 6R.314 — Post-milestone generalization holdouts |
| [parity/6R.43_HOLDOUT_MATRIX.md](parity/6R.43_HOLDOUT_MATRIX.md) | 6R.43 holdout matrix |
| [parity/6R.45_RETENTION.md](parity/6R.45_RETENTION.md) | 6R.45 — `findBestPaths` haplotype retention |
| [parity/6R.46_TRIM.md](parity/6R.46_TRIM.md) | 6R.46 — haplotype trim vs `p12_snp_cluster` |
| [parity/6R.47_EVENTMAP.md](parity/6R.47_EVENTMAP.md) | 6R.47 — EventMap: no 40 bp CIGAR allele cap |
| [parity/6R.48_QD_RNG.md](parity/6R.48_QD_RNG.md) | 6R.48 — QualByDepth RNG stream (process-global Java `Random`) |
| [parity/7R.22_JAVA_EQUIVALENT_KBEST_DEFAULT.md](parity/7R.22_JAVA_EQUIVALENT_KBEST_DEFAULT.md) | 7R.22 — Java-equivalent unbounded K-best is the production default |
| [parity/7R.26_NORMALIZE_ALLELES_PRODUCTION_FIX.md](parity/7R.26_NORMALIZE_ALLELES_PRODUCTION_FIX.md) | 7R.26 — Production fix for normalize_alleles caller ownership |
| [parity/KBEST_RESOURCE_POLICY.md](parity/KBEST_RESOURCE_POLICY.md) | SeqGraph k-best resource policy (6R.140, production default 7R.22) |
| [parity/README.md](parity/README.md) | Index of the parity investigation archive. Not a Yes row. |
| [parity/6R.141_KBEST_ORACLE.tsv](parity/6R.141_KBEST_ORACLE.tsv) | Recorded oracle table. Not a current-main Yes row by itself. |
| [parity/6R.146_FILTER_ORACLE.tsv](parity/6R.146_FILTER_ORACLE.tsv) | Recorded oracle table. Not a current-main Yes row by itself. |
| [parity/6R.148_ORACLE.tsv](parity/6R.148_ORACLE.tsv) | Recorded oracle table. Not a current-main Yes row by itself. |
| [parity-site/README.md](parity-site/README.md) | Static site of historical equivalence and performance tabs. Not a current-main Yes. |
| [parity-site/index.html](parity-site/index.html) | Page of that static site. |
| [parity-site/app.js](parity-site/app.js) | Script of that static site. |
| [parity-site/style.css](parity-site/style.css) | Style of that static site. |
| [parity-site/.nojekyll](parity-site/.nojekyll) | Tells GitHub Pages not to run Jekyll on that site. |

## Explanation

| Document | What it is |
|----------|------------|
| [ARCHITECTURE.md](ARCHITECTURE.md) | Architecture |
| [adr/0001-scope-boundary.md](adr/0001-scope-boundary.md) | ADR 0001: Scope boundary — germline short-variant workflow, not a full GATK4 port |
| [adr/0002-remove-gatk-tools.md](adr/0002-remove-gatk-tools.md) | ADR 0002: Remove the `gatk-tools` stub crate |

Files outside `docs/` that the same reading needs: [CHANGELOG.md](../CHANGELOG.md), [SECURITY.md](../SECURITY.md), [CONTRIBUTING.md](../CONTRIBUTING.md), [LICENSE](../LICENSE), [NOTICE.md](../NOTICE.md), [CITATION.cff](../CITATION.cff).

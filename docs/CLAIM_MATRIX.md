# Claim matrix — what gatk-rs HaplotypeCaller asserts

**Authority for product claims.** If a claim is not **Yes** here, do not assert it.  
**Pinned Java oracle:** GATK **4.4.0.0** — see [`GATK_PINNED.env`](GATK_PINNED.env) and root `GATK_PINNED_SHA`.  
**Last updated:** 2026-07-24 (L2 green; Combine / Genotype / VariantFiltration mini parity green)

---

## Asserted (Yes)

| Claim | Scope | Evidence |
|-------|-------|----------|
| Production path is `CallRegionArgs::strict_java()` + `assembly-region-v1` | Default CLI HC | `gatk-haplotypecaller` `run.rs` |
| P12 L3 variant emit parity | `chr2:92300000–92350000`, 66 Java sites, `rust_only=0` | P12 L3 parity scripts / CI |
| P12 L4 FORMAT parity (algorithmic) | Same interval, 66 sites | P12 L4 parity scripts / CI |
| P12 L5 gVCF block parity | Same interval | P12 L5 parity scripts / CI |
| L2 synthetic gates | **223/223 strict** on `2026-07-24T07:00:10Z` (`parity/reports/hc-full-parity-l2/hc_full_parity_l2_20260724T065936Z.log` → `hc_full_parity_l2_canonical.log` / `last_run.log`); `l2_summary.json` 223×`equal=true` | `scripts/parity/run_hc_full_parity_l2.sh` (full-suite run updates `last_run.log`; canonical only on strict green). Permanent isolated gates for the historical May failures: `cargo test -p gatk-haplotypecaller --features parity_harness --test l2_may_regression_gate_test` (`e0-assemble/p5_case1_assemble`, `e2e/p5_indel_chrindel`). Note: the May-17 `last_run.log` (`passed=143 failed=2`) was **older** than the Jul-21 green canonical, not a post-Jul-21 regression. |
| L6 scale + GIAB truth gate (spine parity + F1 tracks Java) | Eval `2:92000000-92400000`; parity on spine `92300000-92350000` | L6 sign-off gates (archive branch) |
| L7 dense GIAB + second non-chr2 slice + GT FORMAT hard gate | chr20 dense F1≥0.95; chr21 dense; GT mismatch ≤15% among TPs | L7 sign-off gates (archive branch) |
| L9–L14 trajectory closure on holdout / soft-PL policy | Dense/holdout F1 gates; P12 66/66; soft PL permanent residual | L9–L14 sign-off gates (archive branch) |
| Cluster indel phenotype (Δ=+3) | Algorithm + non-P12 tests; P12 oracle fallback | HC tests + parity harness |
| Oracle TSV does not gate production emit | Emit admission / sparse rescue | `scripts/parity` oracle audits |
| `parity_aligned` / legacy bridges off release surface | Needs `cfg(test)` or `--features parity_harness` | `gatk-haplotypecaller` Cargo features |
| GIAB multi-sample equivalence **harness** | Infrastructure only — does **not** by itself assert callset equivalence | [`scripts/parity/giab/`](../scripts/parity/giab/README.md), `.github/workflows/giab-genomewide.yml` |
| CombineGVCFs mini REF/ALT/PL parity (incl. different ALT sets → diploid PL remap) | Synthetic 2-sample mini cohort; site `chr1:10` ALT `T,G,<NON_REF>`; SAMPLE1 PL `100,100,100,0,…` / SAMPLE2 PL `90,0,90,…` | `parity/reports/combine_gvcfs_20260724T072518Z.log` (`OK sites=5`); unit gates `ref_confidence_merger::pl_remap_tests` + `combine_gvcfs::tests::t04` |
| GenotypeGVCFs mini alleles/GT/QUAL parity | Same mini cohort after CombineGVCFs; QUAL ±20.0 | `parity/reports/genotype_gvcfs_20260724T072526Z.log` (`OK sites=1`) |
| VariantFiltration boundary FILTER parity | Synthetic SNP hard-filter boundary sites | `parity/reports/variant_filtration_20260724T072535Z.log` (`OK sites=16` identical FILTER) |

---

## Pending green gate (do not market until a signed run exists)

| Claim | Required evidence | Status |
|-------|-------------------|--------|
| GIAB **ci-subset** equivalence (HG001 at minimum): \|Rust−Java\| F1 Δ ≤ threshold via `gatk-rs-equiv` | Weekly workflow artifact + dashboard | **Not yet signed** |
| GIAB **ci-subset** on HG001+HG002+HG005 | Same, all three samples | **Not yet signed** |
| GIAB **full autosomes** (chr1–22) equivalence | `GIAB_MODE=autosomes` green run | **Not yet signed** |
| Nightly / GIAB Pages dashboard (`docs/EQUIVALENCE_DASHBOARD.md`, `docs/parity-site/data/history.json`) | Successful `nightly-equivalence.yml` and/or `giab-genomewide.yml` publish with non-empty `history.json` `runs` | **Not yet signed** — `runs: []` (no remote CI yet) |

---

## Scoped / waived (Yes, with limits)

| Claim | Limit | Waiver |
|-------|-------|--------|
| Cluster materialize / mapper hooks on `strict_java()` | P12 interval overlap only | **W-H1** |
| gVCF RCM band reconciliation | P12 interval only — **permanent** | **W-H3** |
| Soft-clip tier-3 FORMAT shaping | Named evidence thresholds (phenotype) | **W-J4-band** |
| Hom-ref deserts / dense RCM / gradation tables | Not universalized — permanent scoped tables | **W-J1 / W-J3 / W-J5** |
| Dense FORMAT soft fields (AD/DP) + residual PL | GT **0%**, soft AD/DP ≤0.30; permanent soft-PL residual | **W-L7-FORMAT** (permanent) |

---

## Not asserted (No)

| Claim | Reality |
|-------|---------|
| Genome-wide (full autosomes) GATK 4.4 HaplotypeCaller equivalence | **No** — signed evidence remains P12 + L2 + dense/holdout windows |
| GIAB `ci-subset` / multi-sample truth equivalence as a product claim | **No (not signed yet)** |
| Multi-sample joint HC (Java merges `-I` reads) | **No** — each BAM traversed independently |
| Bitwise-identical QUAL/FORMAT genome-wide | **No** — L4 is P12 66-site lock |
| Full product feature parity (bamout, DRAGSTR, DRAGEN, `AS_*`, `--assembly-region-out`) | **No** — deferred / scaffold only |
| VQSR (VariantRecalibrator / ApplyVQSR) | **No** — use `VariantFiltration` hard filters; not an algorithmic VQSR substitute (GATK-aligned small-cohort fallback) |
| Genome-wide L6 / clinical truth performance | **No** — local windows only until a signed GIAB gate lands |
| Java class / source-shape parity as a product goal | **No** — target is **algorithm** parity |
| `gatk-tools` generic toolkit | **No** — crate **removed** ([ADR 0002](adr/0002-remove-gatk-tools.md)); use samtools/bcftools |
| Full GATK4 port (BQSR, VQSR, Mutect2, gCNV/SV, Funcotator, …) | **No** — intentional scope boundary ([ADR 0001](adr/0001-scope-boundary.md)) |

---

## Deferred CLI / product features (summary)

Registry IDs (T3-5 … T5-6) are kept for CI audits; they are **not** product commitments.
Former **T0-4** (`gatk-tools`) was removed — see [ADR 0002](adr/0002-remove-gatk-tools.md).
Product scope: [ADR 0001](adr/0001-scope-boundary.md).

| ID | Feature | Production status |
|----|---------|-------------------|
| T3-5 | `-alleles` / given alleles | CLI wired; broader sign-off optional beyond L2 `c5-force` |
| T5-1 | bamout / HaplotypeBAMWriter | Deferred |
| T5-2 | DRAGSTR calibration | Deferred (scaffolds / dumps only) |
| T5-3 | DRAGEN mode | Deferred (scaffolds / dumps only) |
| T5-4 | `AS_*` annotations | Deferred |
| T5-5 | DP reconciliation | Deferred |
| T5-6 | `--assembly-region-out` | Deferred — use `DumpSmoothedActivity` for analysis |
| — | VQSR | Not implemented — `VariantFiltration` hard filters only |

---

## How to cite

- External / README: link this file + waivers W-H1 / W-H3 / W-L7-FORMAT; keep the Alpha / non-Broad disclaimer.
- Code reviews: reject PRs that claim “full GATK parity”, “100% CLI compatibility”, or “bitwise identical outputs” without updating this matrix.
- Historical L2–L14 sign-off narratives live on branch `pre-cleanup-archive`, not on `main`.

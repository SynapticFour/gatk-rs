# Claim matrix — what gatk-rs HaplotypeCaller asserts

**Authority for product claims.** If a claim is not **Yes** here with evidence
reachable from **`main`**, do not assert it.

**Pinned Java oracle:** GATK **4.4.0.0** — see [`GATK_PINNED.env`](GATK_PINNED.env) and root `GATK_PINNED_SHA`.

**Last updated:** 2026-08-15 (ci-subset HG001 ΔF1 gate signed on `main` via Finalize
`gate_passed=true`; Peak RSS already wins; product wall still ~1.79× Java)

---

## Asserted (Yes) — evidence on `main`

| Claim | Scope | Evidence on `main` |
|-------|-------|--------------------|
| Production path is `CallRegionArgs::strict_java()` + `assembly-region-v1` | Default CLI HC | `gatk-haplotypecaller/src/run.rs` |
| P12 L3 variant emit parity | `chr2:92300000–92350000`, 66 Java sites, `rust_only=0` | `scripts/parity/run_p12_l3_signoff.sh` (+ P12 tests under `gatk-haplotypecaller/tests/`) |
| P12 L4 FORMAT parity (algorithmic) | Same interval, 66 sites | `scripts/parity/run_p12_l4_signoff.sh` |
| P12 L5 gVCF block parity | Same interval | `scripts/parity/run_p12_l5_gvcf.sh` |
| L2 synthetic gates | **223/223 strict** on `2026-07-24T07:00:10Z` | `parity/reports/hc-full-parity-l2/hc_full_parity_l2_canonical.log` (+ `last_run.log`); `scripts/parity/run_hc_full_parity_l2.sh`; permanent May gates: `cargo test -p gatk-haplotypecaller --features parity_harness --test l2_may_regression_gate_test` |
| Cluster indel phenotype (Δ=+3) | Algorithm + non-P12 tests; P12 oracle fallback | HC tests + parity harness |
| Oracle TSV does not gate production emit | Emit admission / sparse rescue | `scripts/parity` oracle audits |
| `parity_aligned` / legacy bridges off release surface | Needs `cfg(test)` or `--features parity_harness` | `gatk-haplotypecaller` Cargo features |
| GIAB multi-sample equivalence **harness** (prepare → HC matrix → finalize + RTG/`gatk-rs-equiv`) | Infrastructure proven end-to-end on hosted CI; does **not** by itself assert genome-adjacent / autosome callset equivalence | [`scripts/parity/giab/`](../scripts/parity/giab/README.md), `.github/workflows/giab-genomewide.yml`; smoke green: [run 30703069224](https://github.com/SynapticFour/gatk-rs/actions/runs/30703069224) (`gate_passed=true`, `max_\|ΔF1\|=0` on three ~50 kb windows) |
| GIAB **ci-subset** equivalence (HG001): \|Rust−Java\| F1 Δ ≤ 0.02 via `gatk-rs-equiv` | chr20/21 windowed ci-subset concat; **not** full autosomes / multi-sample | Finalize **PASS** `max_\|ΔF1\|=0.0099` on [31903578250](https://github.com/SynapticFour/gatk-rs/actions/runs/31903578250) (`main` post-#116) and [31884483008](https://github.com/SynapticFour/gatk-rs/actions/runs/31884483008); site rate ~0.875. Workflow may still red on **Publish/Deploy** gh-pages push only — cite Finalize gate, not the Actions conclusion bubble. |
| CombineGVCFs mini REF/ALT/PL parity (incl. different ALT sets → diploid PL remap) | Synthetic 2-sample mini cohort; site `chr1:10` ALT `T,G,<NON_REF>` | `parity/reports/combine_gvcfs_20260724T072518Z.log` (`OK sites=5`); unit gates `ref_confidence_merger::pl_remap_tests` + `combine_gvcfs::tests::t04` |
| GenotypeGVCFs mini alleles/GT/QUAL parity | Same mini cohort after CombineGVCFs; QUAL ±20.0 | `parity/reports/genotype_gvcfs_20260724T072526Z.log` (`OK sites=1`) |
| CombineGVCFs → GenotypeGVCFs **cohort scale** (synthetic ladder) | Synthetic N∈{2,10,25,50,100} on chr1 10 kb / 400 SNPs; **recommended ≤ 100 samples** on this gate. Above 100 **untested / not claimed**. No GenomicsDBImport. | `parity/reports/joint_cohort_scale_20260724T184447Z/`; `scripts/parity/run_joint_cohort_scale.sh` |
| VariantFiltration boundary FILTER parity | Synthetic SNP hard-filter boundary sites | `parity/reports/variant_filtration_20260724T072535Z.log` (`OK sites=16` identical FILTER) |

---

## Historical (archive branch only — not “Yes” on `main`)

These narrative gates were signed on branch `pre-cleanup-archive`. They are **not**
reproducible as first-class product claims from current `main` until their
evidence/scripts are restored here. Cite them as historical engineering notes only.

| Claim | Former scope | Where to look |
|-------|--------------|---------------|
| L6 scale + GIAB truth gate | Eval `2:92000000-92400000`; spine `92300000-92350000` | `pre-cleanup-archive` L6 sign-off |
| L7 dense GIAB + second non-chr2 slice + GT FORMAT hard gate | chr20/chr21 dense F1 / GT mismatch bands | `pre-cleanup-archive` L7 sign-off |
| L9–L14 trajectory closure | Dense/holdout F1; P12 66/66; soft-PL residual policy | `pre-cleanup-archive` L9–L14 sign-off |

---

## Pending green gate (do not market until a signed run exists)

| Claim | Required evidence | Status |
|-------|-------------------|--------|
| GIAB **ci-subset** equivalence (HG001 at minimum): \|Rust−Java\| F1 Δ ≤ threshold via `gatk-rs-equiv` | Green `giab-genomewide.yml` with `GIAB_MODE=ci-subset`, artifact + dashboard row | **Signed (HG001)** — see Asserted ([31903578250](https://github.com/SynapticFour/gatk-rs/actions/runs/31903578250)). Remaining: Pages publish hygiene. |
| GIAB **ci-subset** on HG001+HG002+HG005 | Same, all three samples | **Not yet signed** |
| GIAB **full autosomes** (chr1–22) equivalence | `GIAB_MODE=autosomes` green run | **Not yet signed** |
| Nightly / GIAB Pages dashboard (`docs/EQUIVALENCE_DASHBOARD.md`, `docs/parity-site/data/history.json`) | Green `ci-subset`+ finalize → `publish-parity-site` job writing non-empty `history.json` `runs` (smoke never publishes) | **Not yet signed** — wire exists; waits on first successful non-smoke publish |
| GIAB **smoke** as a product / genome-adjacent claim | Smoke is PR/infra hygiene (hybrid P12 from NA12878_20k + chr20/21 30×): three ~50 kb windows only | **Not a product claim** — harness is green ([30703069224](https://github.com/SynapticFour/gatk-rs/actions/runs/30703069224)); do **not** cite smoke as `ci-subset` or autosome equivalence |
| GIAB **full autosomes** on GitHub-hosted runners | Hosted `giab-genomewide.yml` **rejects** `autosomes` (disk + 6 h cap) | Use self-hosted [`.github/workflows/genomewide-validation.yml`](../.github/workflows/genomewide-validation.yml) (operator-provisioned runner; setup notes are internal) |

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
| Genome-wide (full autosomes) GATK 4.4 HaplotypeCaller equivalence | **No** — signed evidence on `main` is P12 + L2 + synthetic joint/filter minis |
| GIAB `ci-subset` / multi-sample truth equivalence as a product claim | **Partial** — HG001 ci-subset ΔF1 **signed**; HG002+HG005 and autosomes **not** |
| Multi-sample joint HC (Java merges `-I` reads) | **No** — each BAM traversed independently |
| CombineGVCFs / GenotypeGVCFs for **large cohorts** (WGS × N≫100, GenomicsDB-class) | **No** — gatk-rs Combine loads full gVCFs in memory and has no GenomicsDBImport path. Signed synthetic scale gate: **recommended ≤ 100 samples** on 10 kb/400-SNP ladder |
| Bitwise-identical QUAL/FORMAT genome-wide | **No** — L4 is P12 66-site lock |
| Full product feature parity (bamout, DRAGSTR, DRAGEN, `AS_*`, `--assembly-region-out`) | **No** — deferred / scaffold only |
| VQSR (VariantRecalibrator / ApplyVQSR) | **No** — use `VariantFiltration` hard filters |
| Clinical / production drop-in replacement for Broad GATK | **No** — Alpha experiment; limited regions |
| Genome-wide L6 / clinical truth performance | **No** — archive L6 is not a `main` claim |
| Java class / source-shape parity as a product goal | **No** — target is **algorithm** parity |
| `gatk-tools` generic toolkit | **No** — crate **removed** ([ADR 0002](adr/0002-remove-gatk-tools.md)); use samtools/bcftools |
| Full GATK4 port (BQSR, VQSR, Mutect2, gCNV/SV, Funcotator, …) | **No** — intentional scope boundary ([ADR 0001](adr/0001-scope-boundary.md)); **not** a launch target |

---

## Deferred CLI / product features (summary)

Registry IDs (T3-5 … T5-6) are kept for CI audits; they are **not** product commitments.
Former **T0-4** (`gatk-tools`) was removed — see [ADR 0002](adr/0002-remove-gatk-tools.md).
Product scope: [ADR 0001](adr/0001-scope-boundary.md).

Utility CLI (`PrintReads`, `CountReadsInRegion`, `Validate`, …) remains callable for
parity harnesses but is **hidden from default `--help`** — prefer **samtools** /
**bcftools** for generic BAM/VCF ops.

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
| — | BQSR / Mutect2 / gCNV / Funcotator | Out of scope — do not stub toward a toolkit clone |

---

## How to cite

- External / README: link this file + waivers W-H1 / W-H3 / W-L7-FORMAT; keep the Alpha / non-Broad disclaimer.
- Code reviews: reject PRs that claim “full GATK parity”, “100% CLI compatibility”, or “bitwise identical outputs” without updating this matrix.
- Historical L6–L14 sign-off narratives live on branch `pre-cleanup-archive`, not as unqualified **Yes** rows on `main`.

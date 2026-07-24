# Parity scripts

**HC full parity:** status and work queue → [`docs/CLAIM_MATRIX.md`](../../docs/CLAIM_MATRIX.md). Governance refs: emitter matrix, fixture manifest, numeric contract in the same directory.

| Script | Purpose |
|--------|---------|
| `compare_semantic_trace.py` | Compare two `gatk_rs.hc.semantic_trace/v1` NDJSON traces; report first divergence. |
| `project_java_to_semantic_trace.py` | Project Java regions TSV + VCF into semantic-trace NDJSON (`impl=java`). |
| `lib_pinned_gatk.sh` | Load pinned GATK ref/SHA/docker defaults from `docs/GATK_PINNED.env`. |
| `verify_pinned_gatk.sh` | **P0:** verify Docker/`gatk`/`GATK_JAR` matches the pin (foundation gate `pinned-gatk-upstream`). |
| `run_hc_full_parity_java_compile.sh` | Compile `HcFullParityGateDump` against pinned GATK Docker jar. |
| `run_hc_full_parity_java_refresh.sh` | Regenerate `parity/fixtures/hc-full-parity/java_dumps/*_<pin>.tsv`. |
| `run_hc_full_parity_l2.sh` | **L2 strict (default):** Rust `hc_full_parity_gate_dump` vs frozen `java_dumps/` (`PARITY_HC_FULL_L2_STRICT=0` for advisory only). Required in CI as `hc-full-parity-l2-strict`. |
| `run_hc_full_parity_j6_truth.sh` | **L6 gate:** NA12878 P12 interval scale + GIAB stratified F1 (`parity/fixtures/hc-full-parity/j6/thresholds.json`). Weekly CI: `.github/workflows/p12-l6-scale.yml`. |
| `run_hc_full_parity_phase_l2.sh` | P1 bundle (probe Java dumps + L2). |
| `run_java_gatk.sh` | Run Java GATK (`GATK_JAR`, `gatk` on `PATH`, or `GATK_DOCKER_IMAGE`; defaults from pin). |
| `run_rust_gatk.sh` | Run `gatk-rs` (`PARITY_RUST_PROFILE=dev\|release`). |
| `run_parity_smoke.sh` | Deterministic smoke matrix (`LC_ALL`, `TZ`, `RUST_LOG`). |
| `run_foundation_gate.sh` | Required Phase 0/1 gate runner (`parity/checks.json`). |
| `run_read_filter_diff.sh` | Runtime Java-vs-Rust read-filter differential checks (synthetic SAM + BAM slice). |
| `run_read_header_semantics_diff.sh` | Runtime Java-vs-Rust SAM header semantics differential checks (`ValidateSamFile` vs `Validate`). |
| `run_p3_region_query_diff.sh` | Phase 3 region-query differential checks (`CountReads` vs `CountReadsInRegion`). |
| `run_p3_indexed_edge_query_diff.sh` | Phase 3 indexed edge-interval matrix (`CountReads` vs `CountReadsInRegion` on boundary/empty/overlap cases). |
| `run_p3_unmapped_supplementary_diff.sh` | Phase 3 supplementary + unmapped matrix (`CountReads` vs `CountReadsInRegion` on split intervals). |
| `run_p3_malformed_corpus_diff.sh` | Phase 3 malformed-corpus differential matrix (`ValidateSamFile`/`ValidateVariants` vs `Validate`). |
| `run_p3_truncation_corruption_diff.sh` | Phase 3 truncation/corruption differential matrix (`ValidateSamFile`/`ValidateVariants` vs `Validate`). |
| `run_p3_freeze_matrix.sh` | Phase 3 Step-50 freeze bundle (smoke profiles + required P3 conformance/runtime checks). |
| `normalize_assembly_region_igv.py` | Normalize Java `HaplotypeCaller --assembly-region-out` IGV lines for stable diffs. |
| `run_hc_full_parity_b1_read_shards.sh` | Rust read-shard TSV vs `parity/fixtures/hc-full-parity/b1/expected/` (uses `hc_full_parity_gate_dump` example; set `PARITY_CARGO_TARGET_DIR` / defaults to repo `target/`). |
| `run_hc_full_parity_b2_assembly_regions.sh` | Rust assembly-region iterator TSV vs `parity/fixtures/hc-full-parity/b2/expected/`. |
| `run_hc_full_parity_b3_apply_summary.sh` | `WalkerApplyStats` (apply count + inactive fast-path) vs `parity/fixtures/hc-full-parity/b3/expected/`. |
| `run_hc_full_parity_b4_walker_traversal.sh` | `traverse_assembly_region_walker` apply stats vs `parity/fixtures/hc-full-parity/b4/expected/` (same schema as B.3). |
| `run_hc_full_parity_phase_b.sh` | **Gate bundle B:** unit tests + B.1–B.4 L1 scripts (run before activity gates). |
| `run_hc_full_parity_c1_raw_activity.sh` | per-locus raw activity profile TSV vs `parity/fixtures/hc-full-parity/c1/expected/`. |
| `run_hc_full_parity_c2_smoothed_activity.sh` | band-pass smoothed activity vs `hc-full-parity/c2/expected/`. |
| `run_hc_full_parity_c3_active_locus.sh` | per-locus `is_active` vs `hc-full-parity/c3/expected/`. |
| `run_hc_full_parity_c4_gl.sh` | genotype-likelihood / MinimalGenotyping activity shortcut vs `hc-full-parity/c4-gl/expected/`. |
| `run_hc_full_parity_phase_c.sh` | **Bundle C:** prior walker gates + C.1–C.5 (incl. `c4-gl`, `c5-multi`, `c5-force`). |
| `run_hc_full_parity_c5_multi.sh` | **C.5:** multisample joint raw-activity TSV vs `c5-multi/` goldens. |
| `run_hc_full_parity_c5_force.sh` | **Tier B C.5.3:** `raw-activity-force` vs `c5-force/` goldens (forced alleles VCF). |
| `run_hc_full_parity_d1_read_filters.sh` | HC read-filter TSV vs `hc-full-parity/d1/expected/`. |
| `run_hc_full_parity_d2_downsample.sh` | positional + allele-biased downsample summaries. |
| `run_hc_full_parity_d3_soft_clip.sh` | HQ soft-clip mean per locus. |
| `run_hc_full_parity_phase_d.sh` | **Bundle D:** prior activity gates + D.1–D.3. |
| `run_p4_active_region_interval_diff.sh` | Phase 4 Step-58 harness: Java assembly-region IGV vs frozen expected corpus listed in `parity/fixtures/p4_assembly_region_cases.tsv` (requires indexed `sample.bam`). |
| `run_p4_freeze_matrix.sh` | Phase 4 Step-62 freeze bundle (smoke profiles + P4 contracts + assembly-region diff + bench smoke). |
| `run_p5_haplotype_candidate_diff.sh` | Phase 5 Step-70 scaffold: frozen Java-export candidate sets vs Rust local-assembly output. |
| `run_p5_runtime_candidate_diff.sh` | Phase 5 runtime candidate diff matrix with per-class match/drift report artifacts. |
| `run_p5_live_java_rust_diff.sh` | Phase 5 live runtime profile: Java HaplotypeCaller on synthetic regions with Rust candidate overlap checks against Java EventMap hap signatures. |
| `compare_haplotype_candidates.py` | Candidate-set comparator with exact-equality and cardinality drift metrics. |
| `run_p5_determinism_matrix.sh` | Phase 5 determinism matrix over thread counts and repeated runs (JSON summary). |
| `run_p5_mismatch_triage_check.sh` | Enforces triage schema categories and no unresolved high-severity mismatches. |
| `build_p5_equivalence_summary.py` | Builds markdown equivalence summary from phase-5 report artifacts. |
| `run_p5_assembly_stability_contract.sh` | Phase 5 Step-75 stability contract: repeated-run deterministic output checks (thread-count matrix). |
| `run_p5_freeze_matrix.sh` | Phase 5 Step-76 freeze bundle (assembly core tests + parity diff + regression + stability + bench smoke). |
| `run_p6_pairhmm_contracts.sh` | Phase 6 Wave-A contract gate (steps 77-79): PairHMM state machine, quality integration, and numeric-stability checks. |
| `run_p6_likelihood_vector_diff.sh` | Phase 6 Step-82 frozen likelihood-vector differential check (PairHMM output vs frozen Java baseline fixture). |
| `run_p6_wave_c_gates.sh` | Phase 6 Wave-C gate bundle (steps 83-86): boundary/artifact contracts + fp policy + PairHMM bench smoke. |
| `run_p6_determinism_matrix.sh` | Phase 6 determinism matrix over thread counts and repeated runs. |
| `run_p6_mismatch_triage_check.sh` | Phase 6 triage schema/disposition checks for mismatches. |
| `run_p6_live_java_refresh.sh` | Phase 6 auditable Java refresh: reruns live Java-vs-Rust runtime profile and fingerprints frozen Step-82 fixture. |
| `run_p6_live_pairhmm_drift.sh` | Phase 6 live PairHMM drift matrix: Java `Log10PairHMM` vs Rust PairHMM over multi-class corpus with JSON/MD delta report. |
| `run_p6_freeze_matrix.sh` | Phase 6 Step-88 freeze bundle for PairHMM hardening gates. |
| `run_p6_watertight_pass.sh` | Phase 6 watertight pass profile runner (`lite|live|full`). |
| `build_p6_equivalence_summary.py` | Builds markdown equivalence summary from Phase-6 report artifacts. |
| `run_p3_region_records_diff.sh` | Phase 3 region record-set differential (`PrintReads -L` vs `ListReadsInRegion`). |
| `diff_outputs.py` | Text diff: strict / normalized, optional regex + presence-only. |
| `compare_sam_parity.py` | Normalized SAM alignment-line parity (ignores `@PG`/`@CO`). |
| `compare_vcf_normalized.py` | **Scaffold** — sorted body comparison after stripping volatile `##` lines. |
| `compare_vcf_strict.py` | Strict VCF text after stripping a small volatile `##` header allowlist. |
| `compare_bam_alignment_parity.py` | BAM/SAM alignment parity via `samtools view` (sorted records + headers). |
| `report.py` | Markdown summary from JSON. |
| `run_p9_cli_contracts.sh` | Phase 9 Step 112: `gatk-cli` HaplotypeCaller integration tests + HC scaffold warmup + `gatk_cli_exit_code` unit tests. |
| `run_p9_hc_scaffold_diff.sh` | Phase 9 Step 113: Rust `HaplotypeCaller` scaffold VCF must match `parity/expected/p9_hc_scaffold_golden.vcf`. |
| `run_p9_java_hc_smoke.sh` | Phase 9 Step 113: Docker `gatk HaplotypeCaller` on parity fixtures (parsable VCF header; transient outputs removed). |
| `run_p9_mismatch_triage_check.sh` | Phase 9 triage schema + high-severity guard (`parity/reports/p9_mismatch_triage.jsonl`). |
| `run_p9_freeze_matrix.sh` | Phase 9 Step 114 freeze bundle (contracts + diff + Java smoke + triage + summary). |
| `build_p9_equivalence_summary.py` | Builds Phase-9 equivalence summary Markdown under gitignored `parity/reports/`. |
| `run_p10_release_readiness.sh` | Phase 10 Steps 115-120 release gate (coverage minimum + triage blocker checks + consolidated readiness summary). |
| `build_p10_release_readiness.py` | Builds `parity/reports/p10_release_readiness.{json,md}` from smoke + P7/P8/P9 summaries/triage artifacts. |
| `run_p11_hc_output_activation_contracts.sh` | Phase 11 scaffold: checks whether Rust HC output transitions from header-only scaffold to non-empty variant records on smoke fixture. |
| `run_p11_hc_output_field_diff_smoke.sh` | Phase 11 scaffold: readiness signal for Java-vs-Rust HC field-level differential promotion once Rust emits non-empty records. |
| `run_p11_hc_output_field_diff_corpus.sh` | Phase 11 Step-124 expanded strict corpus differential (Java-positive synthetic + real no-variant fixture). |
| `run_p11_mismatch_triage_check.sh` | Phase 11 Step-125 triage schema + high-severity guard (`parity/reports/p11_mismatch_triage.jsonl`). |
| `run_p11_freeze_matrix.sh` | Phase 11 scaffold freeze bundle (activation + field-diff readiness + summary). |
| `build_p11_equivalence_summary.py` | Builds `parity/reports/p11_equivalence_summary.{json,md}` with promotion criteria for strict output parity gates. |
| `run_p12_realworld_na12878_20k.sh` | Phase 12 real-world harness: downloads NA12878_20k_b37 public dataset and runs Java/Rust HC differential when `P12_REFERENCE` is set. Optional: `P12_CARGO_RELEASE=1` for faster Rust binary; `RAYON_NUM_THREADS` / `CARGO_BUILD_JOBS` for parallelism. |
| `run_p12_rust_only_na12878_20k.sh` | Re-run **only** Rust HC against the cached Java VCF (no Docker). Defaults to `P12_CARGO_RELEASE=1` and high `RAYON_NUM_THREADS`. Requires `P12_REFERENCE` and existing `P12_JAVA_VCF` (default: `parity/reports/p12_realworld_na12878_20k.java.vcf`). |
| `run_p12_l3_signoff.sh` | **L3 sign-off battery:** unit probes + L3a/L3b/ASM-8 gates + site trace (66/66, `rust_only=0`). Log: `parity/reports/p12_l3_signoff_<timestamp>.log`. |
| `run_p12_l4_signoff.sh` | **L4 sign-off battery:** fixture contracts + L3 regression + algorithmic + harness FORMAT (66/66). Log: `parity/reports/p12_l4_signoff_canonical.log`. |
| `deferred_features_audit.py` | Verify scope ADR 0001 + `gatk-tools` removed (ADR 0002) + deferred registry |
| `oracle_emit_audit.py` | Sprint L-3: `p12_java_only.tsv` / FORMAT TSV must not gate production emit |
| `excellence_gates_audit.py` | Sprint N: N-1…N-6 excellence gates (coords, env, size, bools, unwrap, claims) |
| `coord_allowlist.json` | Sprint N-1 frozen `923*****` allowlist |
| `run_p12_both_modes.sh` | Build release + run P12 with read augment on/off + cluster diagnose. |
| `diagnose_p12_cluster_assembly.sh` | Rust-only cluster assembly dumps (haplotypes, kmer-probe, assembly-stages). |
| `diagnose_p12_cluster_assembly_java_rust.sh` | Java vs Rust side-by-side on Java-active interval `2:92307228-92307400` (ASM-1). |
| `p12_na12878_summarize.py` | Shared JSON/MD summary for P12: `status` (harness) vs `parity_status` (variant set equality). |
| `run_p12_p13_realworld_full.sh` | Downloads hs37d5 + GIAB truth, picks/auto `P12_INTERVAL`, runs P12 then P13. |
| `p13_truth_eval.py` | P13 metrics: stratified SNP/INDEL F1, optional `--thresholds-json` + `--strict-gate` for L6. |
| `j6_truth_summarize.py` | L6 rollup: P12 parity + P13 gate → `parity/reports/hc-full-parity-j6/`. |
| `run_p13_realworld_truth_eval.sh` | Phase 13 truth harness: compares Java/Rust callsets against external truth VCF (`P13_TRUTH_VCF`) and reports precision/recall/F1. |
| `run_p14_multidataset_equivalence.sh` | Phase 14 multi-dataset suite: executes NA12878+GIAB path by default and conditionally runs Syndip + precisionFDA-ready cases when dataset env vars are provided; writes consolidated report. |
| `build_p14_equivalence_summary.py` | Builds `parity/reports/p14_multidataset_equivalence.{json,md}` from per-case P14 artifacts. |

`PARITY_REQUIRE_SAMTOOLS=1` fails the harness if `samtools` is missing (set in CI). Locally omit or set `0` if you do not have `samtools` installed.

## Acceptance modes (Phase 0 / Step 2)

- **Exit-only** — same exit class.
- **Normalized text** — `diff_outputs.py` whitespace-normalized or regex-sliced.
- **File SAM parity** — `compare_sam_parity.py`.
- **Bitwise VCF** — *planned*: strict comparator + files under `parity/expected/` (not yet required by CI).

Canonical acceptance class policy is frozen in `parity/ACCEPTANCE_CLASSES.md`.

`run_parity_smoke.sh` supports `PARITY_SMOKE_PROFILE=smoke|extended`:
- `smoke` (default) runs the fast baseline matrix,
- `extended` adds additional fixture checks for nightly depth.


`FOUNDATION_RUN_ADVISORY=1 ./scripts/parity/run_foundation_gate.sh` also executes non-gating advisory checks listed in `parity/checks.json`.

`run_p6_live_pairhmm_drift.sh` accepts optional tuning env vars for Rust-side matrix calibration:
- `P6_GAP_OPEN_PROB` (default `0.005`)
- `P6_GAP_EXTEND_PROB` (default `0.1`)
- `P6_INSERTION_EMISSION_PROB` (default `0.25`)

`run_p5_live_java_rust_diff.sh` accepts optional manifest override:
- `P5_LIVE_MANIFEST` (defaults to `parity/fixtures/p5_live_regions.tsv`; use `p5_live_regions_extended.tsv` for broader non-blocking coverage checks)

Live Java-vs-Rust overlap semantics:
- If Java emits `EventMap` hap signatures, at least one Rust candidate must overlap.
- If Java emits no `EventMap` signatures, the check is treated as non-applicable only when Java VCF has no variant records in that window.

Each check row supports:

- `acceptance_class` (for parity policy tracing),
- `timeout_s` (runner-enforced timeout), and
- `owner` (module/team accountability).

# Supported feature matrix — HaplotypeCaller

What `gatk-rs HaplotypeCaller` actually does with each user-facing flag.
Recorded from `main` at `756c0cf5bdc3240754ff806ad975a474d28b3f9f`
(PR #152, 2026-10-24). This task did not re-run parity gates and did not
add numerical results.

**This file is not the claim authority.** Agreement with GATK 4.4 stays in
[`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). How to run the command stays in
[`USER_GUIDE.md`](USER_GUIDE.md). A flag that Clap accepts is not, by that
fact, a supported control.

`GenotypeGVCFs --stand-call-conf` is a different command and is not
classified here.

## Classifications

| Label | Meaning in this file |
|-------|----------------------|
| Supported | The production command implements the behavior, and a test or a signed gate exercises it. |
| Partially supported | The behavior exists and is bounded. A missing CLI test, a scoped gate, or a stale comment is named in the row. |
| Deliberately rejected | The command fails before calling. The flag is not a live control. |
| Accepted but not effective | The command accepts the token and the call proceeds. The token does not change that behavior. |
| Unsupported | The feature is not a working HaplotypeCaller control. |
| Not validated | The behavior was not dynamically checked. Source inspection is labeled as such. |

Intent is either **Java-shaped** (the flag exists so a GATK 4.4 command can
be recognized, and the implemented behavior is aimed at that contract) or
**product-specific** (the behavior is a gatk-rs choice).

## Core inputs and outputs

| Flag / feature | Current behavior | Classification | Evidence | Limitations | Intent |
|----------------|------------------|----------------|----------|-------------|--------|
| `-R` / `--reference` | Required. A missing path, an unreadable FASTA, or a FASTA with zero contigs fails before calling. A sibling `.fai` supplies contig names and lengths; otherwise the FASTA is scanned. Contig lookup tolerates an optional `chr` prefix. | Supported | `gatk-haplotypecaller/src/run.rs` `run_haplotype_caller`; `gatk-core/src/reference.rs` `SequenceDictionary::from_fasta_path`, `contig_name_matches_fasta_entry`; CLI tests in `gatk-cli/tests/haplotype_caller_cli_integration.rs` (every HaplotypeCaller invocation passes `parity/fixtures/reference.fa`); [`USER_GUIDE.md`](USER_GUIDE.md); [`README.md`](../README.md) Inputs | No dedicated CLI test for a missing FASTA or a missing `.fai`. Those failures were read from source. The loader does not read a Picard `.dict` (README). | Java-shaped |
| `-I` / `--input` indexed BAM | Required. Each path must be an existing file. Region traversal opens it with `IndexedReader::from_path` and `fetch`. The VCF sample name is the first BAM’s `@RG SM`, otherwise `SAMPLE`. | Supported | `run.rs` input loop and `sample_name_from_bam_inputs`; `gatk-haplotypecaller/src/assembly_region_iterator.rs` `load_records_for_shard_raw`; CLI tests use `parity/fixtures/sample.bam`; [`USER_GUIDE.md`](USER_GUIDE.md) | A missing `.bai` is not a dedicated CLI test. Open/fetch failure is the indexed-reader error path, inspected in source. | Java-shaped |
| Repeated `-I` | Each BAM is traversed on its own. Variant rows are merged into one VCF. More than one `-I` logs that this is not Java joint calling. | Partially supported | `run.rs` `validate_haplotype_caller_input_policy`; [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) No row “Multi-sample joint HC”; [`USER_GUIDE.md`](USER_GUIDE.md); [`README.md`](../README.md) Inputs | No CLI test passes two BAMs. The warning and the per-BAM loop were inspected in source. Joint genotyping of reads from several BAMs is not the behavior. | Product-specific. Java merges the read sets. |
| CRAM | Not a HaplotypeCaller input in the user guide, README, or CLI tests. The caller does not set a CRAM reference on the indexed reader. | Unsupported | [`USER_GUIDE.md`](USER_GUIDE.md) and [`README.md`](../README.md) specify indexed BAM. `gatk-core` CRAM copy helpers (`gatk-core/src/io/hts_copy.rs`, `gatk-core/tests/p3_io_conformance_tests.rs`) are not the HaplotypeCaller command. | CRAM-through-HC was not executed. Those core helpers are not evidence that `HaplotypeCaller -I file.cram` works. | Not offered |
| `-O` / `--output` plain VCF | Writes uncompressed VCF text, including names such as `.g.vcf` and names with no `.vcf` suffix. | Supported | `gatk-core/src/io/vcf.rs` `VcfWriter::new`; CLI test `haplotype_caller_plain_vcf_is_text_and_compressed_suffixes_are_rejected`; [`USER_GUIDE.md`](USER_GUIDE.md) Outputs | The CLI test checks a scaffold header (`GATK_RS_HC_SCAFFOLD_OUTPUT=1`) plus the text magic bytes, not a variant-body oracle. | Java-shaped for the VCF text. Compression is a separate row. |
| `.vcf.gz` / `.bgz` / `.bgzip` / `.bcf` | Rejected before the path is created or truncated. An existing file with that name is left unchanged. Exit code 2. | Deliberately rejected | `gatk-core/src/io/vcf.rs` `reject_unsupported_plain_vcf_output`; same CLI test; [`USER_GUIDE.md`](USER_GUIDE.md) | Suffix match is on the file name only, ASCII case-insensitive. | Product-specific rejection. Java writes BGZF and BCF. |
| Output index (`.tbi` / `.csi`) | No index is written. `OutputConfig.create_output_variant_index` defaults to `true` in `gatk-common/src/config.rs` and has no reader anywhere else in the tree. | Unsupported | `VcfWriter::new` writes a plain `File`. The field is assigned only in `GatkConfig::new`. [`USER_GUIDE.md`](USER_GUIDE.md) Outputs states that an index is not implemented. | The default `true` was inspected statically. It does not create an index. | Product-specific. Indexing is not implemented. |
| `-L` / `--intervals`, repeated | Repeatable. Tokens are joined with `;` and checked against the FASTA dictionary before the run. An unknown contig exits 2. | Supported | `gatk-cli/src/main.rs` `join_interval_cli_args`; `gatk-core/src/reference.rs` `parse_intervals_cli_string`; CLI tests `haplotype_caller_accepts_repeated_interval_flags`, `haplotype_caller_invalid_interval_user_error_exit_2` | The repeated-interval CLI test uses scaffold output and checks that the process exits 0. It does not compare called sites across the two intervals. | Java-shaped |
| `-L` omitted | The interval list is every contig in the dictionary, full length. | Partially supported | `gatk-core/src/reference.rs` `intervals_for_haplotype_caller` (`None` → `whole_genome_interval_specs`); [`USER_GUIDE.md`](USER_GUIDE.md); [`README.md`](../README.md) Inputs | No CLI test omits `-L`. This row is source inspection plus the user-guide sentence. Omitting `-L` on a whole genome is a large run, not a validated genome-wide result. | Java-shaped |

## Output modes

### Unknown `--output-mode` on current `main`

The production command does **not** treat an unknown `--output-mode` as
variant-only VCF.

`gatk-cli/src/main.rs` stores the token with `set_output_mode` and then calls
`GatkConfig::validate` before `run_haplotype_caller`. `validate` in
`gatk-common/src/config.rs` returns `GatkError::configuration` unless the
mode is exactly `VCF`, `GVCF`, or `BP_RESOLUTION`. That error class exits 2
(`gatk-common/src/error.rs` `gatk_cli_exit_code`). `--mode` is a Clap alias
of the same field.

`emit_mode_from_output_mode` in
`gatk-haplotypecaller/src/reference_vcf_emit.rs` still maps every other
string to `EmitMode::Vcf`. `run_haplotype_caller` calls that function only
after CLI validation has passed. A caller that skips `GatkConfig::validate`
still gets the variant-only fallback. That fallback is not a supported CLI
behavior.

The unit test `emit_mode_parsing_matches_gatk_output_modes` asserts `VCF`,
`GVCF`, and `BP_RESOLUTION` only. There is no CLI integration test that an
unknown `--output-mode` exits non-zero and leaves the output path uncreated.
This section is source inspection of `756c0cf5`. It is not a dynamic run
performed for this document.

The Step 4A inventory described a silent CLI fallback. That description
matches the mapper. It does not match the command, which validates first.

| Flag / feature | Current behavior | Classification | Evidence | Limitations | Intent |
|----------------|------------------|----------------|----------|-------------|--------|
| `--output-mode VCF` (default) and alias `--mode VCF` | Variant-only emission. The default pipeline header is `assembly-region-v1`. | Supported | Clap default in `gatk-cli/src/main.rs`; `emit_mode_from_output_mode`; CLI test `haplotype_caller_default_pipeline_is_assembly_region_v1`; [`USER_GUIDE.md`](USER_GUIDE.md) | The CLI test checks the pipeline header on `parity/fixtures/sample.bam` / `chr1:1-32`. It does not prove that a rich interval emits only variant records. | Java-shaped (`EMIT_VARIANTS_ONLY` analogue) |
| `--output-mode GVCF` | Selects `EmitMode::Gvcf`, adds gVCF header lines, and emits reference-confidence blocks on the gVCF path. | Partially supported | `run.rs` (`emit_mode == EmitMode::Gvcf`); `reference_vcf_emit.rs` `inactive_reference_model_to_vcf_records`; `scripts/parity/run_p12_l5_gvcf.sh` invokes `gatk-rs HaplotypeCaller --output-mode GVCF`; [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) P12 L5 Yes row and waiver **W-H3** | The in-tree CLI integration file never asserts a gVCF body. `haplotype_caller_emit_ref_confidence_long_form_accepted` runs scaffold output with `--emit-ref-confidence NONE`. P12 L5 is one interval, `chr2:92300000-92350000`. W-H3 keeps RCM band reconciliation on that interval. This document did not re-run L5. | Java-shaped inside that signed interval |
| `--output-mode BP_RESOLUTION` | Selects `EmitMode::BpResolution`. Locus decisions for reference sites exist in the genotyping module. | Partially supported | `emit_mode_from_output_mode`; `gatk-haplotypecaller/src/genotyping.rs` `decide_locus_emission` tests for `EmitMode::BpResolution`; [`USER_GUIDE.md`](USER_GUIDE.md) Outputs | No CLI integration test writes a BP_RESOLUTION body. [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) has no BP_RESOLUTION Yes row. Java compatibility of the body is not validated. | Java-shaped token. Body parity is not a claim. |
| Unknown `--output-mode` / `--mode` | Configuration error before calling. See the section above. | Deliberately rejected | `GatkConfig::validate`; `gatk_cli_exit_code` | Rejection is source-inspected. No CLI test covers an unknown token. The mapper fallback remains for callers that skip validation. Do not describe that fallback as supported. | Product-specific fail-closed allow-list |
| `--emit-ref-confidence GVCF` or `BP_RESOLUTION` | Replaces `output_mode` with that token. There is no `-ERC` short flag. | Partially supported | `gatk-cli/src/main.rs` match on `emit_ref_conf`; same gVCF evidence as `--output-mode GVCF` for the GVCF token; [`USER_GUIDE.md`](USER_GUIDE.md) | The only CLI integration test passes `NONE` under scaffold output. It does not check a reference-confidence body. `BP_RESOLUTION` has the same gap as the output-mode row. | Java-shaped long form of `-ERC` |
| `--emit-ref-confidence NONE`, or any other token | The token is stored. `output_mode` is left as already set (default `VCF`, or an explicit `--output-mode`). | Accepted but not effective | `gatk-cli/src/main.rs` match arm `_ => {}`; CLI test `haplotype_caller_emit_ref_confidence_long_form_accepted` (`NONE`, scaffold header only) | Unknown tokens are not rejected. A typo does not select gVCF and does not fail the run. `NONE` does not force variant-only output when `--output-mode` is already `GVCF`. | Product-specific. Java `-ERC NONE` is the variant-only mode; this command does not use the token that way. |

## Calling thresholds and controls

Explicit defaults are accepted. They are not knobs. Non-default values of
the four controls below fail closed (PR #152). Wiring them into activity or
assembly was intentionally not done.

| Flag / feature | Current behavior | Classification | Evidence | Limitations | Intent |
|----------------|------------------|----------------|----------|-------------|--------|
| `--stand-call-conf` / `--standard-min-confidence-threshold-for-calling` omitted, or explicit `30` / `30.0` | Accepted. The production activity scorer does not read `HaplotypeCallerConfig.stand_call_confidence`. Active-region confidence stays at `ACTIVE_REGION_STANDARD_CONFIDENCE_CAP` (`4.0`). `from_haplotype_caller_config` would copy the CLI value and cap it at 4.0, and nothing in the crate calls it. | Accepted but not effective | `gatk-cli/src/main.rs` `Option<f64>` plus `unwrap_or(30.0)`; `run.rs` `reject_unsupported_haplotype_caller_cli` constant `HC_DEFAULT_STAND_CALL_CONFIDENCE`; `activity_scoring.rs` `ACTIVE_REGION_STANDARD_CONFIDENCE_CAP` and `from_haplotype_caller_config`; CLI test `haplotype_caller_explicit_production_defaults_still_call`; [`USER_GUIDE.md`](USER_GUIDE.md) | The acceptance test sets `GATK_RS_HC_SCAFFOLD_OUTPUT=1` and checks the scaffold header. It does not show that calling used 30.0. Equality with `30.0` is exact. | Java-shaped default value. The flag is not a live control. |
| `--stand-call-conf` any other explicit value | Configuration error. The VCF is not created. Exit code 2. | Deliberately rejected | `run.rs` `reject_unsupported_haplotype_caller_cli`; unit test `reject_unimplemented_call_confidence_and_base_quality_controls`; CLI test `haplotype_caller_rejects_non_default_stand_call_conf_without_creating_vcf` | Rejection is the product behavior. It is not algorithmic configurability. | Product-specific fail-closed |
| `--stand-emit-conf` / `--standard-min-confidence-threshold-for-emitting` | Parsed (default `10.0`) and passed into production emit as `stand_emit_confidence`, including `passes_hc_variant_emit_biallelic`. | Partially supported | `gatk-cli/src/main.rs`; `run.rs` `assembly_region_variant_records`; `gatk-haplotypecaller/src/hc_genotyping_engine/mod.rs` and `genotype_finalize.rs`; [`USER_GUIDE.md`](USER_GUIDE.md) | No CLI test shows that a non-default value adds or drops a record. Engine tests that use the default do not by themselves prove the CLI wiring for another number. This gap was source-inspected. | Java-shaped |
| `--min-base-quality-score` omitted, or explicit `10` | Accepted. Activity scoring’s default `min_base_quality_score` is 10. `ReadThreadingAssemblerArgs::default` sets `min_base_quality: 10`. Neither copies the CLI field. `from_haplotype_caller_config` is unused. | Accepted but not effective | `gatk-cli/src/main.rs`; `run.rs` `HC_DEFAULT_MIN_BASE_QUALITY_SCORE`; `activity_scoring.rs` `HaplotypeCallerActivityScoringParams::default`; `read_threading_assembler.rs` `ReadThreadingAssemblerArgs::default`; CLI test `haplotype_caller_explicit_production_defaults_still_call` | The acceptance test is scaffold output, as with `--stand-call-conf 30`. PairHMM’s separate base-quality cap is `HC_DEFAULT_BASE_QUALITY_SCORE_THRESHOLD` (18) in `likelihood_engine.rs`. That constant is not this flag. | Java-shaped default value. The flag is not a live control. |
| `--min-base-quality-score` any other explicit value | Configuration error. The VCF is not created. Exit code 2. | Deliberately rejected | `run.rs` `reject_unsupported_haplotype_caller_cli`; unit test `reject_unimplemented_call_confidence_and_base_quality_controls`; CLI test `haplotype_caller_rejects_non_default_min_base_quality_without_creating_vcf` | Same as non-default `--stand-call-conf`. | Product-specific fail-closed |
| `--min-mapping-quality` | Applied. `ReadFilterParams::from_haplotype_caller` copies it. Values above 255 clamp to 255. Default 20. Supplementary reads are kept (`exclude_supplementary: false`), matching the Java standard HC filter comment in that function. | Partially supported | `run.rs` (filters built from `hc`); `read_model.rs` `from_haplotype_caller`, `mapq_passes_minimum`; [`USER_GUIDE.md`](USER_GUIDE.md) | No CLI test shows a higher floor dropping reads. That check is a source reading of `from_haplotype_caller`. Contract tests that fix the floor at 20 do not exercise another CLI value. | Java-shaped |
| `--pair-hmm`, `--pair-hmm-implementation`, `--pairhmm-impl` | Omitted selects `PairHmmImpl::FastestAvailable`. Known tokens parse. An unknown token fails before calling. `score_read_against_haplotypes` uses `resolved_pair_hmm_backend`. `FastestAvailable` resolves to AVX2, NEON, or packed f64 when `best_simd_backend` returns one of those, otherwise logless scalar. `WAVEFRONT` resolves only when requested. | Partially supported | `run.rs` `pair_hmm_impl_from_config`, tests `omitted_pair_hmm_selects_fastest_available`, `accept_pair_hmm_cli_flag`, `reject_unknown_pair_hmm_cli_flag`; `pairhmm_simd/dispatch.rs` `parse_pair_hmm_impl`, `resolve_pair_hmm_impl`; `likelihood_engine.rs` `score_read_against_haplotypes`; CLI test `haplotype_caller_help_pair_hmm_omitted_default_is_fastest_available`; [`USER_GUIDE.md`](USER_GUIDE.md); [`README.md`](../README.md) | The CLI test checks help text, not VCF bytes across kernels. The unknown-token test calls `reject_unsupported_haplotype_caller_cli` directly, not the `gatk-rs` process. Signed parity gates are not a claim that every backend matches Java’s log10 PairHMM. Comments in `resolve_pair_hmm_impl` and on `HcLikelihoodImplementation` still say production scoring stays on log10. Those comments disagree with `pair_hmm_impl_from_config` and with `resolve_pair_hmm_impl`. | Product-specific default (`FASTEST_AVAILABLE`). `LOG10_PAIRHMM` is the explicit Java-shaped scalar kernel. |
| `--max-alternate-alleles` | Hidden. Any value other than 6 is a configuration error and does not create the VCF. The production caller never reads the field except in that check. | Deliberately rejected | `gatk-cli/src/main.rs` (`hide = true`, default 6); `run.rs` `HC_DEFAULT_MAX_ALTERNATE_ALLELES`; unit test `reject_non_default_max_alternate_alleles`; CLI test `haplotype_caller_rejects_non_default_max_alternate_alleles_without_creating_vcf` | Explicit `6` is accepted because it equals the constant. It is not a live cap wired into allele selection. The flag is absent from the user-guide table and hidden in help. | Product-specific fail-closed. Java default 6 is the constant, not a CLI control. |
| `--original-base-qualities` enabled | Configuration error. The VCF is not created. Exit code 2. Omitted stays off. | Deliberately rejected | `gatk-cli/src/main.rs`; `run.rs` `reject_unsupported_haplotype_caller_cli`; unit test `reject_unimplemented_call_confidence_and_base_quality_controls`; CLI test `haplotype_caller_rejects_original_base_qualities_without_creating_vcf`; [`USER_GUIDE.md`](USER_GUIDE.md) | The flag is not applied to PairHMM or activity. | Product-specific fail-closed |
| `--dont-use-soft-clipped-bases` enabled | Configuration error. The VCF is not created. Exit code 2. Omitted stays off. | Deliberately rejected | Same reject function and unit test; CLI test `haplotype_caller_rejects_dont_use_soft_clipped_bases_without_creating_vcf`; [`USER_GUIDE.md`](USER_GUIDE.md) | `read_unclip.rs` is a dump-path helper. It is not this production call. | Product-specific fail-closed |

## Assembly and alleles

| Flag / feature | Current behavior | Classification | Evidence | Limitations | Intent |
|----------------|------------------|----------------|----------|-------------|--------|
| `--alleles` / `-A` | The VCF path is stored and loaded into `CallRegionArgs.given_alleles` and `assemble.given_alleles` on the production run. | Partially supported | `gatk-cli/src/main.rs`; `run.rs` `assembly_region_variant_records`; `given_alleles.rs` `load_given_alleles_from_vcf_path`; test `given_alleles_extend_trim_variants_without_duplicate` in `gatk-haplotypecaller/tests/call_region_given_alleles_order_test.rs`; [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) T3-5; [`USER_GUIDE.md`](USER_GUIDE.md) | No CLI integration test passes `-A`. T3-5 says broader sign-off is optional beyond L2 `c5-force`. The module comment at the top of `given_alleles.rs` still says the CLI flag is not exposed. That comment is false on this `main`. | Java-shaped. Not a full `-alleles` sign-off. |
| Assembly-region calling | The default command runs `assembly-region-v1` on `CallRegionArgs::default()`, which is `strict_java()`. | Supported | `engine.rs` `Default for CallRegionArgs` and `strict_java`; `run.rs` `assembly_region_variant_records`; CLI test `haplotype_caller_default_pipeline_is_assembly_region_v1`; [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) production-path Yes row | The CLI test checks the header token on the tiny fixture. Algorithm agreement is the scoped gates in the parity section, not that header. | Java-shaped on the signed gates |
| `--assembly-region-out` | Not an argument of `HaplotypeCaller`. Help text points at the claim matrix. `DumpSmoothedActivity` is a separate command. | Unsupported | `gatk-cli/src/main.rs` HaplotypeCaller help string and the `DumpSmoothedActivity` command; [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) T5-6 | T5-6 says to use `DumpSmoothedActivity` for analysis. That command is not assembly-region output from HaplotypeCaller. | Not offered on this command |

## Execution

| Flag / feature | Current behavior | Classification | Evidence | Limitations | Intent |
|----------------|------------------|----------------|----------|-------------|--------|
| `--threads` / `-t` / `--nt` | Sets the Rayon pool before calling. `0` uses the logical CPU count. An explicit value overrides `RAYON_NUM_THREADS`. Omitted honors `RAYON_NUM_THREADS` when that value is a positive integer, otherwise the logical CPU count. | Supported | `gatk-cli/src/main.rs` `init_haplotype_caller_threads`; `gatk-cli/tests/rayon_thread_determinism_test.rs` `haplotype_caller_vcf_byte_identical_across_threads_cli_flag` and `haplotype_caller_vcf_byte_identical_across_rayon_thread_counts`; [`USER_GUIDE.md`](USER_GUIDE.md) | The byte-identity test invokes `--threads` at 1, 2, 4, 8, and 16. `-t` and `--nt` are the same Clap argument; that identity was read from the attribute, not from a separate process test. Comparison drops `##fileDate=`, `##source=`, and `##GATKCommandLine`. The fixture is the one that test builds (`parity/fixtures`). | Product-specific pool. The names are Java-shaped. |
| Deterministic calling | On that same fixture, normalized VCF bytes match across those thread counts. Region batches are sorted and deduplicated before write. | Partially supported | `rayon_thread_determinism_test.rs`; [`USER_GUIDE.md`](USER_GUIDE.md) Outputs; `run.rs` `merge_region_emit_batches_is_position_ordered_and_deduped` | The user guide limits the promise to the fixture that test uses. This row does not extend it. | Product-specific |
| `--java-options` | Accepted on the top-level command for launcher familiarity. A stderr note says no JVM is started. `GatkConfig::with_java_options` parses the string. If `-Xmx` is present, validation rejects a memory token that does not start with a digit and end in `M` or `G`. Other tokens are stored and not applied to the Rust process. | Accepted but not effective | `gatk-cli/src/main.rs` (`java_options` field, stderr note, `with_java_options`); `gatk-common/src/config.rs` `JavaOptions::parse` and `validate`; [`USER_GUIDE.md`](USER_GUIDE.md) | No CLI test of a bad `-Xmx` was found in `haplotype_caller_cli_integration.rs`. The parse-and-ignore path was read from source. There is no JVM. | Product-specific compatibility shim |

## Parity scope

These are the validated calling claims already recorded on `main`. This
matrix does not add Yes rows. Logs under `parity/reports/` are gitignored
(`.gitignore`), so the raw L2 log was not re-opened for this document. The
numbers below are the ones [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) and
[`PARITY.md`](PARITY.md) already publish.

| Claim | Bound already on `main` | Evidence | What this is not |
|-------|-------------------------|----------|------------------|
| Production path | `CallRegionArgs::strict_java()` and `assembly-region-v1` | [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md); `engine.rs`; `run.rs` | A statement that every CLI flag is live |
| Canonical mid-B | One ActiveFull region, `2:92317000-92319000`, three oracle SNPs | [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md); [`PARITY.md`](PARITY.md); [`PARITY_MILESTONE_6R.md`](PARITY_MILESTONE_6R.md); `p12_call_none_mid_b_test` | Whole-codebase HC equivalence. [`PARITY.md`](PARITY.md) also says that test returns without asserting unless the gitignored reference and BAM are present. |
| P12 L3 / L4 / L5 | `chr2:92300000-92350000` | [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md); `scripts/parity/run_p12_l3_signoff.sh`, `run_p12_l4_signoff.sh`, `run_p12_l5_gvcf.sh` | A genome-wide gVCF or FORMAT claim. L5 RCM stays under waiver **W-H3**. |
| Strict L2 | **223 passed, 0 failed**, strict, dated **2026-07-24** in the claim matrix | [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md); [`PARITY.md`](PARITY.md) L2 row; `scripts/parity/run_hc_full_parity_l2.sh`; permanent May gate `cargo test -p gatk-haplotypecaller --features parity_harness --test l2_may_regression_gate_test` | A new L2 run. This document did not execute the script. |
| PairHMM | Omitted CLI kernel is `FASTEST_AVAILABLE` (see the `--pair-hmm` row). Signed gates are the rows above. | [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md); [`PARITY.md`](PARITY.md) “Fresh production VCF comparisons” states that the 2026-10-23 local VCFs are not a Yes row and not parity for every PairHMM backend | Proof that AVX2, NEON, packed f64, logless, and log10 emit the same VCF. `WAVEFRONT` is opt-in and outside `FastestAvailable`. |
| HG001 GIAB ci-subset | The claim matrix’s signed ΔF1 row, chr20/21 windowed ci-subset | [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) Asserted and Pending sections | Full autosomes, HG002+HG005, or a clinical callset |

Genome-wide equivalence is not a claim. Full Java GATK feature parity is not
a claim. Both are already No rows in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md).
Historical L6–L14 narratives on `pre-cleanup-archive` are not Yes rows on
`main`.

The 2026-10-23 fresh VCF notes in [`PARITY.md`](PARITY.md) stay what that
page says they are: local comparisons on a detached checkout, not Yes rows.

## Explicitly out of scope

These exclusions are the ones already written in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md)
and the HaplotypeCaller help text. This list does not survey every Java GATK
argument that happens to be absent.

| Feature | Classification | Evidence |
|---------|----------------|----------|
| bamout / HaplotypeBAMWriter | Unsupported | Claim matrix T5-1. Named in `gatk-cli/src/main.rs` long help as deferred. Not a field of the `HaplotypeCaller` command. |
| DRAGSTR | Unsupported | Claim matrix T5-2. `HcLikelihoodEngineConfig` defaults `dragstr_params_loaded` to false (`likelihood_engine.rs`). No HaplotypeCaller CLI flag. |
| DRAGEN | Unsupported | Claim matrix T5-3. Help text defers it. |
| `AS_*` annotations | Unsupported | Claim matrix T5-4. |
| DP reconciliation | Unsupported | Claim matrix T5-5. |
| `--assembly-region-out` | Unsupported | Claim matrix T5-6. See the assembly section. |
| BQSR (`BaseRecalibrator`) | Unsupported | Claim matrix No row; [`docs/adr/0001-scope-boundary.md`](adr/0001-scope-boundary.md); [`USER_GUIDE.md`](USER_GUIDE.md) states this project does not recalibrate base qualities. |
| VQSR (`VariantRecalibrator` / `ApplyVQSR`) | Unsupported | Claim matrix No row. `VariantFiltration` hard filters are a different command. |
| Other Java HaplotypeCaller options with no row above | Unsupported | They are not arguments of the `HaplotypeCaller` command in `gatk-cli/src/main.rs`. Absence from this matrix is not an implied implementation. |

## Validation gaps

These are open on purpose. This document records them. It does not close them.

| Gap | Disposition |
|-----|-------------|
| No CLI test that an unknown `--output-mode` exits 2 and does not create a file | Later validation. The command path is fail-closed in `GatkConfig::validate`. The residual `emit_mode_from_output_mode` fallback is a documented limitation for callers that skip validation. |
| No in-tree CLI test that `--output-mode GVCF` or `--emit-ref-confidence GVCF` writes reference-confidence blocks | Later validation outside the signed P12 L5 interval. P12 L5 already invokes the production CLI with `--output-mode GVCF` on one interval. |
| No CLI test that a non-default `--stand-emit-conf` changes records | Later validation. The value is passed into emit. |
| No CLI test that a higher `--min-mapping-quality` drops reads | Later validation. The value is copied onto `ReadFilterParams`. |
| No CLI test of two `-I` BAMs | Documented limitation. Repeated `-I` is per-BAM traversal, not Java joint calling. |
| CRAM as `HaplotypeCaller -I` | Unsupported. Not a product input. |
| BGZF, BCF, and an output index | Deliberately rejected (compressed suffixes, tested) or unsupported (index, no writer). |
| Genome-wide or full-Java-feature parity | Not a claim. Bounded by [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). |
| Unknown `--emit-ref-confidence` tokens | Documented limitation: accepted and not effective. |
| Stale comments in `given_alleles.rs`, `resolve_pair_hmm_impl`, and `HcLikelihoodImplementation` | Documented. User-facing help and [`USER_GUIDE.md`](USER_GUIDE.md) match the CLI behavior for PairHMM and `-A`. The comments do not. |

## Reading order

1. This file for what a HaplotypeCaller flag does.
2. [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) before citing GATK agreement.
3. [`USER_GUIDE.md`](USER_GUIDE.md) for the command line a user should type.

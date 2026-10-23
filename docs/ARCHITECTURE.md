# Architecture

Native Rust workspace for a GATK 4.4–aligned germline short-variant path.
Product claims: [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). How those claims are
tested: [`PARITY.md`](PARITY.md). Independence and trademarks:
[`../NOTICE.md`](../NOTICE.md).

The implementation is Rust modules and data types aimed at observable Java
behavior. It is not a class-by-class port of GATK.

## Crates

| Crate | Role |
|-------|------|
| [`gatk-cli`](../gatk-cli/) | `gatk-rs` binary. Clap tools, thread-pool setup, process exit codes. |
| [`gatk-haplotypecaller`](../gatk-haplotypecaller/) | Calling engine: activity through VCF/gVCF, plus CombineGVCFs and GenotypeGVCFs. |
| [`gatk-core`](../gatk-core/) | FASTA, BAM, and VCF I/O via HTSlib; intervals; sequence dictionary; VariantFiltration. |
| [`gatk-common`](../gatk-common/) | Shared errors and configuration types. |
| [`gatk-rs-equiv`](../gatk-rs-equiv/) | Separate CLI that scores Rust vs Java callsets with hap.py or RTG vcfeval. |
| [`fuzz/`](../fuzz/) | LibFuzzer target excluded from the workspace members list. Shares scenarios with `gatk-rs-equiv`. |

`gatk-cli` enables the haplotypecaller `dev-dumps` feature so
`DumpSmoothedActivity` can build. That does not turn dump modules into the
default HaplotypeCaller path. Feature `parity_harness` is for tests and L2
scripts. It is not the release CLI surface.

There is no `gatk-tools` crate
([ADR 0002](adr/0002-remove-gatk-tools.md)).

## Call path

Default `HaplotypeCaller` uses `CallRegionArgs::strict_java()` and the
pipeline id `assembly-region-v1` (`HC_PIPELINE_ASSEMBLY_REGION_V1` in
`region_vcf_emit`). `GATK_RS_HC_SCAFFOLD_OUTPUT=1` switches that command to
a header-only scaffold used by tests. It is not the production path.

Within one assembly region the engine stages are:

1. **Activity and regions** — `activity_profile`, `active_region`,
   `assembly_region_iterator`. Reads are fetched from an indexed BAM over
   padded spans.
2. **Local assembly** — `read_threading_graph`, `read_threading_assembler`,
   `seq_graph`, dangling-tail recovery, pruning. Haplotypes come from a
   k-best search (`kbest_haplotype`, `seq_kbest_haplotype`).
3. **Events** — `event_map` builds variation events from haplotypes.
   `assembly_region_trimmer` trims haplotypes before the likelihood engine.
4. **Likelihoods** — `pairhmm`. When `--pair-hmm` is omitted,
   `pair_hmm_impl_from_config` selects `FASTEST_AVAILABLE`, which resolves
   to a supported host SIMD kernel (AVX2, NEON, or packed f64) or else the
   Logless scalar kernel. `pairhmm_log10` runs when the flag is
   `LOG10_PAIRHMM`. `pairhmm_logless` and `pairhmm_simd` are the other
   kernels. Smith-Waterman realignment is `smith_waterman`.
5. **Genotyping** — `hc_genotyping_engine` / `genotyping`: allele
   likelihoods, genotype likelihoods, QUAL, and the decision to emit.
6. **Annotations and records** — `variant_site_hc_annotations`,
   `region_vcf_emit`, and `ref_confidence` for gVCF blocks.

`combine_gvcfs` and `genotype_gvcfs` are separate from per-sample
HaplotypeCaller. Hard filters are `gatk_core::variant_filtration`, invoked
by the `VariantFiltration` subcommand, not inside the caller.

Region work can run in parallel (`run.rs`, Rayon). Results are merged by
genomic position and `region_index`, then deduplicated, before the VCF write.
See [Concurrency](#concurrency-and-determinism).

## K-best resource policy

Completed-path **K** (`num_best_haplotypes_per_graph`, 128 on
`CallRegionArgs::strict_java`) is not the frontier limit.

When `GATK_RS_EXPERIMENTAL_KBEST_POLICY` is unset, the production policy is
`unbounded_diagnostic`: the search does not refuse at a 1024-long frontier.
`legacy_1024` is selected only by setting that variable. `byte_budget` is an
experimental ceiling and needs `GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET`.
Neither env var is a public CLI flag or a claim-matrix Yes row.
Detail: [`parity/KBEST_RESOURCE_POLICY.md`](parity/KBEST_RESOURCE_POLICY.md).

Older 6R prose that still calls `legacy_1024` the production default is
historical. The investigation log is the second half of [`PARITY.md`](PARITY.md).

## Native dependencies and external tools

- **HTSlib** through `rust-htslib` / `hts-sys` for BAM, VCF, and FASTA.
  The crate builds HTSlib; the project does not shell out to GATK.
- **samtools** for indexes and for parity comparators. Not linked into
  `gatk-rs`.
- **Java GATK 4.4.0.0** only in the parity and equivalence harnesses
  ([`GATK_PINNED.env`](GATK_PINNED.env)).
- **hap.py or RTG** only for `gatk-rs-equiv` and GIAB scoring.

## Concurrency and determinism

`--threads` sizes the global Rayon pool (`worker_threads`). Active regions
inside one process run in parallel. The merge in `run.rs` orders batches
before writing.

The CI test `gatk-cli/tests/rayon_thread_determinism_test.rs` requires a
HaplotypeCaller VCF to be byte-identical for `RAYON_NUM_THREADS` and
`--threads` in `{1, 2, 4, 8, 16}` after volatile header lines are
normalized. That is the determinism contract enforced in
`.github/workflows/ci.yml`, `quality.yml`, and `pr-check.yml`. It applies
to the fixture that test runs. It does not certify every external BAM.

`GATK_RS_HC_SEQUENTIAL=1` forces one region at a time. It is a memory
control, not a different calling algorithm.

QualByDepth’s Java `Random` stream is part of the mid-B numeric contract
(seed `47382911`). A different set of annotated sites can move later QD
values even when the QD formula matches.

## Extension points

| Change | Where it belongs |
|--------|------------------|
| New INFO/FORMAT annotation on the HC path | `variant_site_hc_annotations` and the emit path that already consumes those values. Add a regression that states the Java contract. |
| Assembly, k-best, PairHMM, or genotyping behavior | The stage module above. Do not add a parallel Java-shaped type tree. |
| Hard-filter expression | `gatk-core` `variant_filtration`. |
| A new GATK tool outside HC → Combine → Genotype → hard filter | Out of scope unless [ADR 0001](adr/0001-scope-boundary.md) is revised. Do not stub it. |
| A wider product claim | [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) in the same change, with evidence already on `main`. |

P12-only behavior is gated by interval overlap and the waivers in the claim
matrix (W-H1, W-H3, and the W-J tables). Do not copy a P12 branch onto every
locus to turn a fixture green.

## Design constraints and debt

- Algorithm parity with pinned GATK 4.4 on the signed gates. Java source
  layout is not a goal.
- `parity_harness` and the `*_dump` modules exist so L2 and forensic tests
  can see intermediate state. They are not a second production caller.
  With that feature, `P12_BASELINE_EMIT_FILTER` restricts strict emit to
  Java baseline VCF rows for comparison. Production `strict_java` does not
  enable it. These harness env flags are read only under `parity_harness`
  or tests:

  | Flag | Effect when set |
  |------|-----------------|
  | `P12_PHASE_E` | L3 gate harness. Also enables the baseline emit filter. |
  | `P12_BASELINE_EMIT_FILTER` | Restrict strict emit to Java baseline VCF rows. |
  | `GATK_RS_P12_EVENT_REGISTRY` | Inject the Phase-E gap-SNP registry into the graph-only EventMap. |
  | `GATK_RS_P12_ENSURE_BRIDGES` | Enable P12 list-inject and post-HMM `ensure_*` bridges. Forced off when `GATK_RS_ASM8_ONLY` is set. |
  | `P12_L4_JAVA_FORMAT` | Overlay FORMAT, QUAL, and AF from `parity/fixtures/p12-java-format/all_sites.tsv`. Not implied by `P12_PHASE_E`. |
  | `GATK_RS_ENABLE_READ_SUPPLEMENT` | On `CallRegionArgs::parity_aligned`, enable read-event supplement. |
  | `GATK_RS_ENABLE_REF_MOTIF` | Inject reference-motif indels when that path has no CIGAR events. |
  | `GATK_RS_ENABLE_CLUSTER_INJECT` | Run `ensure_assembly_cluster_indel_events` outside `strict_java_assembly`. |
  | `GATK_RS_ASM8_ONLY` | Explicit graph-only toggle. Forces `GATK_RS_P12_ENSURE_BRIDGES` off. Production graph-only does not require this variable. |
  | `GATK_RS_HC_GIVEN_VCF` | VCF path of given alleles on the harness path. Production `strict_java` ignores it. |
  | `GATK_RS_DIAGNOSTIC_SKIP_SITE_RESHAPE` | Skip Class-A / A2 / A3 AD/PL reshape. Default off. |
- Several `GATK_RS_*` variables (`P12_PHASE_E`, cluster inject, scaffold
  output, RSS abort, experimental k-best) change or observe a run. The user
  guide lists the ones an operator can hit. Unset is the production
  configuration unless a variable’s own contract says otherwise.
- Compiler warnings for unused items are still present in the haplotypecaller
  crate. Hooks deny rustc warnings except `dead_code`, and they do not deny
  Clippy lints. That is an engineering allowance, not a behavior change.
- Generic file tools remain in the binary for the smoke harness and are
  hidden from `--help`.

Further contracts belong in Rustdoc on the type that owns them, and in the
claim matrix when they are something the project asserts.

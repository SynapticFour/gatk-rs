# Parity Fixtures

Small, deterministic smoke fixtures used by `scripts/parity/run_parity_smoke.sh`.

These files are intentionally tiny and human-auditable. Larger truth datasets
will be added in later parity phases.

## Contents

- `reference.fa`, `reference.fa.fai`, `reference.dict`: minimal `chr1` reference
  (32 bp) plus sidecar index files required by Java GATK for
  `ValidateVariants` and similar tools.
- `sample.sam`: one aligned read with `@RG` including `PL`, and an `RG:Z:` field on the
  alignment line (required for Java `PrintReads` parity — otherwise GATK filters the read).
  Used with `ValidateSamFile` / Rust `Validate` SAM checks.
- `sample.bam`: BAM equivalent of `sample.sam`, written with Rust `PrintReads`
  (HTSlib) for `ValidateSamFile` / Rust `Validate` BAM checks.
- `sample.bam.bai`: BAI index for `sample.bam`, required for interval traversal in
  Java `HaplotypeCaller` (Phase-4 assembly-region parity harness).
- `sample.vcf`: one variant record consistent with `reference.fa` at that locus.
- `regions.interval_list`: interval list for future interval-driven checks.
- `p4_region_boundaries_known_a.tsv`, `p4_region_boundaries_known_b.tsv`:
  tiny per-locus activity probability tracks used by Phase-4 active-region
  deterministic boundary contract tests (steps 55-57 scaffolding).
- `p5_assembly_case1_reads.tsv`: tiny deterministic local-assembly read corpus
  for Phase-5 candidate-set parity contracts (Java debug export vs Rust assembly).
- `p5_equivalence_regions.tsv`: Phase-5 equivalence matrix manifest listing
  region classes/case IDs and fixture-to-expected mappings used by runtime diff reports.
- `p5_live_reference.fa`: tiny unique reference used by the live Java runtime profile.
- `p5_live_case_snp.sam`: small mixed ref/alt SAM corpus on `chrLive` for live
  Java-vs-Rust candidate overlap checks.
- `p5_live_regions.tsv`: manifest for the battery-friendly live runtime matrix.
- `p5_live_regions_extended.tsv`: broader live runtime manifest (edge intervals) for
  non-blocking coverage expansion and drift discovery.
- `p5_live_reference_indel.fa`, `p5_live_case_indel.sam`: indel-oriented corpus with
  mixed `20M`, `10M1I9M`, and `10M1D10M` reads for live candidate overlap checks.
- `p5_live_reference_repeat.fa`, `p5_live_case_repeat.sam`: repeat-rich corpus for
  low-complexity repeat branch behavior.
- `p5_live_reference_lowcomplex.fa`, `p5_live_case_lowcomplex.sam`: low-complexity
  homopolymer-style corpus for fragile-context live differential checks.
- `p4_assembly_region_cases.tsv`: manifest for frozen Java assembly-region interval
  differential checks (`run_p4_active_region_interval_diff.sh`).
- `p6_pairhmm_case1_reads.tsv`: tiny PairHMM likelihood-vector corpus (read, quals,
  mapq, haplotypes) for Phase-6 differential contracts.

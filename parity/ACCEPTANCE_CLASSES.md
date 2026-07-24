# Acceptance Classes (Phase 0 / Step 2)

This policy freezes the acceptance vocabulary used by parity checks.

## Classes

- `exit-only`
  - Contract: Java and Rust must return the same exit class.
  - Typical use: `--help`, `--version`, clearly expected failures.

- `normalized-equivalent`
  - Contract: outputs are compared after deterministic normalization (line extraction, stable sort, whitespace normalization where configured).
  - Typical use: summary/help text, stable scalar reports.

- `sam-file-parity`
  - Contract: SAM header/alignments equivalent under parity normalization (`compare_sam_parity.py`, ignores `@PG/@CO`).
  - Typical use: `PrintReads`-style output checks.

- `vcf-strict`
  - Contract: strict VCF text equivalence after volatile-header allowlist stripping (`compare_vcf_strict.py`).
  - Typical use: locked VCF fixture validation.

- `bam-alignment-parity`
  - Contract: BAM/SAM headers and sorted alignments equivalent via `samtools view` (`compare_bam_alignment_parity.py`).
  - Typical use: BAM output parity.

- `statistical-equivalent`
  - Contract: distributions/metrics must fall within explicitly frozen tolerances and be reproducible under fixed seed and thread profile.
  - Typical use: stochastic or floating-point-sensitive slices.
  - Note: no current required check uses this class yet; it is reserved and now formally defined.

## Profile matrix

- `PARITY_SMOKE_PROFILE=smoke`: default fast gate.
- `PARITY_SMOKE_PROFILE=extended`: includes additional fixture rows (e.g. interval-list and whole-reference count checks).

CI/nightly can run one or both profiles depending on runtime budget.

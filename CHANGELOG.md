# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
once a first numbered release is cut. Until then, versions remain `0.1.x` / Alpha.

## [Unreleased]

### Added

- Canonical mid-B HaplotypeCaller path vs GATK 4.4 (graph through QUAL/QD):
  [`docs/PARITY.md`](docs/PARITY.md),
  [`docs/PARITY_MILESTONE_6R.md`](docs/PARITY_MILESTONE_6R.md).

### Changed

- Scope ADR:
  [`docs/adr/0001-scope-boundary.md`](docs/adr/0001-scope-boundary.md)
  (germline short-variant workflow only).
- Removed stub crate `gatk-tools` (ADR
  [`docs/adr/0002-remove-gatk-tools.md`](docs/adr/0002-remove-gatk-tools.md));
  use samtools/bcftools for generic BAM/VCF ops.
- Product-facing repository surface: README “Why does this exist?”, slim
  [`CONTRIBUTING.md`](CONTRIBUTING.md), and docs reduced to
  [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) +
  [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) (+ `docs/GATK_PINNED.env`).
- Development process docs (L2–L14 sign-offs, audits, sprint plans) retained only
  on git branch `pre-cleanup-archive`.
- Removed `legacy_assembly_graph` feature and orphan root `examples/`.
- Parity dump modules compile only with `--features parity_harness` (or under
  `cfg(test)`); they are not part of the default release build of `gatk-cli`.
- CLI independence / trademark disclaimer; honest `--version` banner (not Broad GATK).

### Fixed

- Active-region reference confidence: empty genotyping evidence returns
  zero-evidence loci (GQ=0/dp=0) without touching the reference FASTA.
- Unit tests no longer pass uninitialized `bam::Record::new()` into pileup paths.
- Mid-B dense-cluster depth-cap and pre-TTC zero-stripe DP expectations aligned
  with production helpers.

### Added

- [`NOTICE.md`](NOTICE.md), [`LICENSE`](LICENSE) (Apache-2.0).
- [`gatk-rs-equiv`](gatk-rs-equiv/) GIAB / differential-fuzz tooling and
  [`tools/equivalence/README.md`](tools/equivalence/README.md) index.

## [0.1.0] — Alpha

Initial Alpha workspace: HaplotypeCaller focused on scoped algorithm parity
against pinned GATK 4.4. See [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) for
what is and is not asserted.

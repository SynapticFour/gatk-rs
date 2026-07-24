# GATK-RS Parity Harness

This folder contains the differential-testing harness used to compare Java GATK
behavior and Rust GATK-RS behavior in a deterministic way.

## Implementation plan (where this fits)

| Phase | Scope | Status |
|-------|--------|--------|
| **0** | CLI/help/version smoke, harness layout, optional “Java missing” skip in CI | Done |
| **1** | File-based Validate parity (SAM/BAM/VCF/FASTA), HTSlib-backed `PrintReads` and BAM validation, richer `diff_outputs` diagnostics, CI against real GATK in Docker | Done |
| **2+** | HaplotypeCaller engine + **scoped** output parity vs Java (P12 / L2 / dense gates) | **Signed on gates** — see [`docs/CLAIM_MATRIX.md`](../docs/CLAIM_MATRIX.md). Not a genome-wide bitwise lock. |

**Current checkpoint:** HC substrate is implemented; assert only what the claim matrix lists as Yes.

## Layout

- `fixtures/`: small deterministic fixtures for smoke parity checks
- `reports/`: generated comparison reports (JSON + markdown summary)

## Deterministic runner conventions

The parity scripts force a deterministic environment:

- `LC_ALL=C`
- `TZ=UTC`
- `RUST_LOG=error`
- stable output directories under `parity/reports/`

## Required tools

- Rust toolchain and workspace dependencies
- Java runtime (for running Java GATK)
- One of:
  - `GATK_JAR` environment variable (path to `gatk-package-*-local.jar`), or
  - `gatk` on PATH, or
  - `GATK_DOCKER_IMAGE` (for example `us.gcr.io/broad-gatk/gatk:4.4.0.0`)

Optional for Docker:

- `GATK_DOCKER_PLATFORM` (for example `linux/amd64` on Apple Silicon hosts)

## Commands

Run the smoke parity suite:

```bash
./scripts/parity/run_parity_smoke.sh
```

This generates:

- JSON + Markdown summary under the gitignored `parity/reports/` output directory

## Current scope

Smoke checks validate wrapper behavior and argument-level parity for:

- `--help`
- `--version`
- `HaplotypeCaller --help`
- `PrintReads --help`

The suite also includes initial content and mapped-command parity checks:

- `--version` banners (Java GATK string vs Rust independent-project string; presence check, not identical text)
- `HaplotypeCaller --help` summary text
- `PrintReads --help` tool-name presence
- mapped help parity: Java `ValidateSamFile --help` vs Rust `Validate --help`
- invalid-interval user error parity for `HaplotypeCaller -L chr999:1-10`

File-based parity for validation-style runs:

- Matching **exit codes** (both must succeed).
- Matching **success signals** in stdout via `diff_outputs.py` (`--extract-regex` +
  `--presence-only`): Picard/GATK log lines differ from Rust, so the harness checks
  that each side emits an expected success substring (for example `No errors found`
  vs `SAM validation passed`, `Processed N total variants` vs `VCF validation passed`,
  and `Processed N total bases` vs `FASTA validation passed`).

Concrete tool pairs:

- Java `ValidateSamFile` on `parity/fixtures/sample.sam` vs Rust `Validate` (`-t SAM`) with the same reference
- Java `ValidateSamFile` on `parity/fixtures/sample.bam` (binary BAM produced from the SAM via Rust `PrintReads`) vs Rust `Validate` (`-t BAM`) with the same reference
- Java `ValidateVariants` on `parity/fixtures/sample.vcf` vs Rust `Validate` (`-t VCF`)
- Java `CountBasesInReference` on `parity/fixtures/reference.fa` vs Rust `Validate` (`-t FASTA`)

When a check uses `--extract-regex`, per-check JSON from `diff_outputs.py` also includes
`java_extracted` and `rust_extracted` (post-regex, capped) to simplify debugging failures.

**PrintReads file parity:** after the stdout-based checks, the harness runs Java and Rust
`PrintReads` on `parity/fixtures/sample.sam` and compares the written SAM files with
`scripts/parity/compare_sam_parity.py` (header `@SQ`/`@RG`/`@HD` and alignment lines, ignoring
`@PG`/`@CO`). Alignment lines must include an `RG` tag so Java GATK does not drop reads under
`WellformedReadFilter`.

As tool implementations mature, this harness is designed to add richer
behavioral comparisons for each tool incrementally.

## Foundation Gate (Phase 0/1)

Use the canonical foundation gate runner (required checks first):

```bash
./scripts/parity/run_foundation_gate.sh
```

Config lives in `parity/checks.json` with:

- `required`: hard-gating checks
- `advisory`: non-gating checks (run with `FOUNDATION_RUN_ADVISORY=1`)
- per-check metadata: `acceptance_class`, `timeout_s`, `owner`

This keeps the baseline explicit and CI-aligned while we iterate toward broader parity matrices.

# Canonical parity targets (Phase 0 / Step 1)

**Status:** historical Phase 0 target list. It is not the current claim set.
The P1 row that says HaplotypeCaller VCF body parity is “future” is stale
relative to the signed P12 and L2 rows. Use [`docs/CLAIM_MATRIX.md`](../docs/CLAIM_MATRIX.md)
and [`docs/PARITY.md`](../docs/PARITY.md) for what is asserted now. This file
remains the original P0 surface so the smoke harness still has a written
target list. Do not add a tool here without updating the claim matrix.

Frozen list for **strict parity-first** work. Scope is **non-Spark HaplotypeCaller** and **shared I/O** used on the HC path.

Why this list is retained:
- It is the original Phase 0 smoke-harness target list, so those checks still have a written surface.
- It is not the current claim set. What the project asserts now is [`docs/CLAIM_MATRIX.md`](../docs/CLAIM_MATRIX.md).
- Changing this historical list still requires a claim-matrix update if that change would make the project assert something new. See also [`docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md).
- Spark tools and non-HC workflows were outside this Phase 0 list.

| Priority | Java (GATK / Picard) | Rust (`gatk-rs`) | Notes |
|----------|----------------------|------------------|--------|
| P0 | CLI `--help`, `--version` | Same surface | Exit codes + banner where applicable |
| P0 | `ValidateSamFile`, `ValidateVariants`, `CountBasesInReference` | `Validate` (SAM/BAM/VCF/FASTA) | Mapped checks in `run_parity_smoke.sh` |
| P0 | `PrintReads` | `PrintReads` | SAM file parity (`compare_sam_parity.py`) |
| P0 | `HaplotypeCaller` (args, intervals, failures) | `HaplotypeCaller` | Invalid interval / missing input parity; **no** bitwise VCF yet |
| P1 | `HaplotypeCaller` VCF body | HC pipeline | **Future** — golden under `parity/expected/` when ready |

## Explicit out-of-scope (until Phase 3+)

- Spark tool variants and distributed execution semantics.
- End-to-end genotyping/annotation parity outputs beyond Phase-2 foundation contracts.
- Large truth-set benchmarking and performance dashboards (Phase 10).

Changes to this table require a PR that updates [`docs/CLAIM_MATRIX.md`](../docs/CLAIM_MATRIX.md).

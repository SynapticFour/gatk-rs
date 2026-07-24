# Canonical parity targets (Phase 0 / Step 1)

Frozen list for **strict parity-first** work. Scope is **non-Spark HaplotypeCaller** and **shared I/O** used on the HC path.

Freeze policy:
- This list is the canonical Step-1 target surface for P0/P1/P2 completion.
- Any addition/removal requires updating [`docs/CLAIM_MATRIX.md`](../docs/CLAIM_MATRIX.md) and [`docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md) in the same change.
- Spark tools and non-HC workflows are explicitly out of scope until Phase 3+.

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

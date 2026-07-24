# ADR 0002: Remove the `gatk-tools` stub crate

- **Status:** Accepted
- **Date:** 2026-07-23
- **Related:** [ADR 0001 — scope boundary](0001-scope-boundary.md)
- **Decision:** Option **(a)** — delete `gatk-tools`; do not ship a generic BAM/VCF/FASTA toolkit.

## Context

`gatk-tools` contained only stub I/O (`bam.rs`, `fasta.rs`, `vcf.rs`, `io.rs`) that
returned “not yet implemented”. Nothing in the production path depended on it
(`gatk-cli` / `gatk-haplotypecaller` use `gatk-core`). Stub crates that look like
product surface erode trust more than they communicate intent.

Two options were considered:

| Option | Summary |
|--------|---------|
| **(a)** Delete `gatk-tools` | Scope the project to HC + joint-genotyping/filtering; point users at samtools/bcftools for generic file ops |
| **(b)** Slim `gatk-tools` | Keep only workflow-needed helpers (e.g. VCF sort/index) and delete the rest |

## Decision

Choose **(a)**.

Rationale:

1. **Credibility** — empty TODOs in a crate named like a toolkit read as unfinished product, not honest scope.
2. **Capacity** — maintaining generic BAM/VCF utilities duplicates battle-tested **samtools** / **bcftools**; effort belongs in joint genotyping, filtering, and HC parity.
3. **Existing I/O** — workflow-needed parse/write already lives in `gatk-core` and is used by HC / CombineGVCFs / GenotypeGVCFs / VariantFiltration.
4. **Alignment with project posture** — this is not a Broad GATK drop-in toolkit; it is an HC-focused algorithm-parity project (see [ADR 0001](0001-scope-boundary.md) and `docs/CLAIM_MATRIX.md`).

Option (b) was rejected for now: any sort/index needs that arise can call out to
samtools/bcftools in scripts, or land as thin helpers next to the consuming tool
in `gatk-core` if they are truly HC-workflow-specific — not as a parallel
toolkit crate.

## Consequences

- Workspace member `gatk-tools` removed; Docker/CI/publish steps updated.
- README and architecture docs state the scope boundary explicitly.
- `scripts/parity/deferred_features_audit.py` asserts the crate stays gone and
  that this ADR remains the record of the decision.

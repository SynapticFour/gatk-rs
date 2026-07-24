# ADR 0001: Scope boundary — germline short-variant workflow, not a full GATK4 port

- **Status:** Accepted
- **Date:** 2026-07-23
- **Decision:** gatk-rs is **intentionally not** a complete GATK4 reimplementation. Product
  effort concentrates on the germline short-variant calling and cohort workflow:
  **HaplotypeCaller → CombineGVCFs → GenotypeGVCFs → hard-filtering**
  (`VariantFiltration`).

## Context

GATK4 is a large toolkit spanning pre-processing, germline and somatic calling,
CNV/SV, recalibration, and annotation. Porting all of it would dilute algorithm
parity work on the path that most clinical/research germline short-variant
pipelines actually run after alignment. This project already asserts scoped
parity (see [`docs/CLAIM_MATRIX.md`](../CLAIM_MATRIX.md)), not genome-wide
drop-in equivalence.

A clear “in / out of scope” record prevents stub crates, fake CLI surface, and
scope creep that looks like an unfinished Broad clone.

## Decision

**In scope (product focus):**

| Stage | Tool | Role |
|-------|------|------|
| Call | `HaplotypeCaller` | Germline SNPs/indels via local reassembly (pinned GATK 4.4 algorithm parity on agreed gates) |
| Merge gVCFs | `CombineGVCFs` | Multi-sample reference-confidence merge |
| Joint genotype | `GenotypeGVCFs` | Cohort AF / QUAL / genotypes from combined gVCF |
| Filter | `VariantFiltration` | GATK-style hard filters (not VQSR) |

Supporting I/O for that workflow lives in `gatk-core`. Generic file utilities
are out of scope ([ADR 0002](0002-remove-gatk-tools.md): use samtools/bcftools).

**Not planned (with rationale):**

| Area | Why out of scope |
|------|------------------|
| **BaseRecalibrator / BQSR** | Pre-call BAM recalibration. Users should bring BAMs already recalibrated (or calibrated by their aligner/pipeline). Reimplementing BQSR would own a large, separate math surface without improving HC joint-genotyping parity. |
| **VariantRecalibrator / VQSR** | Needs a trained Gaussian mixture model and large cohorts. Hard-filtering (`VariantFiltration`) is the pragmatic, GATK-documented fallback for smaller callsets — see README / `CLAIM_MATRIX` (VQSR not asserted). |
| **Mutect2 / somatic analysis** | Different likelihood model, tumor–normal semantics, and evaluation culture. Worth a **separate** project if pursued — not a side quest inside germline HC. |
| **gCNV / structural & copy-number variants** | Entirely different algorithms (read-depth / segmentation / graph SV), not a variant of HaplotypeCaller local assembly. |
| **Funcotator / functional annotation** | Downstream of calling. Mature external tools already exist (**VEP**, **snpEff**, etc.); duplicating them does not advance HC algorithm parity. |

## Consequences

- README and architecture docs lead with this boundary; new tools must justify
  fit to the germline short-variant spine or be rejected as out of scope.
- Equivalence claims stay tied to HC / Combine / Genotype / hard-filter gates
  in `CLAIM_MATRIX.md` — not “full GATK4”.
- Requests for BQSR, VQSR, Mutect2, gCNV/SV, or Funcotator are answered by
  pointing here (and to external tools where appropriate), not by adding stubs.

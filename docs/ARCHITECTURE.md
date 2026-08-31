# Architecture

gatk-rs is a native Rust workspace focused on a **GATK 4.4–aligned HaplotypeCaller**.
It is an independent community project (see [`NOTICE.md`](../NOTICE.md)). Product claims
live only in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). Canonical mid-B HC Java 4.4
contracts: [`PARITY.md`](PARITY.md).

## Workspace layout

| Crate / path | Role |
|--------------|------|
| [`gatk-cli`](../gatk-cli/) | `gatk-rs` binary — HaplotypeCaller and workflow tools |
| [`gatk-haplotypecaller`](../gatk-haplotypecaller/) | HC engine: activity → assembly regions → PairHMM → genotyping → VCF/gVCF emit |
| [`gatk-core`](../gatk-core/) | BAM/VCF/FASTA I/O, intervals, reference helpers, VariantFiltration |
| [`gatk-common`](../gatk-common/) | Shared errors / config helpers |
| [`gatk-rs-equiv`](../gatk-rs-equiv/) | GIAB / hap.py / vcfeval equivalence + differential fuzz driver |
| [`fuzz/`](../fuzz/) | LibFuzzer target sharing scenarios with `gatk-rs-equiv` |
| [`scripts/parity/`](../scripts/parity/) | L2/P12/GIAB harness scripts that back the claim matrix |
| [`parity/fixtures/`](../parity/fixtures/) | Tracked synthetic / P12 fixtures for gates |

Index: [`tools/equivalence/README.md`](../tools/equivalence/README.md).

## What is implemented

Default CLI HaplotypeCaller runs the production path:

`CallRegionArgs::strict_java()` → assembly-region pipeline → variant / gVCF emission.

That path covers local reassembly, PairHMM likelihoods, genotyping, and core INFO/FORMAT
fields needed for the signed P12 / L2 / dense-window gates in the claim matrix.

## Post-call tools (beyond HC)

| Tool | Role |
|------|------|
| `CombineGVCFs` | Merge per-sample gVCFs |
| `GenotypeGVCFs` | Joint genotype a combined gVCF |
| `VariantFiltration` | Soft hard-filters via FILTER tags (`gatk-core::variant_filtration`) |

**VariantFiltration vs VQSR:** hard-filtering does **not** replace VQSR
algorithmically. It is the pragmatic fallback for smaller cohorts where VQSR
cannot be trained cleanly — the same recommendation GATK publishes when VQSR is
out of reach. Official SNP/indel expression tables live in
`gatk_core::variant_filtration::{GATK_HARD_FILTER_SNP, GATK_HARD_FILTER_INDEL}`.

## What is not implemented (or not product)

- Genome-wide clinical drop-in equivalence to Broad GATK (not asserted).
- Joint multi-sample HC across merged `-I` inputs.
- **VQSR** (Gaussian mixture recalibration) — use `VariantFiltration` hard filters instead for small cohorts.
- A generic BAM/VCF/FASTA **toolkit** crate (`gatk-tools` removed; see [`adr/0002-remove-gatk-tools.md`](adr/0002-remove-gatk-tools.md)). Use **samtools** / **bcftools** for sort, index, and generic file ops.
- Broader GATK4 surface (BQSR, VQSR, Mutect2, gCNV/SV, Funcotator) — see [`adr/0001-scope-boundary.md`](adr/0001-scope-boundary.md).
- bamout, DRAGSTR calibration, DRAGEN mode, allele-specific `AS_*`, Java `--assembly-region-out`.
- Bitwise-identical QUAL/FORMAT everywhere.

Feature flag `parity_harness` (includes `dev-dumps`) exposes dump/oracle surfaces used by
L2 scripts; it is **not** part of the default release CLI surface. Pure `*_dump` modules
compile only with `dev-dumps` / `parity_harness` (or under `cfg(test)`). The CLI enables
`dev-dumps` so `DumpSmoothedActivity` keeps working.

Harness env flags (ignored unless built with `--features parity_harness`):
`P12_PHASE_E`, `P12_BASELINE_EMIT_FILTER`, `GATK_RS_P12_EVENT_REGISTRY`,
`GATK_RS_P12_ENSURE_BRIDGES`, `P12_L4_JAVA_FORMAT`, `GATK_RS_ENABLE_READ_SUPPLEMENT`,
`GATK_RS_ENABLE_REF_MOTIF`, `GATK_RS_ENABLE_CLUSTER_INJECT`, `GATK_RS_ASM8_ONLY`,
`GATK_RS_HC_GIVEN_VCF`.

## Equivalence proof

Scientific evidence is runnable, not narrative:

1. **Unit / integration tests** — `cargo test --workspace`
2. **L2 synthetic + P12 real-window gates** — `scripts/parity/` (CI workflows under `.github/workflows/`)
3. **GIAB equivalence** — `gatk-rs-equiv` + `scripts/parity/giab/`
4. **Differential fuzz** — `gatk-rs-equiv differential-fuzz` / `fuzz/run_hc_differential.sh`

Pinned Java oracle: [`GATK_PINNED.env`](GATK_PINNED.env) (GATK 4.4.0.0).

## Design posture

Prefer Rust-native modules and algorithm parity with the pinned Java behavior over cloning
Java class trees. Observable contracts and waivers are recorded in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md).
Further detail belongs in Rustdoc and code comments, not additional markdown sprawl.

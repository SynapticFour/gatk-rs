# Architecture

gatk-rs is a native Rust workspace focused on a **GATK 4.4–aligned HaplotypeCaller**.
It is an independent research reimplementation (see [`NOTICE.md`](../NOTICE.md)). Product claims
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
`GATK_RS_HC_GIVEN_VCF`, `GATK_RS_DIAGNOSTIC_SKIP_SITE_RESHAPE`.

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

SeqGraph k-best **K** (completed paths) is separate from the live-frontier **resource
policy**. Production default remains `legacy_1024`. Experimental env
`GATK_RS_EXPERIMENTAL_KBEST_POLICY` (`unbounded_diagnostic` / `byte_budget`) is not a
product contract: [`parity/KBEST_RESOURCE_POLICY.md`](parity/KBEST_RESOURCE_POLICY.md).
Chr20_tiny 6R.130–6R.185 (default USE_PLS GT/PL/GQ, first AD write after
retainEvidence, Class-A3 skip of pileup GT replacement so an assigned
calculator genotype is preserved; hom-ref calculator 0/0 is not
VCF-emitted on either diagnostic path; `2:92316347` FORMAT AD/PL now
matches Java after removing Rust-only dangling-merge haplotype
materialization; 6R.165 proved INFO FS/MQ/SOR still used `region.reads`
pileup rather than Java `AlleleLikelihoods`; 6R.166 switched FS/SOR onto
that post-filter likelihood table; 6R.167 switched MQ onto Java
`sampleEvidence`; 6R.168 applied Java RMS so MQ=40.25; 6R.169 inventoried
Rust-only `ReadPosRankSum`/`InbreedingCoeff` as Java-suppressed default
annotators; 6R.170 proved ReadPosRankSum pileup vs informative
likelihoods; 6R.171 proved empty/NaN RankSum was collapsed to `0.0`;
6R.172 restored `Option<f64>` so finite `0.0` still emits and undefined
omits; 6R.173 inventoried the next genuine remaining VCF field as INFO
DP at `2:92305634`; 6R.174 switched INFO DP onto Java
`Coverage.evidenceCount` so that site is 3 while FORMAT DP stays 2;
6R.175 proved remaining SOR 0.693 vs 1.179 is extra informative
membership of mate-on-other-contig FLAG=145, not `calculateSOR`;
6R.176 applied that mate-contig gate at PairHMM construction so
SOR is 0.693 and INFO DP stays 3;
6R.177 proved remaining InbreedingCoeff=1.0 is always-insert vs Java
`MIN_SAMPLES=10` emptyMap, not a sample-count mismatch;
6R.178 gated record emission on `n_genotypes >= 10` and left the
`1 - het/n` formula stacked vs Java HWE;
6R.179 proved the next INFO split at `2:92307324 TTC/T` is region-wide
vs per-variant annotation membership for DP/MQ/SOR;
6R.180 bound INFO DP/MQ/SOR to that per-variant genotyping
`AlleleLikelihoods` so the site is DP=1 MQ=44 SOR=1.609;
6R.181 proved `2:92307333 T/G` still falls back to region-wide
evidence because the cluster-TG path never constructs the
per-variant object;
6R.182 proved Java’s object is `filterPoorlyModeledEvidence` 8→2 then
`marginalize`+`retainEvidence` 2→1, so overlap-of-six is not that
object;
6R.183 proved Rust’s poorly-modeled pass already matches Java n=2 and
the stored n=8 is an unfiltered P12-cluster refresh overwrite;
6R.184 proved that last refresh is the final membership write and that
replaying the existing normalize+filter on it restores Java n=2 then
diagnostic `retainEvidence` n=1;
6R.185 restored that Java-order lifecycle after the last P12 refresh
so stored evidence is n=2, without constructing `annotation_likelihoods`)
is engineering
discovery in [`PARITY.md`](PARITY.md), not a claim-matrix Yes row.

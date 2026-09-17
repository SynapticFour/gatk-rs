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
Chr20_tiny 6R.130–6R.221 (default USE_PLS GT/PL/GQ, first AD write after
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
so stored evidence is n=2, without constructing `annotation_likelihoods`;
6R.186 constructed that per-variant object on the cluster-TG
early-template path from the stored n=2 matrix via `marginalize` then
`retainEvidence`, so INFO DP=1 MQ=44 SOR=1.609;
6R.187 proved the next remaining INFO field is DP at `2:92305635 A/G`
from a genotyping-subset object n=1 versus Java stored-hap
`retainEvidence` n=3, not the 6R.186 empty-object cause;
6R.188 proved that n=1 object is the FORMAT-narrowed SiteScore
subset, and that 6R.186 cluster-TG construction is not selected;
6R.189 bound SiteScore annotation to the stored-hap loc-loop object
so `2:92305635 A/G` is INFO DP=3 MQ=41.96 SOR=0.693;
6R.190 proved the next remaining INFO field is rust-only
`ReadPosRankSum` at `2:92305716 A/C` from `region.reads` pileup
versus Java `fillQualsFromLikelihood` on the annotation
`AlleleLikelihoods`;
6R.191 bound ReadPosRankSum to that annotation object so the
target omits the key (REF=0 ALT=3, 6R.172 undefined semantics);
6R.192 proved the next remaining INFO field is Java-only
`BaseQRankSum`/`MQRankSum` at `2:92307359 CT/C` (header-only emit
predicate, not another RankSum source-object or QD formula);
6R.193 emitted BaseQRankSum from that same fillQuals object
(REF=1 ALT=1 → finite `0.0`; MQRankSum still header-only);
6R.194 emitted MQRankSum from that same fillQuals membership
(REF=1 ALT=1 unequal MAPQ → finite `-0.674`; QUAL/QD not retuned);
6R.195 proved QD at `2:92307359` is QUAL/`getDepth` with no
independent formula or rounding arrow (Java 15.80 = Java QUAL/2;
Rust 15.82 = Rust QUAL/2; production unchanged);
6R.196 bound biallelic QUAL AF to Java’s length-based alt prior so
that indel is QUAL 31.60 / QD 15.80 while SNPs with the same PL stay
31.64;
6R.197 proved the next remaining INFO field is DP at `2:92307403 C/A`
(Java 6 vs Rust 4; Java-only BaseQ/ReadPos; MQRankSum already
matches; 6R.196 QUAL pin stays 31.60; production unchanged);
6R.198 proved that site is the cluster-downstream early-template
with empty `annotation_likelihoods` (stored n=6, diagnostic retain
n=4; Java INFO DP=6 is evidenceCount of the n=6 stored set, not the
6R.186 retain helper; production unchanged);
6R.199 attached that stored-unique n=6 object on the
cluster-downstream early-template (not the 6R.186 retain subset;
FORMAT/QUAL unchanged; INFO DP still 4 via Coverage overlap,
deferred to 6R.200);
6R.200 proved that remaining DP 6 vs 4 is a second Coverage overlap
filter on the n=6 object, not a wrong source object; Java
`Coverage.annotate` is `evidenceCount()` with no overlap;
production unchanged;
6R.201 removed that second Coverage overlap so INFO DP is the
attached unique evidenceCount (n=6) at `2:92307403 C/A`;
FORMAT/QUAL/MQ unchanged; RankSums not retuned;
6R.202 proved remaining Java-only BaseQ/ReadPos is empty-REF
fillQuals `getElementForRead` on the two extra stored unique rows
(invocation and n=6 object already reach RankSum; MQRankSum already
matches; production unchanged);
6R.203 retried the pre-realign covering CIGAR in shared
`getElementForRead` so those REF rows yield BaseQ 30,30 and
ReadPos 65,91 (`BaseQRankSum=-1.834`, `ReadPosRankSum=1.282`);
MQRankSum stays 1.834; INFO DP stays 6; FORMAT/QUAL unchanged;
6R.204 proved the next INFO DP split at `2:92316296 A/T` is
region-wide stored unique n=3 versus Java loc-loop retainEvidence
n=2 because the two-read hom-alt early-template attaches no
annotation object; production unchanged;
6R.205 attached the stored-haplotype loc-loop object on that
two-read hom-alt arm so annotation is retainEvidence n=2
(INFO DP=2, MQ=47.00, SOR=2.303; FORMAT/QUAL unchanged);
6R.206 proved the next INFO DP split at `2:92316416 C/A` is
region-wide stored unique n=3 versus Java loc-loop retainEvidence
n=1 because the one-read hom-alt early-template attaches no
annotation object; production unchanged;
6R.207 attached the stored-haplotype loc-loop object on that
one-read hom-alt arm so annotation is retainEvidence n=1
(INFO DP=1, MQ=21.00, SOR=1.609; FORMAT/QUAL unchanged);
6R.208 proved the next INFO DP/SOR split at `2:92317399 C/A` is
Java loc-loop retainEvidence n=2 versus Rust attached n=1 after
the 6R.180 same-QNAME collapse; production unchanged;
6R.209 skipped that collapse on the loc-loop object so both
same-QNAME mates survive (INFO DP=2, SOR=0.693; FORMAT/QUAL/MQ
unchanged; 6R.180/6R.186 stay on FORMAT-subset / stored-unique);
6R.210 attached the stored-haplotype loc-loop object on the
gap-sparse shaped-early FORMAT path so annotation is retainEvidence
n=1 (INFO DP=1, MQ=24.00, SOR=1.609; MQRankSum omitted; FORMAT/QUAL
unchanged);
6R.211 attached the same object on the gap-tail het early-template
FORMAT path so annotation is retainEvidence n=3 at `2:92325193 C/T`
(INFO DP=3, MQ=28.03, SOR=0.223; FORMAT/QUAL unchanged);
6R.212 attached the same object on the weak-sparse het early-template
FORMAT path so annotation is retainEvidence n=3 at `2:92325268 C/T`
(INFO DP=3, MQ=28.03, SOR=1.179; FORMAT/QUAL unchanged);
6R.213 proved the next remaining common-site split is FORMAT PL at
`20:29455015 G/T` (Java `69,0,2140` vs Rust `122,0,2304`); QUAL/GQ/QD
are downstream of those GLs; GT/AD/INFO match; production unchanged;
6R.214 proved the pre-marginalization haplotype matrix at that site is
already a different object — Java 30 trimmed 130 bp haplotypes (mapper
18/12) versus Rust 84 of 181/212 bp (mapper 66/18, 13 extra window
sequences); stop before `marginalize`; production unchanged;
6R.215 proved the first extra haplotypes are materialized by
`merge_rt_kbest_pre_remove_paths` RT k=10 after SeqGraph k=25 already
matches Java's 78; production unchanged;
6R.216 proved Java never re-enters the ReadThreadingGraph after
SeqGraph `findBestPaths` and `createGraph(k=10)` returns null on
cycles, so the 44 k=10 paths have no Java counterpart; global RT
disable is too broad; production unchanged;
6R.217 proved a test-only `createGraph` cycle abort at RT extract
drops those 44 and restores unique merge=78 while p11/indel4 acyclic
k=10 ALTs remain; production extract still abort=false;
6R.218 applied that abort at SeqGraph-path RT extract so assemble=78,
hap_n=30, and PL/GQ/QUAL/QD match Java at `20:29455015 G/T`;
6R.219 proved the next FORMAT split is AD at `20:29455379 G/A`
(Java 42,5 vs Rust 44,5) because FORMAT is not the Java
retainEvidence `DepthPerAlleleBySample` object;
6R.220 proved the 52-row remarg is isolated `try_genotype`
(`AD 47,5`) while `call_region` FORMAT is a later object
(`AD 44,5`); no BAM-level 52→49 filter; Case A; production
unchanged;
6R.221 proved `assign_genotype_likelihoods_for_region` is that
FORMAT object while production-arg `try_genotype` stays 47,5;
colocated merge does not fire; annotation n=52 both; production
unchanged).
This is engineering
discovery in [`PARITY.md`](PARITY.md), not a claim-matrix Yes row.

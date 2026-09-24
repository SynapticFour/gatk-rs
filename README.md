# gatk-rs

[![CI](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

Built by **[Synaptic Four](https://synapticfour.com)** · [contact@synapticfour.com](mailto:contact@synapticfour.com).

gatk-rs is an independent research reimplementation (single-steward) and is not
affiliated with, endorsed by, or supported by the Broad Institute.
"GATK" is a trademark of the Broad Institute; this project's name and
branding will be revisited if requested. Parity tests **call** a pinned
GATK 4.4 jar as an oracle; they do not ship Broad source (see [`NOTICE.md`](NOTICE.md)).

> **Maturity: Alpha** — validated on limited genomic regions and fixtures, not as a
> genome-wide clinical drop-in. Independent research reimplementation; not Broad;
> not required for the Ferrum join. Authoritative claims and non-claims:
> [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). Trademark and third-party data notes:
> [`NOTICE.md`](NOTICE.md).

### One screen: what this proves / does not prove

| Proves (on `main`) | Does **not** prove |
|--------------------|--------------------|
| Germline spine vs pinned **GATK 4.4**: HC → CombineGVCFs → GenotypeGVCFs → hard `VariantFiltration` | Full GATK4 toolkit (BQSR, VQSR, Mutect2, gCNV/SV, Funcotator, …) |
| **L2** synthetic gates (223/223) + **P12** L3/L4/L5 on `chr2:92300000–92350000` | Genome-wide / full-autosome HC equivalence |
| Canonical **mid-B** HC path (assembly through QUAL/QD) on `2:92317000-92319000` — [`docs/PARITY.md`](docs/PARITY.md) | Whole-codebase HaplotypeCaller equivalence |
| Synthetic joint-genotype cohort ladder (**≤100** samples on a tiny interval) | WGS × large-N / GenomicsDB-class joint calling |
| Scoped algorithm parity with honest waivers (W-H1 / W-H3 / W-L7-FORMAT) | Clinical drop-in, bitwise-identical QUAL/FORMAT everywhere, or a product launch |
| Equivalence **harness** green on GIAB **smoke** (hosted CI, RTG F1 Δ=0 on three ~50 kb windows) | Signed vs unsigned GIAB / autosome rows — **[`CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) only**; do not infer from this table |

Authority: [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). Canonical mid-B is one
ActiveFull region ([`docs/PARITY.md`](docs/PARITY.md)), not a genome-wide product claim.
Chr20_tiny 6R.130–6R.265 holdouts in that file are engineering discovery through
default USE_PLS GT/PL/GQ, the first AD write, Class-A3 preservation of
an assigned calculator genotype (not a `pl_gt == 0` guard), the
hom-ref emission boundary (Java also does not VCF-emit that 0/0),
closure of the `2:92316347` FORMAT AD/PL split by not materializing
dangling path_bases as a standalone haplotype, the 6R.165
proof-only FS/MQ/SOR annotation-membership split, the 6R.166
FS/SOR evidence-source fix, the 6R.167 MQ evidence-source
fix, the 6R.168 RMS formula so MQ=40.25, the 6R.169
proof that remaining INFO extras are Java-suppressed
`ReadPosRankSum`/`InbreedingCoeff`, and the 6R.170
proof that ReadPosRankSum still uses pileup versus
Java informative likelihoods, and the 6R.171 proof that
undefined RankSum is collapsed to numeric `0.0` (Java omits the key;
a legitimate finite `0.0` must remain emittable), and the 6R.172
change to `Option<f64>` so undefined omits while finite `0.0` still
emits, and the 6R.173 reconnaissance that the next genuine remaining
VCF field is INFO DP at `2:92305634` (Java 3 vs Rust 2), and the 6R.174
switch of INFO DP onto Java `Coverage.evidenceCount` (FORMAT DP stays 2),
and the 6R.175 proof that remaining SOR at that site is extra
informative membership of a mate-on-other-contig read, not the SOR
formula, and the 6R.176 PairHMM mate-contig membership so FLAG=145
is `addEvidence(0)` only (SOR 0.693, INFO DP stays 3),
and the 6R.177 proof that InbreedingCoeff=1.0 is always-insert versus
Java `MIN_SAMPLES=10` emptyMap, and the 6R.178 emission gate so
one-sample records omit InbreedingCoeff (formula still `1 - het/n`),
and the 6R.179 proof that `2:92307324 TTC/T` INFO DP/MQ/SOR still used
the region-wide PairHMM matrix versus Java’s per-variant n=1
annotation likelihoods, and the 6R.180 bind of those INFO fields onto
that per-variant object (DP=1 MQ=44 SOR=1.609),
and the 6R.181 proof that `2:92307333 T/G` still has no
per-variant annotation `AlleleLikelihoods` on the cluster-TG
early-template path,
and the 6R.182 proof that Java constructs that object as
`filterPoorlyModeledEvidence` 8→2 then `marginalize`+`retainEvidence`
2→1 (direct dump n=1 FLAG=83 MAPQ=44; Rust overlap candidate is n=6),
and the 6R.183 proof that Rust’s poorly-modeled pass already drops to
n=2 and stored n=8 is a later unfiltered cluster refresh,
and the 6R.184 proof that replaying the existing normalize+filter
on that last refreshed matrix restores Java n=2 then diagnostic
`retainEvidence` n=1,
and the 6R.185 restore of that Java-order filter after the last
P12 refresh (stored n=2),
and the 6R.186 construction of the cluster-TG per-variant annotation
object (n=2 → `marginalize` → `retainEvidence` n=1; INFO DP=1 MQ=44
SOR=1.609),
and the 6R.187 proof-only reconnaissance that the next remaining
INFO field is DP at `2:92305635 A/G` (attached genotyping-subset
n=1 vs Java stored-hap `retainEvidence` n=3; not the 6R.186 cause),
and the 6R.188 proof that the n=1 object is the FORMAT-narrowed
SiteScore subset (not 6R.186 cluster-TG construction),
and the 6R.189 production bind of SiteScore annotation to the
stored-hap loc-loop object (INFO DP=3 MQ=41.96 SOR=0.693 at
`2:92305635 A/G`; FORMAT/QUAL unchanged),
and the 6R.190 proof-only reconnaissance that the next remaining
INFO field is rust-only `ReadPosRankSum` at `2:92305716 A/C`
(pileup vs Java `fillQualsFromLikelihood`; FORMAT/QUAL/DP/MQ/SOR
closed; 6R.189 site stays closed),
and the 6R.191 production bind of ReadPosRankSum to that
annotation `AlleleLikelihoods` object (`2:92305716 A/C` omits
the key; REF=0 ALT=3; formula unchanged),
and the 6R.192 proof-only reconnaissance that the next remaining
INFO field is Java-only `BaseQRankSum`/`MQRankSum` at
`2:92307359 CT/C` (header-only; not another RankSum source or QD),
and the 6R.193 production emit of BaseQRankSum from that same
fillQuals object (`2:92307359 CT/C` finite zero; MQRankSum still
header-only),
and the 6R.194 production emit of MQRankSum from that same
fillQuals membership (`2:92307359 CT/C` finite `-0.674`; QUAL/QD
not retuned),
and the 6R.195 proof-only that QD at that site follows QUAL/AD-depth
(Java 15.80 = Java QUAL/2; Rust 15.82 = Rust QUAL/2; no independent
QD arrow; production unchanged),
and the 6R.196 production bind of biallelic QUAL AF to Java’s
length-based alt prior (`2:92307359 CT/C` QUAL 31.60; QD follows),
and the 6R.197 proof-only reconnaissance that the next remaining
INFO field is DP at `2:92307403 C/A` (Java 6 vs Rust 4; Java-only
BaseQ/ReadPos; production unchanged),
and the 6R.198 proof-only that this is the cluster-downstream
early-template skip (empty annotation object; stored n=6 vs
retain n=4; production unchanged),
and the 6R.199 production attach of stored-unique n=6 annotation
evidence on that arm (not the 6R.186 retain subset; INFO DP still
4 is deferred),
and the 6R.200 proof-only that INFO DP 6 vs 4 is a second Coverage
overlap filter on that n=6 object (production unchanged),
and the 6R.201 production removal of that second Coverage overlap
so INFO DP consumes attached unique evidenceCount (n=6) at
`2:92307403 C/A` (FORMAT/QUAL/MQ unchanged; RankSums not retuned),
and the 6R.202 proof-only that remaining Java-only BaseQRankSum /
ReadPosRankSum is empty-REF fillQuals `getElementForRead` on the two
extra stored unique rows (invocation and n=6 object already reach
RankSum; MQRankSum already 1.834; production unchanged),
and the 6R.203 production shared RankSum `getElementForRead`
covering-CIGAR retry at `2:92307403 C/A` so BaseQRankSum=−1.834
and ReadPosRankSum=1.282 (MQRankSum stays 1.834; INFO DP stays 6;
FORMAT/QUAL unchanged),
and the 6R.204 proof-only that remaining INFO DP at `2:92316296 A/T`
is Java loc-loop retainEvidence n=2 versus Rust region-wide stored
unique n=3 (empty two-read hom-alt annotation; production unchanged),
and the 6R.205 production attach of the stored-haplotype loc-loop
object on that two-read hom-alt arm (INFO DP=2, MQ=47.00, SOR=2.303;
FORMAT/QUAL unchanged),
and the 6R.206 proof-only that remaining INFO DP at `2:92316416 C/A`
is Java loc-loop retainEvidence n=1 versus Rust region-wide stored
unique n=3 (empty one-read hom-alt annotation; production unchanged),
and the 6R.207 production attach of the stored-haplotype loc-loop
object on that one-read hom-alt arm (INFO DP=1, MQ=21.00, SOR=1.609;
FORMAT/QUAL unchanged),
and the 6R.208 proof-only that remaining INFO DP/SOR at
`2:92317399 C/A` is Java loc-loop retainEvidence n=2 versus Rust
attached n=1 after the 6R.180 same-QNAME collapse (FORMAT/QUAL match;
production unchanged),
and the 6R.209 production skip of that same-QNAME collapse on the
loc-loop retainEvidence object so both mates survive (INFO DP=2,
SOR=0.693; FORMAT/QUAL/MQ unchanged),
and the 6R.210 production attach of the stored-haplotype loc-loop
object on the gap-sparse shaped-early FORMAT path (INFO DP=1,
MQ=24.00, SOR=1.609; MQRankSum omitted; FORMAT/QUAL unchanged),
and the 6R.211 production attach of the same object on the gap-tail
het early-template FORMAT path (INFO DP=3, MQ=28.03, SOR=0.223;
FORMAT/QUAL unchanged),
and the 6R.212 production attach of the same object on the weak-sparse
het early-template FORMAT path (INFO DP=3, MQ=28.03, SOR=1.179;
FORMAT/QUAL unchanged),
and the 6R.213 proof that the earliest remaining common-site split is
FORMAT PL at `20:29455015 G/T` (Java `69,0,2140` vs Rust `122,0,2304`;
QUAL/GQ/QD follow those GLs; GT/AD/INFO match; production unchanged),
and the 6R.214 proof that those GLs already sit on a different
pre-marginalization haplotype matrix (Java 30 trimmed 130 bp haplotypes
vs Rust 84 of 181/212 bp; production unchanged),
and the 6R.215 proof that the extra haplotypes are materialized by
post-SeqGraph RT k=10 merge after Java-equivalent k=25 findBestPaths
(production unchanged),
and the 6R.216 proof that Java never k-bests that cyclic k=10
ReadThreadingGraph after SeqGraph findBestPaths (production
unchanged; global RT disable is not the Java rule),
and the 6R.217 diagnostic that createGraph cycle abort at RT extract
drops the 44 extras (production unchanged),
and the 6R.218 production SeqGraph-path extract abort (cyclic k=10
contributes 0; assemble 78; PL/GQ/QUAL/QD match Java at
`20:29455015 G/T`),
and the 6R.219 proof that FORMAT AD at `20:29455379 G/A` is not
Java `DepthPerAlleleBySample` on retainEvidence (Rust 44,5 vs Java
42,5; production unchanged),
and the 6R.220 proof that the 52-row remarg is isolated
`try_genotype` while `call_region` FORMAT is a later object
(production unchanged),
and the 6R.221 proof that `assign_genotype_likelihoods_for_region`
is that FORMAT object while production-arg `try_genotype` stays
the 52-row state (production unchanged),
and the 6R.222 proof that predecessor `20:29455375 T/A` REJECT
poisons TLS region-likelihood rows so target SiteScore emits
`AD 44,5` (production unchanged),
and the 6R.223 proof that the first inner operation is
`with_region_likelihood_rows` cache HIT (production unchanged),
and the 6R.224 proof that the HIT is allocator reuse of a dropped
2808-cell subset Vec (production unchanged),
and the 6R.225 proof that the cache must memoize logical
likelihood-row identity rather than a dropped Vec pointer
(production unchanged),
and the 6R.226 production change that keys the cache by exact
sparse-cell identity + `n_haps`,
and the 6R.227 proof that isolated and production loc-loop both
consume P2 so remaining vs Java is the 52-row remarg (production
unchanged),
and the 6R.228 proof that Java's event-local object is 47 reads
and the five extra Rust REF rows are absent from Java hap_ll
(production unchanged),
and the 6R.229 proof that the first Java drop of those five is
`filterPoorlyModeledEvidence` (production unchanged),
and the 6R.230 proof that those five already diverge at PairHMM
read-sequence input (Rust clip to `29455355` vs Java padded-window
bases; production unchanged),
and the 6R.231 proof that the first clipping mutation is
`hard_clip_to_region` on trim padded `29455355=29455375−20`
versus Java `29455294=29455314−20` (production unchanged),
and the 6R.232 proof that the missing `20:29455314 G>C` event is a
Java-only assembled haplotype, not an EventMap drop (production
unchanged),
and the 6R.233 proof that that haplotype is a SeqGraph path both
sides already have and that production k-best `K=128` is the first
miss (Java rank 116 vs Rust rank 129 for FNV `c7acc50dfb9f9ecc`;
production unchanged),
and the 6R.234 proof that those ranks come from different k-best
path states (Java 18 edges vs Rust 19; extra `log10(63/78)` sink
split), not from the `K=128` cap itself (production unchanged),
and the 6R.235 proof that the 30+8 vs 38 bp sink is SeqGraph
topology (pre-zip join at `TGTTTCTT`; zip emits the 8 bp sink;
production unchanged),
and the 6R.236 proof that the join is dangling-tail
`addEdge(971,428,weight=1)` onto the reference k-mer after prune
(absent at raw thread; Java has no equivalent splice; production
unchanged),
and the 6R.237 proof that Java `mergeDanglingTail` rejects the same
`3I9M` candidate (`refIndexToMerge=0`) while Rust `saturating_sub`
yields 1 (production unchanged),
and the 6R.238 proof that Java `0` is a path-index sentinel (LCA /
no merge), not a graph vertex id, while Rust `saturating_sub` maps
the underflow onto path index 1 (production unchanged),
and the 6R.239 production fix that checked
`(last_ref_idx + 1).checked_sub(matching_suffix)` preserves Java's
path-index-0 no-splice sentinel (carrier rank 116 at unchanged
`K=128`),
and the 6R.240 proof that covering `20:29455314 G>C` is in both
EventMaps and omitted by Java `stand-call-conf=30` while Rust
emitted at `stand_emit=10` (production unchanged),
and the 6R.241 production change that strict-Java emit uses
Java `standardConfidenceForCalling=30` so covering G>C is omitted
while GLs stay `PL 21,0,1461` QUAL ~13.63,
and the 6R.242 proof that the next genuine remaining VCF split is
common-site INFO DP at `20:29455649 T/TGTTTG` (Java 123 vs Rust 230;
production unchanged),
and the 6R.243 proof that empty annotation_likelihoods there is
colocated-merge construction after retainEvidence n=123 already exists
(production unchanged),
and the 6R.244 proof that merge drops the in-scope 123-read `subset`
instead of Java `prepareReadAlleleLikelihoodsForAnnotation` reuse
(production unchanged),
and the 6R.245 proof that the minimal Java-equivalent attach is
`subset.into_owned()` onto the Call (production unchanged),
and the 6R.246 production change that performs that attach so INFO DP
at `20:29455649 T/TGTTTG` is 123,
and the 6R.247 proof that remaining PL 3518 vs Java 3517 is calculator
1/1 log10 GL divergence, not rounding (production unchanged),
and the 6R.248 proof that the 1/1 GL is merged 2/2 (`TGTTTG/TGTTTG`)
and the homozygous calculator is Java-equivalent; the 3517.5 crossing
is already in the allele-row `L(read|TGTTTG)` column (production unchanged),
and the 6R.249 proof that that column is `max` over five EventMap
`T/TGTTTG` haplotypes matching Java `marginalize`; remaining inputs are
the five PairHMM columns (production unchanged),
and the 6R.250 proof that those columns share a Java-equivalent
read/quality plane while the known GKL-float residual cannot cross
3517.5 (production unchanged),
and the 6R.251 proof that the five PairHMM 161-mers are EventMap-exact
trim subsequences sharing insertion `GTTTG` (Java bytes unverified;
production unchanged),
and the 6R.252 proof that Java 4.4.0.0 assembleReads parents match
while trim spans differ (`29455560-29455728` vs `29455569-29455724`;
production unchanged),
and the 6R.253 proof that the trim split is Java STR padding 84 on
`A/AT` @ 29455644 versus Rust indel pad 75 (production unchanged),
and the 6R.254 measurement that Java 174-mers on frozen Rust 123-read
evidence move 2/2 *away* from 3517.5 (diagnostic PL 3519; production
unchanged),
and the 6R.255 measurement that Java-span hard-clip of those reads
against Rust 161-mers moves 2/2 further away (diagnostic PL 4280;
production unchanged),
and the 6R.256 measurement that Java 174-mers **and** Java clip together
still land at diagnostic PL 3519 (continuous 3519.157, above 3517.5;
production unchanged),
and the 6R.257 measurement that those joint primitive PairHMM arrays
match Java 4.4 `modifyReadQualities`+GCP (remaining difference is
kernel configuration, not trim/clip; production unchanged),
and the 6R.258 measurement that the first backend value is GKL float
`1/174` vs f64 `1/174` and is not PL-causal (production unchanged),
and the 6R.259 measurement that the first remaining ph2pr-chain value
is GKL float `q/10` vs f64 `q/10` and is not PL-causal (production
unchanged),
and the 6R.260 measurement that the first isolated `powf` split is Q=20
(exponent −2.0 exact) and is not PL-causal (production unchanged),
and the 6R.261 measurement that GKL float `1.f − ph2pr` at Q=20 is not
PL-causal (production unchanged),
and the 6R.262 measurement that GKL float `ph2pr / 3.f` at Q=20 is not
PL-causal (production unchanged),
and the 6R.263 measurement that GKL AVX distm blend is bit-preserving
and not PL-causal (production unchanged),
and the 6R.264 measurement that GKL AVX `M*distm` f32 multiply is not
PL-causal (production unchanged),
and the 6R.265 measurement that GKL AVX `X` first primitive
`M_t_1 * pMX` f32 multiply is not PL-causal (production unchanged),
not a new Yes row.
Historical L6–L14 narratives live on `pre-cleanup-archive` only — not unqualified **Yes** rows here.

## This is an AI-implementation stress test

The point of this repository is not “a GATK product in Rust.” It is a **deliberate experiment**: can agent-written code survive a mature, scientifically loaded Java codebase if a human owns the claim boundary?

**Agents produced the implementation.** **A human specified the experiment, redirected agents, rejected extras, and signed every claim and non-claim.** The evidence of that discipline is [`CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) — including what we do **not** claim. That combination (AI speed + human-owned gates) is the skill under test. It is not a substitute for a clinical drop-in, and it is not “the AI did it, so trust it.”

## Why does this exist?

A Rust reimplementation of the GATK4 HaplotypeCaller germline spine
(HC → CombineGVCFs → GenotypeGVCFs → hard-filter VariantFiltration).

On the **signed scopes** in [`CLAIM_MATRIX`](docs/CLAIM_MATRIX.md) there is
reproducible gate evidence — and that file is the only list of what this tree
does **not** claim. That is not genome-wide equivalence and not a clinical drop-in.

This is also an **AI-implementation stress test**: every line of implementation
in this repository was written by AI. I did not type the Rust. What I *did*
do, for months, was specify, review, redirect agents that were confidently
heading the wrong way, and rein things in when they got creative in ways
nobody asked for. The scarce skill is **owning the claim boundary while
agents write code** — not pretending the agents were not there, and not
pretending they do not need a human who can tell when they are wrong.

Authorship notes live in [`NOTICE.md`](NOTICE.md). Validity is the claim
matrix and the tests, not who typed the source.

If you know genomics and something here is wrong, open an issue.

### Live equivalence dashboard

**[Equivalence dashboard (GitHub Pages)](https://synapticfour.github.io/gatk-rs/)**
(source: [`docs/parity-site/`](docs/parity-site/)) — Chart.js view of hap.py
metrics **only after** a signed publish lands in `history.json`.
Until then the UI shows **“No published runs yet.”** Treat the site as
**instrumentation**, not GIAB evidence. Signed vs unsigned rows live only in
[`CLAIM_MATRIX`](docs/CLAIM_MATRIX.md).

## Validated Scope

gatk-rs targets the **germline short-variant** workflow against pinned GATK
**4.4** — a focused experiment, not a toolkit clone:

**HaplotypeCaller → CombineGVCFs → GenotypeGVCFs → VariantFiltration** (hard filters).

Claims and non-claims for what that path has actually proven:
[`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). Structure:
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md). Mid-B HC algorithm contracts:
[`docs/PARITY.md`](docs/PARITY.md).

**External pilots:** run the same spine on your own BAMs and compare Java↔Rust
yourself — [`docs/PILOT_GUIDE.md`](docs/PILOT_GUIDE.md) +
[`scripts/pilot/compare_callsets.py`](scripts/pilot/compare_callsets.py).

**Architecture decision (read this first for scope):**
[**ADR 0001 — Scope boundary**](docs/adr/0001-scope-boundary.md)
(why BQSR, VQSR, Mutect2, gCNV/SV, and Funcotator are intentionally out of
scope). Related: [ADR 0002](docs/adr/0002-remove-gatk-tools.md) (no generic
`gatk-tools` crate — use **samtools** / **bcftools** for sort/index/view).

## What this is

A native Rust workspace focused on **HaplotypeCaller** plus the post-call
helpers above. Generic BAM/VCF utilities belong in **samtools** / **bcftools**
([ADR 0002](docs/adr/0002-remove-gatk-tools.md)); leftover utility subcommands
stay callable for parity harnesses but are hidden from default `--help`.

### VariantFiltration vs VQSR

`gatk-rs variant-filtration` implements GATK-compatible **hard-filtering**
(`--filter-expression` / `--filter-name`, plus `--preset snp|indel` for the
official Best Practices tables). **This does not replace VQSR algorithmically** —
it is the recommended pragmatic fallback for smaller cohorts where VQSR is not
cleanly trainable. That matches official GATK guidance (VQSR is recommended only
once the cohort is large enough to support model training). See
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) and [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md).

## Build and run

```bash
cargo build -p gatk-cli --release
# `--threads N` / `--nt N` / `-t N` sizes the Rayon pool (default: CPU count).
# Active regions are processed in parallel inside one process; VCF rows are
# sorted/deduped before write (byte-identical across thread counts).
./target/release/gatk-rs HaplotypeCaller -R ref.fa -I reads.bam -O out.vcf \
  -L chr2:92300000-92350000 --threads 8
# Official SNP hard filters (GATK Best Practices table):
./target/release/gatk-rs variant-filtration -V snps.vcf -O snps.filtered.vcf --preset snp
# Or explicit Java-compatible pairs:
./target/release/gatk-rs variant-filtration -V snps.vcf -O snps.filtered.vcf \
  --filter-expression "QD < 2.0" --filter-name QD2 \
  --filter-expression "FS > 60.0" --filter-name FS60
./target/release/gatk-rs --help
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for tests and contribution workflow.

## Performance (fair comparison)

End-to-end HaplotypeCaller timings are published only from the dedicated quiet
host ([`docs/ci/PERF_BENCHMARK_HOST.md`](docs/ci/PERF_BENCHMARK_HOST.md)), never
from GitHub-hosted runners or the genomewide correctness VM.

| Contract | Detail |
|----------|--------|
| Workflow | [`.github/workflows/benchmark.yml`](.github/workflows/benchmark.yml) |
| Harness | [`scripts/perf/run_fair_hc_comparison.sh`](scripts/perf/run_fair_hc_comparison.sh) |
| Configs | gatk-rs `LOGLESS_HMM` (scalar), gatk-rs `AVX`/SIMD, Java **`FASTEST_AVAILABLE`**, Java `LOGLESS_CACHING` |
| Stats | ≥5 repeats → **median ± sample stdev** (not best-of-N) |
| Regions | small / medium / large nested GIAB windows |
| Metrics | wall, user, sys, Peak-RSS; RAPL energy when `perf` allows |
| Primary baseline | Java **`FASTEST_AVAILABLE`** (native AVX verified before timing) |
| Raw report | [`docs/perf/FAIR_HC_COMPARISON.md`](docs/perf/FAIR_HC_COMPARISON.md) · JSON · host [`HOST_SPECS.md`](docs/perf/HOST_SPECS.md) |
| Dashboard | Performance tab on the [public site](docs/parity-site/) (`data/perf_history.json`) |

**Production PairHMM default remains `LOG10_PAIRHMM`** until a signed GIAB run
clears SIMD/Logless. SIMD code path: `--pair-hmm AVX` (unit gate
`pairhmm_simd_vs_scalar_test`).

Headline HC speedups appear in `FAIR_HC_COMPARISON.md` after the first successful
dedicated-host run. Until then, treat any laptop Criterion numbers (e.g. local
NEON microbench in [`docs/perf/PAIRHMM_SPEEDUP.md`](docs/perf/PAIRHMM_SPEEDUP.md))
as **dev-host only** — not interchangeable with the fair suite.

## Memory profile (Peak-RSS)

Two labeled Peak-RSS profiles via
[`scripts/perf/run_hc_memory_profile.sh`](scripts/perf/run_hc_memory_profile.sh)
vs pinned Java GATK **4.4.0.0** (`-Xms1g -Xmx4g`). Full tables, commands, and
raw logs: [`docs/perf/HC_MEMORY_PROFILE.md`](docs/perf/HC_MEMORY_PROFILE.md).

### A. Trivial smoke — not a public memory claim

Checked-in fixture `parity/fixtures/`, interval `chr1:1-32` (32 bp). Peak-RSS
on this window is dominated by JVM/runtime fixed cost. **Do not quote a
“X% less memory” figure from it.** Raw table: [`docs/perf/HC_MEMORY_PROFILE.md`](docs/perf/HC_MEMORY_PROFILE.md).

### B. Realistic GIAB-dense window — public-claim basis

Multi-Mb NA12878 window on the known-dense chr20 locus
(default `20:10000000-12000000`, 2 Mb). **Only this profile** may back a public
memory claim, and **only** when measured on the dedicated
`gatk-rs-benchmark` host ([`docs/ci/PERF_BENCHMARK_HOST.md`](docs/ci/PERF_BENCHMARK_HOST.md))
with [`docs/perf/HOST_SPECS.md`](docs/perf/HOST_SPECS.md) populated.

| Engine | Peak RSS |
|--------|----------|
| gatk-rs / Java GATK 4.4 | *Pending dedicated-host run* — see `HC_MEMORY_PROFILE.md` |

Until that host run lands, do **not** advertise a genome-wide memory savings %.

## Equivalence

| Path | Purpose |
|------|---------|
| [`gatk-rs-equiv/`](gatk-rs-equiv/) | GIAB / hap.py / vcfeval + differential fuzz |
| [`scripts/parity/`](scripts/parity/) | L2 / P12 / GIAB harness scripts |
| [`tools/equivalence/README.md`](tools/equivalence/README.md) | Index of the proof surface |
| [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) | What those tests do and do not claim |

## Changelog

See [`CHANGELOG.md`](CHANGELOG.md).

## License

Apache License 2.0 — see [`LICENSE`](LICENSE). Trademark, GIAB, and oracle-jar notes: [`NOTICE.md`](NOTICE.md).

## Security and conduct

- Vulnerability reporting: [`SECURITY.md`](SECURITY.md)
- Code of conduct: [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md)

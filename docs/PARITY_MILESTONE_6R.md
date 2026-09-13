# 6R milestone — canonical mid-B HaplotypeCaller path

**Date:** 2026-08-30  
**Pinned Java:** GATK **4.4.0.0** SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`  
**Docker:** `broadinstitute/gatk:4.4.0.0`  
**Interval:** `2:92317000-92319000` (ActiveFull `2:92317262–92317491`)

This is a **short** public note. Forensic 6R.23–6R.41 reports are not in the public
tree.

## What was proven

On the canonical mid-B region, Rust matches GATK 4.4 through the full HC path that
produces the three oracle SNPs:

`92317399 C/A`, `92317407 T/C`, `92317412 G/C`

with FORMAT `GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0`, MLEAC=1, MLEAF=0.500, QUAL 78.32
(printed), QD 25.36 / 28.73 / 30.97.

**Canonical mid-B: CONVERGED. Whole-codebase GATK 4.4 parity: NOT ESTABLISHED.**

## What was fixed (classes, not a diary)

Java-faithful contracts, each with an observable test: reference-haplotype k-mer
uniqueness boundary; dangling-head mismatch-cap abort; `getBasesForPath` any-source
expansion; allele-keep (no unique-supporter collapse on default HC); trimmer
`maxEnd`; AF EM last Dirichlet update; MLEAC from EM not called GT; QualByDepth
`Random(47382911)` jitter.

Philosophy: Java contract → observable test → first divergence → smallest general
fix → regression. See [`PARITY.md`](PARITY.md).

## What remains unknown

Genome-wide / autosome HC equivalence; other samples; intervals where reservoir
`nextInt` interleaves with QD; remaining claim-matrix waivers (W-H1, W-H3, …).

## Later independent holdouts (not a Yes row)

Chr20_tiny genotype-boundary discovery through default USE_PLS GT/PL/GQ, the
first AD write, Class-A3 preservation of an assigned calculator genotype, and
the hom-ref emission boundary (6R.130–6R.159) is recorded in
[`PARITY.md`](PARITY.md), as is the next genuine emitted-VCF split at
`2:92316347` (6R.160), the extra after-assemble deletion hap (6R.161),
its dangling-merge object (6R.162), Java `mergeDanglingTail`=`addEdge`
vs Rust haplotype materialization (6R.163), the 6R.164 removal of
that Rust-only standalone haplotype so `2:92316347` FORMAT matches Java,
the 6R.165 proof that remaining INFO FS/MQ/SOR are annotation
pileup-vs-likelihoods membership splits, and the 6R.166 production
switch of FS/SOR onto Java `getContingencyTable` post-filter
likelihoods, the 6R.167 MQ switch onto Java `sampleEvidence`,
the 6R.168 RMS formula (`sqrt(sum(mq²)/n)` → 40.25),
and the 6R.169 proof that remaining INFO extras are
`ReadPosRankSum` / `InbreedingCoeff` (Java-default annotators
suppressed at this site, not missing Java values),
and the 6R.170 proof that ReadPosRankSum still uses `region.reads`
pileup versus Java informative likelihoods (REF=0/ALT=3),
and the 6R.171 proof that undefined RankSum is collapsed to `0.0` and
inserted (Java `emptyMap`), while a legitimate finite `0.0` must still
emit,
and the 6R.172 production change to `Option<f64>` so undefined omits
and finite `0.0` still emits (live pileup RankSum remains 6R.170),
and the 6R.173 proof-only live VCF reconnaissance that the next genuine
remaining field is INFO DP at `2:92305634 G/T` (Java 3 vs Rust 2),
and the 6R.174 production switch of INFO DP onto Java `Coverage.evidenceCount`
(FORMAT DP stays 2),
and the 6R.175 proof that remaining SOR at `2:92305634` is
`[0,0;1,2]` vs Java `[0,0;1,1]` because FLAG=145 (mate on contig 9)
is still PairHMM-informative (Java `filterNonPassingReads` +
`addEvidence(0)`),
and the 6R.176 production PairHMM mate-contig membership so that
read is `addEvidence(0)` only: SOR `[0,0;1,1]` → 0.693, INFO DP stays 3,
and the 6R.177 proof that InbreedingCoeff=1.0 at that site is
always-insert vs Java `MIN_SAMPLES=10` `emptyMap`,
and the 6R.178 production gate so `hc_info_values` omits
InbreedingCoeff when `n_genotypes < 10` (formula still `1 - het/n`),
and the 6R.179 proof that the next INFO split (`2:92307324 TTC/T`)
is region-wide vs per-variant annotation membership (DP/MQ/SOR),
and the 6R.180 production bind of INFO DP/MQ/SOR onto that
per-variant genotyping `AlleleLikelihoods` (n=1, MAPQ=44),
and the 6R.181 proof that `2:92307333 T/G` still uses region-wide
annotation evidence because the cluster-TG early-template path never
constructs a per-variant `AlleleLikelihoods`,
and the 6R.182 proof that Java’s annotation object at that site is a
**new** `marginalize` of the post-`filterPoorlyModeledEvidence` hap
matrix (n=2) then `retainEvidence` (n=1), reused for annotation; Rust
naive overlap is n=6 because the stored hap matrix is n=8,
and the 6R.183 proof that extras KEEP does **not** fire at that site:
the Java-equivalent poorly-modeled pass is already n=2, then
`refresh_region_read_likelihoods(apply_normalize=false)` restores n=8,
and the 6R.184 proof that the last of those refreshes is the final
membership write and that replaying the existing normalize+filter on
it restores Java n=2 (then diagnostic `retainEvidence` n=1),
and the 6R.185 production restore of that Java-order lifecycle after
the last P12 refresh (stored n=2; annotation object still later).
It is **not** a
claim-matrix Yes row and does not
establish chr20 VCF allele-set closure.
Production SeqGraph k-best remains `legacy_1024`.

## Test status (6R.41 / 6R.42 hygiene)

After this documentation cleanup: `gatk-haplotypecaller --lib` **576** passed
(`--test-threads=1`); `p12_call_none_mid_b_test` **1** passed;
`cargo fmt --all -- --check` ok. Dangling-recovery unit tests were moved out of
the production module so the N-3 size gate still holds; algorithm unchanged.

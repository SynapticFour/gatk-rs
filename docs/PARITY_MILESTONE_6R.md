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
the last P12 refresh (stored n=2),
and the 6R.186 construction of the missing per-variant annotation
`AlleleLikelihoods` on the cluster-TG early-template path
(stored n=2 → `marginalize` → `retainEvidence` n=1; INFO DP=1 MQ=44
SOR=1.609),
and the 6R.187 proof-only reconnaissance that the next remaining
INFO field is DP at `2:92305635 A/G` (Rust attached genotyping-subset
n=1 vs Java stored-hap `retainEvidence` n=3; not the 6R.186 cause),
and the 6R.188 proof that this n=1 object is the FORMAT-narrowed
SiteScore subset (6R.186 cluster-TG is not selected),
and the 6R.189 production bind of SiteScore annotation to the
stored-hap loc-loop object (`marginalize` then `retainEvidence`;
INFO DP=3 MQ=41.96 SOR=0.693 at `2:92305635 A/G`; FORMAT unchanged),
and the 6R.190 proof-only reconnaissance that the next remaining
INFO field is rust-only `ReadPosRankSum` at `2:92305716 A/C`
(Java `fillQualsFromLikelihood` empty-REF omit vs Rust
`region.reads` pileup emit; 6R.189 site stays closed),
and the 6R.191 production bind of ReadPosRankSum onto the same
per-variant annotation `AlleleLikelihoods` (target `2:92305716 A/C`
REF=0 ALT=3 → key omitted; formula/undefined semantics unchanged;
BaseQRankSum/MQRankSum still header-only),
and the 6R.192 proof-only reconnaissance that the next remaining
INFO field is Java-only `BaseQRankSum`/`MQRankSum` at
`2:92307359 CT/C` (`hc_info_values` never inserts; ReadPos already
matches; QD is QUAL/2),
and the 6R.193 production emit of BaseQRankSum from the same
fillQuals membership (target REF=1 ALT=1 → `Some(0.0)`; MQRankSum
untouched),
and the 6R.194 production emit of MQRankSum from the same
fillQuals membership (target REF=1 ALT=1 unequal MAPQ →
`Some(-0.674)`; QUAL/QD not retuned),
and the 6R.195 proof-only that QD at `2:92307359 CT/C` follows
QUAL/`getDepth` (Java 15.80 = Java QUAL/2; Rust 15.82 = Rust
QUAL/2; not an independent QD arrow; production unchanged),
and the 6R.196 production bind of biallelic QUAL AF to Java’s
length-based alt Dirichlet prior (indel QUAL 31.60; same-PL SNPs
stay 31.64; QD follows QUAL),
and the 6R.197 proof-only reconnaissance that the next remaining
INFO field is DP at `2:92307403 C/A` (Java 6 vs Rust 4; Java-only
BaseQRankSum/ReadPosRankSum; MQRankSum already matches; 6R.196 pin
stays QUAL 31.60),
and the 6R.198 proof-only that this site takes the cluster-downstream
early-template (no loc-loop annotation object; stored hap n=6,
retainEvidence n=4; Java INFO DP=6 is evidenceCount of the n=6
stored set; 6R.186 helper is the wrong 6-row object),
and the 6R.199 production attach of the stored-unique n=6 annotation
object on that cluster-downstream arm (not 6R.186 retain n=4;
FORMAT/QUAL unchanged; INFO DP still 4 is a Coverage overlap
6R.200 arrow),
and the 6R.200 proof-only that this remaining DP split is a second
Coverage overlap filter (`java_alignment_read_overlaps_interval`)
on the n=6 object; Java `Coverage.annotate` is `evidenceCount()`
with no overlap; production unchanged,
and the 6R.201 production removal of that second Coverage overlap
so INFO DP consumes attached unique evidenceCount (n=6) at
`2:92307403 C/A` (FORMAT/QUAL/MQ unchanged; RankSums not retuned),
and the 6R.202 proof-only that remaining Java-only BaseQRankSum /
ReadPosRankSum is empty-REF fillQuals `getElementForRead` on the two
extra stored unique rows (RankSums are invoked; MQRankSum already
1.834; production unchanged),
and the 6R.203 production shared RankSum `getElementForRead`
covering-CIGAR retry so those REF rows yield Java-compatible
BaseQ/ReadPos elements at `2:92307403 C/A`
(`BaseQRankSum=-1.834`, `ReadPosRankSum=1.282`; MQRankSum stays
1.834; INFO DP stays 6; FORMAT/QUAL unchanged),
and the 6R.204 proof-only that remaining INFO DP at `2:92316296 A/T`
is Java loc-loop retainEvidence n=2 versus Rust region-wide stored
unique n=3 (empty two-read hom-alt annotation; production unchanged),
and the 6R.205 production attach of
`annotation_likelihoods_from_stored_haplotypes` on that two-read
hom-alt arm so the annotation object is retainEvidence n=2 (INFO DP=2,
MQ=47.00, SOR=2.303; FORMAT/QUAL unchanged),
and the 6R.206 proof-only that remaining INFO DP at `2:92316416 C/A`
is Java loc-loop retainEvidence n=1 versus Rust region-wide stored
unique n=3 (empty one-read hom-alt annotation; production unchanged),
and the 6R.207 production attach of
`annotation_likelihoods_from_stored_haplotypes` on that one-read
hom-alt arm so the annotation object is retainEvidence n=1 (INFO DP=1,
MQ=21.00, SOR=1.609; FORMAT/QUAL unchanged),
and the 6R.208 proof-only that remaining INFO DP/SOR at
`2:92317399 C/A` is Java loc-loop retainEvidence n=2 versus Rust
attached n=1 after the 6R.180 same-QNAME collapse (FORMAT/QUAL match;
production unchanged),
and the 6R.209 production skip of that same-QNAME collapse on the
loc-loop retainEvidence object so both mates survive (INFO DP=2,
SOR=0.693; FORMAT/QUAL/MQ unchanged; 6R.180/6R.186 collapse stays
on FORMAT-subset / stored-unique callers),
and the 6R.210 production attach of
`annotation_likelihoods_from_stored_haplotypes` on the gap-sparse
shaped-early FORMAT path so the annotation object is retainEvidence
n=1 at `2:92318199 C/T` (INFO DP=1, MQ=24.00, SOR=1.609; MQRankSum
omitted; FORMAT/QUAL unchanged),
and the 6R.211 production attach of the same helper on the gap-tail
het early-template FORMAT path so the annotation object is
retainEvidence n=3 at `2:92325193 C/T` (INFO DP=3, MQ=28.03,
SOR=0.223; FORMAT/QUAL unchanged),
and the 6R.212 production attach of the same helper on the weak-sparse
het early-template FORMAT path so the annotation object is
retainEvidence n=3 at `2:92325268 C/T` (INFO DP=3, MQ=28.03,
SOR=1.179; FORMAT/QUAL unchanged),
and the 6R.213 proof that the earliest remaining common-site split is
FORMAT PL at `20:29455015 G/T` (`69,0,2140` vs `122,0,2304`; QUAL/GQ/QD
downstream; GT/AD/INFO match; production unchanged),
and the 6R.214 proof that the first divergent object behind those GLs is
the pre-marginalization haplotype matrix (Java 30 trimmed haplotypes vs
Rust 84; 13 extra window sequences; production unchanged),
and the 6R.215 proof that those extra haplotypes are inserted by
`merge_rt_kbest_pre_remove_paths` (RT k=10) after SeqGraph k=25 already
matches Java's 78 (production unchanged),
and the 6R.216 proof that Java's default SeqGraph assembler never
k-bests the cyclic k=10 ReadThreadingGraph those 44 come from
(production unchanged; do not disable RT merge globally),
and the 6R.217 diagnostic that honoring `createGraph` cycle abort at
RT extract restores unique merge=78 without suppressing p11/indel4
acyclic k=10 ALTs (production unchanged),
and the 6R.218 production SeqGraph-path extract abort so cyclic k=10
contributes 0, assemble=78, hap_n=30, and PL/GQ/QUAL/QD match Java
at `20:29455015 G/T`,
and the 6R.219 proof that FORMAT AD at `20:29455379 G/A` is a
FORMAT-specific subset (44,5) rather than Java
`DepthPerAlleleBySample` on retainEvidence (production unchanged),
and the 6R.220 proof that the 52-row remarg is isolated
`try_genotype` (`AD 47,5`) while `call_region` FORMAT is a later
object (`AD 44,5`) with no Java post-retainEvidence subset
(production unchanged),
and the 6R.221 proof that `assign_genotype_likelihoods_for_region`
is that FORMAT object while production-arg `try_genotype` stays
`AD 47,5` (production unchanged),
and the 6R.222 proof that predecessor `20:29455375 T/A` REJECT
poisons TLS region-likelihood rows so the target SiteScore emits
`AD 44,5` / `PL 78,0,1811` (production unchanged),
and the 6R.223 proof that the first inner operation is
`with_region_likelihood_rows` cache HIT (production unchanged),
and the 6R.224 proof that the HIT is allocator reuse of a dropped
2808-cell subset Vec so `key(P1)==key(P2)` while P1≠P2 (production
unchanged),
and the 6R.225 proof that the TLS row cache keys storage identity
rather than sparse-cell identity, so the smallest safe contract is
HIT iff the same logical likelihood-row population (production
unchanged),
and the 6R.226 production change that keys the TLS row cache by
exact sparse-cell identity + `n_haps` so `20:29455379 G/A` rebuilds
P2 rather than consuming predecessor dense rows,
and the 6R.227 proof that isolated and production both consume that
P2 52×2, so remaining vs Java is the 52-row remarg versus Java's
47-read object (production unchanged),
and the 6R.228 proof that Java's live post-retainEvidence set is
47 reads and the five extra Rust REF rows are missing from Java
hap_ll (production unchanged),
and the 6R.229 proof that those five enter Java PairHMM / hap_ll
and DROP together at `filterPoorlyModeledEvidence` (production
unchanged),
and the 6R.230 proof that the first likelihood-value split is
PairHMM read-sequence input (Rust clip to `29455355` vs Java
padded-window bases; production unchanged),
and the 6R.231 proof that the first clipping mutation is
`clip_finalized_reads_in_place` / `hard_clip_to_region` on the
trimming padded span `29455355=29455375−20` versus Java
`29455294=29455314−20` (production unchanged),
and the 6R.232 proof that Java `20:29455314 G>C` exists because
assembly materializes haplotype `c7acc50dfb9f9ecc` while Rust's 128
haplotypes are all `G` there (production unchanged),
and the 6R.233 proof that the C-bearing SeqGraph path already
exists in both (edge multiplicity 3) and that production k-best
`K=128` is the first miss (Java rank 116 vs Rust rank 129 for the
same FNV; production unchanged),
and the 6R.234 proof that the rank shift is a k-best path-state
split (18 vs 19 edges; extra Rust `log10(63/78)`), so `K=128` is
a downstream cutoff (production unchanged),
and the 6R.235 proof that the 30+8 vs 38 bp sink is SeqGraph
topology (pre-zip join at `TGTTTCTT`; initial zip emits the 8 bp
sink; production unchanged),
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
and the 6R.240 proof that covering `20:29455314 G>C` is present
in both EventMaps and omitted by Java `calculateGenotypes` at
`stand-call-conf=30` while Rust emitted at `stand_emit=10`
(production unchanged),
and the 6R.241 production change that
`HcGenotypingConfig::strict_java` uses Java
`standardConfidenceForCalling=30` so covering G>C is omitted
(`PL 21,0,1461` QUAL ~13.63 unchanged; non-strict emit stays 10),
and the 6R.242 proof that the next genuine remaining VCF split is
common-site INFO DP at `20:29455649 T/TGTTTG` (Java 123 vs Rust 230;
covering G>C omitted by both; production unchanged),
and the 6R.243 proof that empty annotation_likelihoods at that site is
colocated-merge construction (`Vec::new()`) after retainEvidence n=123
already exists (production unchanged),
and the 6R.244 proof that merge drops the in-scope 123-read `subset`
instead of Java annotation-object reuse (production unchanged),
and the 6R.245 proof that the minimal attach is `subset.into_owned()`
onto the Call (production unchanged),
and the 6R.246 production change that performs that attach so INFO DP
at `20:29455649 T/TGTTTG` is 123,
and the 6R.247 proof that remaining PL 3518 vs Java 3517 is the
calculator 1/1 log10 GL (continuous 3517.516…), not PL rounding
(production unchanged),
and the 6R.248 proof that merged 2/2 is `TGTTTG/TGTTTG`, the homozygous
calculator equals `Σ L(read|allele2)`, and the split is the allele-row
column (production unchanged),
and the 6R.249 proof that the TGTTTG column is `max` over five EventMap
haplotypes and that Java `marginalize` is the same max (production
unchanged),
and the 6R.250 proof that the five PairHMM columns share a
Java-equivalent read/quality plane and that the known GKL-float residual
cannot cross 3517.5 (production unchanged),
and the 6R.251 proof that the five 161-mers are EventMap-exact trim
subsequences sharing insertion `GTTTG` (Java bytes unverified;
production unchanged),
and the 6R.252 proof that recovered GATK 4.4.0.0 has identical
untrimmed TGTTTG parents but a different trim span
(`29455560-29455728` vs `29455569-29455724`; production unchanged),
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
and the 6R.263 measurement that GKL AVX distm blend
(`_mm256_blendv_ps`) is bit-preserving and not PL-causal (production
unchanged),
and the 6R.264 measurement that GKL AVX `M*distm` `_mm256_mul_ps` is
not PL-causal (production unchanged),
and the 6R.265 measurement that GKL AVX `X` first primitive
`VEC_MUL(M_t_1, pMX)` is not PL-causal (production unchanged).
It is **not** a
claim-matrix Yes row and does not
establish chr20 VCF allele-set closure.
Production SeqGraph k-best remains `legacy_1024`.

## Test status (6R.41 / 6R.42 hygiene)

After this documentation cleanup: `gatk-haplotypecaller --lib` **576** passed
(`--test-threads=1`); `p12_call_none_mid_b_test` **1** passed;
`cargo fmt --all -- --check` ok. Dangling-recovery unit tests were moved out of
the production module so the N-3 size gate still holds; algorithm unchanged.

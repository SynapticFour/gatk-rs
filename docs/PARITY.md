# HaplotypeCaller parity (GATK 4.4)

Pinned Java: GATK **4.4.0.0**, SHA [`2dbc025821bc5f686c423ff332a41e6cef892a77`](https://github.com/broadinstitute/gatk/commit/2dbc025821bc5f686c423ff332a41e6cef892a77).  
See also [`GATK_PINNED.env`](GATK_PINNED.env) and root `GATK_PINNED_SHA`.

Product claims still live only in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). This page is the
**engineering** record of how HaplotypeCaller algorithm contracts are proven.

## Current achievement

The implementation has undergone source-backed GATK 4.4 parity investigation on the
canonical **mid-B** ActiveFull region (`2:92317262–92317491`, interval
`2:92317000-92319000`, NA12878 20k b37).

That path is **converged** from assembly through VCF INFO/FORMAT:

| Stage | Mid-B status |
|-------|----------------|
| Read-threading graph (k-mer selection, uniqueness on the extended-region haplotype) | Converged |
| Dangling recovery / `best_prefix_match` mismatch-cap | Converged |
| `path_bases` / `getBasesForPath` source expansion | Converged |
| `removePathsNotConnectedToRef` | Converged |
| SeqGraph zip / simplification | Converged |
| k-best haplotypes | Converged |
| EventMap | Converged |
| Haplotype trimming (`maxEnd`) | Converged |
| PairHMM / FORMAT GT, AD, DP, GQ, PL | Converged |
| AF calculator QUAL | Converged (printed QUAL 78.32) |
| MLEAC / MLEAF | Converged |
| QualByDepth / QD (GATK `Random` seed `47382911`) | Converged |

Oracle sites on that region:

```
2:92317399 C>A
2:92317407 T>C
2:92317412 G>C
```

with Java-equivalent FORMAT `GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0`, MLEAC=1, MLEAF=0.500,
QUAL 78.32 (Rust 78.323 before VCF 2-decimal print), and QD 25.36 / 28.73 / 30.97.

**Canonical mid-B parity: CONVERGED.**

**Whole-codebase GATK 4.4 HaplotypeCaller parity: NOT YET ESTABLISHED.**

## Frozen 6R canonical checkpoint

One witness is closed. It is not GATK-wide parity, not every region, not every
HaplotypeCaller mode, and not every GATK version.

Witness: `20:29455649 T/TGTTTG`, GATK 4.4.0.0, Java source SHA
`2dbc025821bc5f686c423ff332a41e6cef892a77`.

Rust matches that Java genotype-likelihood vector and the continuous PL vector
bit for bit. Emitted integer PL is `570,148,3485,0,2762,3517`, and the
biallelic subset is `570,0,3517`.

Production semantics on this path:

- anchored, alternate-aware tandem-repeat trim padding (`29455560–29455728` for `20:29455644 A>AT`);
- failed-mate reads are absent from PairHMM and genotyping evidence;
- clipped reads use Java `ReadCoordinateComparator`;
- heterozygote genotype combination uses the Java-compatible Jacobian helper.

The 122×3 allele matrix is bit-identical to Java, and the 122 genotyping reads
are an identity permutation of Java's order. Forcing the pre-6R.312
`(tid, pos, qname)` sort leaves that matrix unchanged and recreates the
historical GL gaps of 2 ULP and 1 ULP. Java order removes those gaps. No
tolerance or rounding workaround is used.

Cumulative 6R gates through 6R.312: 182. `HOLDOUT_6R243` was not run.
Details: [`parity/6R.312_FINAL_CLOSURE.md`](parity/6R.312_FINAL_CLOSURE.md).

6R.314 holds out the four production semantics on synthetic inputs that are
not `20:29455649 T/TGTTTG` and not `20:29455644 A>AT`. All five executed
holdouts passed. Cumulative 6R gates through 6R.314: 187. `HOLDOUT_6R243`
was not run. Details: [`parity/6R.314_HOLDOUT.md`](parity/6R.314_HOLDOUT.md).

This is not genome-wide equivalence, not a clinical drop-in, and not a claim that every
interval, sample, or annotation matches Java. Signed product scopes remain P12 / L2 /
synthetic joint gates in the claim matrix.

## Independent holdouts (6R.43–6R.50)

A frozen 10-region panel (`scripts/parity/6r43_holdout_panel.json`) is discovery, not a
genome-wide score. After 6R.45–6R.50:

| Class | Regions |
|-------|---------|
| **A** (7) | `ctrl_mid_b`, `p12_het_tail`, `p12_desert`, `p12_snp_cluster`, `p12_indel_mix`, `p12_post`, `p12_mid_a` |
| **C** (3) | `chr20_tiny`, `chr20_w47`, `chr21_w10` |

Previously green holdouts remain green. Canonical mid-B remains CONVERGED.

Remaining **C** regions are predominantly Stage E haplotype-content differences (Java
internal graph topology UNKNOWN) and Stage D/G/H cases where an allele exists in EventMap
but is not emitted.

## Independent chr20_tiny genotype-boundary holdouts (6R.130–6R.312)

Engineering discovery on `20:29455000-29456500` (target SNP `20:29456196 A/T`, HG001).
**Not** a claim-matrix Yes row and **not** chr20 VCF allele-set closure.

Diagnostic k-best uses `GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic` so the
Java top-128 haplotype population can be compared. Production SeqGraph k-best remains
`legacy_1024`. Resource policy is experimental only:
[`parity/KBEST_RESOURCE_POLICY.md`](parity/KBEST_RESOURCE_POLICY.md). Forensic
`docs/parity/*_REPORT.md` files stay local (gitignored).

Closed on this fixture, within the established Java-float / Rust-f64 residual:

haplotype list / trim / assemble → SeqGraph k-best (Java top-128 under diagnostic
policy) → EventMap genomic identity → PairHMM inputs / normalization → poorly-modeled
KEEP **201×25** → `filterAlleles` / realign (two best-haplotype winner-index differences
are non-causal) → haplotype→allele mapping **25/25** → max-marginalize → raw diploid GLs.

6R.150: `retainEvidence` interval is `SimpleInterval(mergedVC).expandWithinContig(2)`
then `target.overlaps(read)` on alignment coordinates (`20:29456194-29456198`).
KEEP 201 → overlap **41**. The alignment predicate matches.

6R.151: production Java-strict genotyping no longer runs
`dedupe_likelihood_subset_by_qname` on that overlap set. Java keeps overlapping
mates independently (`AlleleLikelihoods.retainEvidence`; Mutect `groupEvidence`
is a different API). Helper still exists for Rust-specific sparse/P12 augment
and empty-subset pileup rescue. Membership **41/41**. The previously dropped
mate `HWI-D00360:7:H88WKADXX:2:2107:6787:30989` flags=99 is retained.

`calculateGLsForThisEvent` is NO_CALL + PL. Default `USE_PLS_TO_ASSIGN` does not
apply `assumingHW` priors to GT/PL.

6R.152: after the equivalent 41-read object, Java `assumingHW` SNP log10 priors
(`0`, `log10(het)−log10(3)`, `2·log10(het)−log10(3)`) and Rust remainder-mass
simplex (`log10(1−het−hom_var)`, `log10(het)`, `log10(hom_var)`) are **not** a
common scale. Pairwise relative preferences diverge (0/1 vs 1/1 is 3 log10 vs
`log10(2)`). Ranking of the prior vector itself is still `0/0 > 0/1 > 1/1` on
both sides. On this fixture, USE_PLS GT/PL remain `0/0` and `0,5,1174`; swapping
priors does not change the winner. `production_change: NONE`.

6R.153: default `USE_PLS_TO_ASSIGN` converts log10 GLs to integer PL with
`Math.round(-10*(GL−max))`, selects GT as first `maxElementIndex` of the GLs
(Rust production: first min PL; unique winner agrees), and GQ as
`round(10*(best−second))` / PL gap. This fixture: GT `0/0`, PL `0,5,1174`,
GQ `5` both sides. PairHMM residual does not change integer PL.
`classification: NO_DIVERGENCE`. `production_change: NONE`.

6R.154: first FORMAT/AD write is Java `DepthPerAlleleBySample.annotateWithLikelihoods`
(`searchBestAllele` + `isInformative` with `confidence > 0.2`) on the 41-read
remaining A,T object. Rust `InformativeAd::from_marginalized_rows` matches
**row-by-row**. Live AD **37,4** (not the 6R.101 `36,19` from a different site).
Both restored mates contribute A. `classification: NO_DIVERGENCE`.
`production_change: NONE`. QUAL, L9, and VCF emission are **not** claimed.

6R.155: the first operation after that 37,4 write that mutates production
`genotyped_calls` is Rust-only `SiteReshape::apply_class_a_family` Class-A3
(`sparse_snp_genotype_from_read_depths` on QNAME-deduped `region.reads`
pileup **22,24**, not the 41-row likelihood object). Java has no analogue
after `annotateWithLikelihoods` (reverse-trim no-op; physical phasing off).
Live result: GT `0/1`, PL `81,0,36`, AD `22,24`, GQ 36, DP 46 versus Java
annotation `0/0`, `0,5,1174`, `37,4`. Diagnostic skip restores the
Java-equivalent genotype object; the site then fails the hom-ref emit gate.
`classification: RUST_ONLY_MULTI_FIELD_MUTATION`. `production_change: NONE`
(do not globally delete sparse infrastructure; Class-A3 GT-flip on default
HC is the candidate later lossless bypass). QUAL, L9, and VCF stay closed.

6R.156: Class-A3 is a Rust-only site heuristic after a **valid** PairHMM
genotype. Predicate: SNP, pileup both alleles and neither >2× the other,
informative REF ≥ 3× ALT. It does **not** test existing GT. Because
PL-argmin is 0/0, the sparse arm (`sparse_snp_genotype_from_read_depths`)
replaces GT/PL/AD/GQ/DP from QNAME-first-wins untrimmed `region.reads`
(22,24; 276 QNAMEs claim a slot without a locus base). Comments describe
GIAB-het AD restore / hom-alt GT-flip / “no alt-hap” rescue; the
implementation overrides Java-equivalent 0/0. Single caller
(`apply_class_a_family`). Other `sparse_snp_*` callers (L9, empty-mapper)
are separate. Counterfactual C/D restore 0/0 0,5,1174 37,4 without
deleting sparse infrastructure. `classification: A3_RUST_ONLY_SEMANTIC_OVERRIDE`.
`production_change: NONE`.

6R.157: `pl_gt == 0` is the **min-PL genotype index**, not a validity bit.
Valid 0/1 has `pl_gt==1`; valid 1/1 has `pl_gt==2`; every assigned PL vector
has some slot 0 after normalization. Class-A3 still does not consult the
existing GT: A3 predicates plus `pl_gt != 1` sparse-replace a valid 0/0 or
1/1 (the 1/1 row is overdetermined with Class-A2). Valid 0/1 keeps GT/PL and
still rewrites AD/DP. `RegionGenotypeResult` has no NO_CALL field at this
caller (`SiteScore` always assigns). Java 4.4 `USE_PLS_TO_ASSIGN` never
pileup-replaces an assigned 0/0, 0/1, or 1/1; `call == null` drops the site.
Guard A is not Java-equivalent. Guard B is: do not mutate when the calculator
already assigned a genotype. Coordinate-free forensic
`forensic_6r157_class_a3_calculator_boundary_contract`. Next-round candidate
(not implemented): in `apply_class_a_family` only, if `class_a3`, return the
incoming genotype. `classification: CANDIDATE_B_JAVA_EQUIVALENT_AT_A3_CALLER`.
`production_change: NONE`.

6R.158: `SiteReshape::apply_class_a_family` returns the incoming calculator
genotype when Class-A3 holds. Valid 0/0, 0/1, and 1/1 FORMAT fields are
unchanged (not a `pl_gt == 0` guard). The generic
`sparse_snp_genotype_from_read_depths` helper and L9 / empty-mapper /
cluster callers are untouched. Canonical site: A3 no longer writes
`0/1` `81,0,36` `22,24`; calculator `0/0` `0,5,1174` `37,4` GQ 5 DP 41
is preserved through reshape. `genotyped_calls` is then **absent**
(emit eligibility not chased in 6R.158).
`classification: A3_CALCULATOR_GENOTYPE_PRESERVATION`.
`production_change: YES`.

6R.159: after the preserved calculator 0/0 `PL=0,5,1174` `AD=37,4` `GQ=5`,
the first Rust emission predicate is `passes_java_emit_not_hom_ref`
(best diploid index `== 0`) inside `java_emit_would_pass`. GQ is not read
on that path (0/0 with GQ=40 still drops; 1/1 with GQ=5 still emits).
Java 4.4 default HC (`EMIT_VARIANTS_ONLY`, `stand-call-conf=30`, ERC off)
uses AF `siteIsMonomorphic` as `bestGuessIsRef`, not sample GT; this
vector is monomorphic at both 10 and 30 (QUAL 32.12, still dropped because
hom-ref under `EMIT_VARIANTS_ONLY`). Both sides omit the VCF record.
`genotyped_calls` is the `returnCalls` / `call != null` list, not all
calculator results. `classification: NO_VCF_EMISSION_DIVERGENCE`.
`production_change: NONE`.

6R.160: first genuine emitted-VCF difference after 6R.158/6R.159, in genomic
order on the runnable 6R.43 panel, is `2:92316347 G/A` (`p12_mid_a`).
Java `1/1 AD=0,3 PL=135,9,0`; Rust `1/1 AD=0,1 PL=45,3,0`. First proven
arrow: Java `DepthPerAlleleBySample` AD `0,3` ≠ Rust informative AD `2,5`
(PairHMM GLs het-best) after remarg; `retainEvidence` does not change that
AD; the sparse `45,3,0` overwrite is downstream. AF 10 vs 30 is not causal.
Canonical `20:29456196` stays closed on the diagnostic emit path (both
absent); production `run_haplotype_caller` still emits a later rust-only
het there. `classification: TRUE_CALLSET_DIFFERENCE`.
`production_change: NONE`.

6R.161: at `2:92316347 G/A`, Java AD `0,3` vs Rust pre-overwrite informative AD
`2,5` is a vote-membership consequence of a haplotype-population split.
Java `after_assemble` n=2 (REF `9b4092f3b20a5de7`, ALT `8c9a6fc8f302f4a3`);
Rust n=3 (same two plus deletion hap `3a53b2a941bbdc43`, cigar
`40D2M155D172M107D`). Untrimmed shared haplotypes match; trimmed PairHMM
inputs already differ (203 bp vs 395 bp). PairHMM float residual and the
0.2 informative cut are not causal. The later Rust `45,3,0` overwrite is
not this round. `classification: HAPLOTYPE_POPULATION_DIVERGENCE`.
`production_change: NONE`.

6R.162: extra after_assemble hap `3a53b2a941bbdc43` (`40D2M155D172M107D`,
len=174) is a dangling-merge haplotype object (`score=25`, `kmer_size=0`)
from `apply_dangling_merge_haplotypes` after RT extract at expanded k=35.
SeqGraph k-best at k=35 is n=2 and matches Java `assembleReads` (REF
`9b4092f3b20a5de7`, ALT `8c9a6fc8f302f4a3`). Configured k=10/25 skip
non-unique ref on both sides. k-best / `legacy_1024` is not causal.
`classification: HAPLOTYPE_CONSTRUCTION_DIVERGENCE`.
`production_change: NONE`.

6R.163: Java `mergeDanglingTail` only `addEdge(..., weight=1)`. Rust also
pushes `DanglingMergeHaplotype` (`path_bases` of the dangling alt walk)
and `apply_dangling_merge_haplotypes` inserts it (`score=25`). Extra
`3a53b2a941bbdc43` is a proper substring of SeqGraph ALT
`8c9a6fc8f302f4a3` at offset 195; not source→sink; `extra_in_kbest=false`.
`classification: HAPLOTYPE_MATERIALIZATION_DIVERGENCE`.
`production_change: NONE`.

6R.164: Java-exact dangling-tail recovery keeps the graph splice and no
longer materializes dangling `path_bases` as a standalone haplotype.
Target assemble n=2/2 (hashes `8c9a6fc8f302f4a3` / `9b4092f3b20a5de7`);
extra `3a53b2a941bbdc43` removed. FORMAT at `2:92316347 G/A` matches
Java `GT=1/1 AD=0,3 PL=135,9,0` QUAL 121.84. Remaining INFO FS/MQ/SOR
is not this round. `classification: HAPLOTYPE_MATERIALIZATION_DIVERGENCE`.
`production_change: YES`.

6R.165: proof-only. At `2:92316347 G/A` FORMAT/QUAL stay closed. INFO
FS/MQ/SOR diverge because Rust annotates from `region.reads` pileup
(`FS`/`SOR` table `[0,2;3,2]`, MQ = mean of five ALT MAPQs = 36.4)
while Java uses post-filter `AlleleLikelihoods` (informative strand
table `[0,0;2,1]` → FS=0 SOR=1.179; RMS of MAPQs 47,47,21 → MQ=40.25).
FS and SOR share that 2×2 membership split; MQ is a different
membership list. `production_change: NONE`.

6R.166: FS/SOR consume Java `StrandBiasTest.getContingencyTable` on
post-filter allele likelihoods (informative best allele + strand).
Target table `[0,0;2,1]` → FS=0 SOR=1.17865 (prints 1.179). FORMAT/QUAL
unchanged. MQ still the 6R.165 pileup mean 36.4.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: YES`.

6R.167: MQ consumes Java `RMSMappingQuality.calculateRawData` sampleEvidence
(unique post-filter likelihood-matrix reads, `MQ != 255`). Target MAPQs
`[47,47,21]` match Java. Production aggregation remains the arithmetic mean
(38.333…); Java RMS 40.25 is 6R.168. FORMAT/FS/SOR unchanged.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: YES`.

6R.168: MQ aggregation is Java `makeFinalizedAnnotationString`
(`sqrt(sum(mq²)/n)` then `%.2f`) on the frozen `[47,47,21]` list.
Target MQ=40.25. FORMAT/FS/SOR unchanged.
`classification: ANNOTATION_FORMULA_DIVERGENCE`.
`production_change: YES`.

6R.169: proof-only live INFO inventory at `2:92316347 G/A`. FORMAT/QUAL/FS/SOR/MQ
stay closed. Java record keys `AC,AF,AN,DP,ExcessHet,FS,MLEAC,MLEAF,MQ,QD,SOR`
match Rust at Java print precision. Rust-only record keys are
`ReadPosRankSum=0` and `InbreedingCoeff=1`. Both annotators are Java
`StandardAnnotation` but suppressed here (`RankSumTest` empty/NaN
informative REF; `InbreedingCoeff.MIN_SAMPLES=10`). First remaining
INFO arrow is ReadPosRankSum pileup vs informative likelihoods.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: NONE`.

6R.170: proof-only. ReadPosRankSum still uses `region.reads` pileup
(REF=2 ALT=5) while Java `RankSumTest.fillQualsFromLikelihood` on the
reused 6R.166 post-filter likelihoods is REF=0 ALT=3 → NaN → no INFO key.
Empty/NaN→insert `0` is a stacked later predicate, not this arrow.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: NONE`.

6R.171: proof-only. At the RankSum boundary, Java `RankSumTest` maps
NaN/empty → `emptyMap` (no INFO key) while a legitimate finite `0.0` is
still emitted (`%.3f`). Rust `rank_sum_z` returns `f64` and maps empty
either list (and MWU NaN) to `0.0` before `hc_info_values` always inserts
`ReadPosRankSum`. Case B identical lists `[10,20]` vs `[10,20]` is a real
`0.0`; omit-on-zero would be wrong. First downstream arrow is
`empty OR → return 0.0`.
`classification: ANNOTATION_EMIT_PREDICATE_DIVERGENCE`.
`production_change: NONE`.

6R.172: production. `read_pos_rank_sum` / `rank_sum_z` return `Option<f64>`:
empty or NaN → `None` (no INFO key); finite z including `Some(0.0)` → insert.
`hc_info_values` inserts ReadPosRankSum only for `Some(z)` (not a generic
skip-zero; FS=0 still emits). Evidence source remains 6R.170 pileup, so the
live target can still show a finite pileup RankSum. Java-equivalent REF=0
ALT=3 is omitted.
`classification: ANNOTATION_EMIT_PREDICATE_DIVERGENCE`.
`production_change: YES`.

6R.173: proof-only reconnaissance of the post-6R.172 live VCF vs frozen
6R.43 Java oracles. Record set 126/129 (Java-only 0, Rust-only 3 on
`chr20_tiny`). Chr2 FORMAT/QUAL closed. Known extras remain 6R.170
ReadPosRankSum at `2:92316347` and 6R.169 InbreedingCoeff (one sample).
First genuine remaining field is INFO DP at `2:92305634 G/T` (Java 3 vs
Rust 2); FORMAT DP already matches (2). Production unchanged.
`classification: RECONNAISSANCE_ONLY`.
`production_change: NONE`.

6R.174: production. INFO DP is Java `Coverage.annotate` →
`likelihoods.evidenceCount()` on post-filter retainEvidence (unique overlapping
reads), not FORMAT/`DepthPerSampleHC`. At `2:92305634 G/T` Java INFO DP=3 vs
previous Rust FORMAT copy=2; FORMAT DP stays 2. Third read
`H06JUADXX130110:1:1101:10018:4569 FLAG=145`. SOR left stacked.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE` (wrong source object:
Coverage vs FORMAT DP).
`production_change: YES`.

6R.175: proof-only. At `2:92305634 G/T` FORMAT/QUAL/INFO DP stay closed.
SOR 0.693 vs 1.179 is table membership, not `calculateSOR`. Java
`getContingencyTable` is `[0,0;1,1]` (the two 53752 mates). Rust 6R.166
helper is `[0,0;1,2]` because `H06JUADXX130110:1:1101:10018:4569 FLAG=145`
(mate on contig 9) still has real PairHMM likelihoods and is informative
ALT_REV. Java `filterNonPassingReads` (`MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`)
drops that read before PairHMM and only `addEvidence(..., 0)`s it.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: NONE`.

6R.176: production. PairHMM membership in `compute_region_read_likelihoods`
honors Java `MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE` on original BAM mate
contig (clipped copies lose `mtid`). FLAG=145 stays in Coverage so INFO DP=3
and FORMAT DP stays 2; it receives `addEvidence(..., 0)` and is not PairHMM
scored. SOR table `[0,0;1,1]` → 0.693.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: YES`.

6R.177: proof-only. At `2:92305634 G/T` FORMAT/QUAL/INFO DP/SOR stay closed.
InbreedingCoeff is Java-omitted vs Rust `1.0` because
`InbreedingCoeff.annotate` returns `emptyMap` when `n < MIN_SAMPLES=10`,
while `hc_info_values` always inserts `1 - het/n` from the single 1/1 GT.
Same one-sample population (`NA12878`). Not a formula split at this site
(Java never computes F). MQ 41.96 recorded, not investigated.
`classification: ANNOTATION_EMIT_PREDICATE_DIVERGENCE`.
`production_change: NONE`.

6R.178: production. At `2:92305634 G/T` `hc_info_values` inserts
InbreedingCoeff only when `n_genotypes >= 10` (Java
`InbreedingCoeff.MIN_SAMPLES`). One-sample NA12878 omits the key.
The `1 - het/n` formula in `annotate_hc_variant_site` is unchanged
(Java HWE `1 - het/(2pq n)` remains stacked). Header still declares
the INFO key. FORMAT/QUAL/INFO DP/SOR/FS/MQ unchanged.
`classification: ANNOTATION_EMIT_PREDICATE_DIVERGENCE`.
`production_change: YES`.

6R.180: production. At `2:92307324 TTC/T` annotations consume the per-variant
genotyping `AlleleLikelihoods` (`GenotypedSiteCall.annotation_likelihoods`)
instead of the region-wide PairHMM matrix. INFO DP=1 MQ=44 SOR=1.609.
Formulas unchanged. FORMAT/QUAL unchanged.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: YES`.

6R.181: proof-only. At `2:92307333 T/G` FORMAT/QUAL already match
(`GT=1/1 AD=0,1 DP=1 GQ=3 PL=45,3,0 QUAL=35.48`). Remaining INFO
DP/MQ/SOR is region-wide fallback because the cluster-TG early-template
path never constructs a per-variant annotation `AlleleLikelihoods`
(Java loc-loop n=1 MAPQ=44; Rust `annotation_likelihoods` empty; overlap
candidate n=6). Formulas not reached.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE` (B — missing
object construction).
`production_change: NONE`.

6R.182: proof-only. Direct Java `genotype-emit-at-loc` dump at `2:92307333 T/G`
replaces 6R.181 INFO inference. Java constructs a **new**
`AlleleLikelihoods<GATKRead,Allele>` via `hap_ll.marginalize(alleleMapper)`
from the **stored** hap matrix after `filterPoorlyModeledEvidence` (8→2),
then `retainEvidence` 2→1, then
`prepareReadAlleleLikelihoodsForAnnotation` **reuses** that object.
Retained read confirmed from the likelihood row:
`H06JUADXX130110:1:1101:10052:88682 FLAG=83 MAPQ=44`. Rust cluster-TG still
has empty `annotation_likelihoods`. Naive overlap on Rust n=8 is n=6; a
diagnostic Java-recipe (static threshold then overlap) on Rust LLs is n=1.
First causal arrow is poorly-modeled extras KEEP (stored n=8 vs Java n=2).
6R.181 B remains a later missing construction. Do not attach overlap-of-six.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE` (D — wrong upstream
genotyping subset).
`production_change: NONE`.

6R.183: proof-only. At `2:92307333 T/G` Rust `filter_poorly_modeled_region_read_likelihoods`
already matches Java (`java_equiv_keep=2`, extras KEEP does not fire). Stored n=8
because later `refresh_region_read_likelihoods(apply_normalize=false)` in
`strict_java_p12_cluster_span` overwrites that n=2 matrix. LL residual is not causal.
6R.181 B remains later. Do not attach overlap-of-six.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE` (D — wrong upstream
genotyping subset; first arrow is unfiltered post-filter refresh).
`production_change: NONE`.

6R.184: proof-only. Completes the `read_likelihoods` lifecycle at the same site:
Java-order normalize+filter is n=2, then three unfiltered P12 refreshes (hap
columns 9→11) replace it with stored n=8; cluster-TG does not write a later
object. Diagnostic replay of the existing normalize+filter on that last
refreshed matrix restores n=2 (same two Java survivors); then diagnostic
`marginalize`+`retainEvidence` is n=1 FLAG=83 MAPQ=44 best=G. Normalize does
not flip KEEP/DROP (`max_ll` unchanged). The refresh must stay; the missing
step is re-applying the existing Java-order lifecycle after the last refresh.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE` (D).
`production_change: NONE`.

6R.185: production. After the last strict-java P12-cluster
`refresh_region_read_likelihoods(..., false)`, reuse
`apply_java_order_normalize_and_filter` (existing `normalize` +
`filter_normalized_region_read_likelihoods`). Stored evidence at
`2:92307333 T/G` is n=2 (same two Java survivors). Refresh itself is
unchanged. `annotation_likelihoods` is still empty (6R.181 B, later).
FORMAT/QUAL unchanged; INFO DP/MQ/SOR/FS now reflect region-wide
fallback on n=2 (DP=1 MQ=42.05 SOR=0.693 FS=0.0), not Java n=1 MAPQ=44.
Closed `2:92307324` and `2:92305634` unchanged.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE` (D, upstream closed).
`production_change: restore Java-order filter after final P12 refresh`.

6R.186: production. On the strict-Java cluster-TG early-template path,
construct the per-variant annotation `AlleleLikelihoods` from the 6R.185
stored n=2 hap matrix: `marginalize(alleleMapper)` then
`retainEvidence(mergedVC ±2)`. At `2:92307333 T/G` the attached object
is n=1 (MAPQ=44 FLAG=83 best=G). INFO DP=1 MQ=44.00 SOR=1.609.
FORMAT/QUAL unchanged. Formulas unchanged. Closed `2:92307324`,
`2:92305634`, and `2:92316347` unchanged.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE` (B, missing
object construction).
`production_change: construct cluster-TG annotation_likelihoods`.

6R.187: proof-only. Fresh post-6R.186 INFO reconnaissance. It does not
assume that remaining annotation differences share the 6R.181/6R.186
cause. FORMAT/QUAL stay closed on chr2. First remaining genuine field
is INFO DP at `2:92305635 A/G` (Java 3 vs Rust 1). MQ 41.96 vs 53 and
SOR 0.693 vs 1.609 share that evidence object. Neighbor `2:92305634`
still has empty `annotation_likelihoods` (region-wide n=3, matches
Java). The target attaches a genotyping-subset object unique n=1
(MAPQ=53 FLAG=99) instead of Java stored-hap `retainEvidence` n=3.
Closed `2:92307333` still has `annotation_likelihoods` n=1.
`classification: A — WRONG SOURCE OBJECT`.
`production_change: NONE`.

6R.188: proof-only. At `2:92305635 A/G` the attached
`annotation_likelihoods` is the FORMAT-narrowed genotyping subset
(keep_qnames n=2 mates → 6R.180 same-QNAME collapse n=1 FLAG=99
MAPQ=53), not Java’s loc-loop stored-hap object after `marginalize` +
`retainEvidence` (n=3, includes MAPQ=25). Neighbor `2:92305634` stays
empty and luckily matches via stored fallback. 6R.186 cluster-TG is
not selected (`is_cluster_tg_snp` is `92307333` only).
`classification: A2 — WRONG SUBSET CONSTRUCTED AND ATTACHED`.
`production_change: NONE`.

6R.189: production. SiteScore annotation binds to the stored-hap
loc-loop object (`annotation_likelihoods_from_stored_haplotypes`:
`marginalize` then `retainEvidence(±2)`), not the FORMAT genotyping
subset. Target `2:92305635 A/G` is INFO DP=3 MQ=41.96 SOR=0.693.
FORMAT/QUAL stay `GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0` / `78.32`.
Neighbor `2:92305634` empty-fallback and 6R.186 cluster-TG stay.
`classification: A2 — WRONG SUBSET (closed)`.
`production_change: SiteScore stored-hap loc-loop attach`.

6R.190: proof-only. Fresh post-6R.189 INFO/VCF reconnaissance. It does
not assume another FORMAT-subset, annotation-object, or 6R.181–6R.189
cause, and it does not assume DP/MQ/SOR is next. FORMAT/QUAL stay
closed on chr2. First remaining genuine field is rust-only
`ReadPosRankSum` at `2:92305716 A/C`. Java
`RankSumTest.fillQualsFromLikelihood` on the annotation
`AlleleLikelihoods` has empty REF (hom-alt) and omits the key. Rust
still fills both lists from `region.reads` pileup and emits. Closed
`2:92305635 A/G` stays FORMAT `GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0`
QUAL `78.32` INFO DP=3 MQ=41.96 SOR=0.693 annotation n=3. Later
`2:92307359` Java-only BaseQ/MQ RankSums are not this arrow.
`classification: A — WRONG SOURCE OBJECT`.
`production_change: NONE`.

6R.191: production. ReadPosRankSum evidence binds to the same
per-variant annotation `AlleleLikelihoods` already used by SiteScore
(`GenotypedSiteCall.annotation_likelihoods` / `fillQualsFromLikelihood`).
Target `2:92305716 A/C` is REF=0 ALT=3 → undefined → key omitted.
FORMAT/QUAL stay `GT=1/1 AD=0,3 DP=3 GQ=9 PL=130,9,0` / `116.84`.
Mann-Whitney / 6R.172 `Option<f64>` / empty-list omit stay.
BaseQRankSum and MQRankSum remain header-only. Closed `2:92305635 A/G`
stays FORMAT/QUAL/DP=3/MQ=41.96/SOR=0.693 annotation n=3; 6R.186
cluster-TG stays n=1.
`classification: A — WRONG SOURCE OBJECT (closed)`.
`production_change: ReadPosRankSum from annotation_likelihoods`.

6R.192: proof-only. Fresh post-6R.191 INFO/VCF reconnaissance. It does
not assume BaseQRankSum, MQRankSum, another RankSum source-object,
`2:92305759`, or QD formatting. FORMAT/QUAL stay closed on chr2.
First remaining genuine field is Java-only `BaseQRankSum` /
`MQRankSum` at `2:92307359 CT/C` (header-only in Rust;
`hc_info_values` never inserts). ReadPosRankSum already matches.
QD 15.80 vs 15.82 is QUAL/2, not the first causal arrow. Closed
`2:92305635` / `2:92305716` / 6R.186 cluster-TG stay.
`classification: F — WRONG EMISSION PREDICATE`.
`production_change: NONE`.

6R.193: production. BaseQRankSum is emitted from the same
`fillQualsFromLikelihood` membership as ReadPosRankSum, with
`getReadBaseQualityAtReferenceCoordinate` as the element and 6R.172
`Option<f64>` (finite `0.0` emits). Target `2:92307359 CT/C` is
REF=1 ALT=1 → `Some(0.0)` → INFO insert. FORMAT/QUAL/ReadPos stay.
MQRankSum stays header-only. Closed hom-alts still omit BaseQ.
`classification: F — WRONG EMISSION PREDICATE (closed)`.
`production_change: BaseQRankSum from annotation likelihoods`.

6R.194: production. MQRankSum is emitted from the same
`fillQualsFromLikelihood` membership as ReadPos/BaseQ, with
`read.getMappingQuality()` as the element and 6R.172 `Option<f64>`.
Target `2:92307359 CT/C` is REF=1 ALT=1 unequal MAPQ → `Some(-0.674)`
→ INFO insert. FORMAT/QUAL/ReadPos/BaseQ stay. QD/QUAL are not
retuned. Closed hom-alts still omit MQRankSum.
`classification: ANNOTATION_EMIT / EVIDENCE-MEMBERSHIP PARITY (closed)`.
`production_change: MQRankSum from annotation likelihoods`.

6R.195: proof-only. QD at `2:92307359 CT/C` is Java
`QualByDepth.annotate`: `-10 * vc.getLog10PError()` / `getDepth` AD
sum, then `fixTooHighQD` (no-op below 35) then `%.2f`. Rust is the
same unrounded `ann.qual` / AD-depth. Java 15.80 is Java QUAL/2;
Rust 15.82 is Rust QUAL/2. RankSums stay closed. QUAL 31.60 vs 31.64
is not retuned.
`classification: E — NO QD DIVERGENCE: QD correctly follows an already-divergent QUAL`.
`production_change: NONE`.

6R.196: production. Biallelic QUAL AF uses Java length-based alt
Dirichlet prior (`alt.length()==refLength` → SNP else indel). Target
`2:92307359 CT/C` (PL `39,0,39`) is indel prior → unrounded QUAL
`31.60158` → 31.60; QD follows to 15.80. SNPs with the same PL stay
31.64. Emit-threshold AF unchanged.
`classification: D — WRONG PROBABILITY / QUAL FORMULA`.
`production_change: length-based AF alt prior at QUAL`.

6R.197: proof-only. Fresh post-6R.196 INFO/VCF reconnaissance. It does
not assume another QUAL prior, QD formula, RankSum emit predicate, or
6R.174 DP-source regression. FORMAT/QUAL stay closed on chr2. First
remaining genuine field is INFO DP at `2:92307403 C/A` (Java 6 vs
Rust 4), with co-observed Java-only `BaseQRankSum`/`ReadPosRankSum`.
`MQRankSum` already matches. Closed `2:92307359` stays QUAL 31.60 /
QD 15.80 / RankSums present. Rust `annotation_likelihoods` is empty
(n=0); INFO DP=4 is the region-wide emit fallback. Causal object
construction is not proven in this round.
`classification: G — UPSTREAM OBJECT / LIFECYCLE (provisional)`.
`production_change: NONE`.

6R.198: proof-only. At `2:92307403 C/A` Rust takes the cluster-downstream
early-template (`is_cluster_downstream_snp`) and returns shaped FORMAT
without constructing the Java loc-loop annotation AlleleLikelihoods.
Stored hap unique n=6 (MQ RMS 40.58 = Java MQ). Diagnostic
`marginalize` then `retainEvidence(±2)` is n=4. Java INFO DP=6 is
`Coverage.evidenceCount()` of that n=6 stored set, not the 6R.186
retain helper (which would regress MQ to ~44.15). BaseQ/ReadPos stay
omitted because the two extra stored rows are realigned off-locus.
Controls `2:92305634`/`635`/`716`/`7324`/`7333`/`7359`/`16347` stay.
`classification: E — DIFFERENT BRANCH / SPECIAL-PATH SEMANTICS`.
`production_change: NONE`.

6R.199: production. Cluster-downstream early-template attaches the
Java-equivalent loc-loop annotation object from stored unique haplotype
evidence (n=6 at `2:92307403 C/A`), not the 6R.186 retainEvidence
subset (n=4). FORMAT/QUAL stay `0/1` `2,4` `6` `72` `162,0,72` /
154.64. MQ of the attached object is 40.58. INFO DP remains 4 because
Rust Coverage still overlap-filters the n=6 object (Java
`evidenceCount` is 6); that is a new 6R.200 arrow, not combined here.
`classification: E — DIFFERENT BRANCH / SPECIAL-PATH SEMANTICS`.
`production_change: cluster-downstream stored-unique annotation attach`.

6R.200: proof-only. After 6R.199 the attached object is n=6. Java
`Coverage.annotate` is `likelihoods.evidenceCount()` (list cardinality;
no overlap). Rust `coverage_evidence_count` unique-counts then applies
`java_alignment_read_overlaps_interval(±2)`, dropping
`H06HDADXX130110:1:1101:10061:17286` FLAG=83 and
`H06JUADXX130110:1:1101:10011:51168` FLAG=81 (6→4). FORMAT/QUAL/MQ
stay. First remaining INFO DP 6 vs 4 is that second filter.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE / WRONG SECONDARY FILTER / DOUBLE OVERLAP FILTER`.
`production_change: NONE`.

6R.201: production. Coverage consumes attached unique-evidence
cardinality (`evidenceCount`) with no second ±2 overlap filter.
At `2:92307403 C/A` annotation n=6, Coverage n=6, INFO DP=6. The two
rows above remain counted. FORMAT/QUAL/MQ unchanged. 6R.199 object
construction unchanged. Remaining INFO at this site is Java-only
BaseQRankSum/ReadPosRankSum (record only).
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE / WRONG SECONDARY FILTER / DOUBLE OVERLAP FILTER`.
`production_change: Coverage unique evidenceCount; no second overlap`.

6R.202: proof-only. RankSums **are invoked** on the 6R.199 n=6 object
(MQRankSum already 1.834). BaseQ and ReadPos share empty REF fillQuals
(`getElementForRead` none on the two extra stored unique rows whose
Rust CIGARs miss the locus; those rows are MQ REF MAPQ 40 and 22).
Java emits BaseQ −1.834 / ReadPos 1.282. FORMAT/QUAL/DP/MQ stay.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: NONE`.

6R.203: production. Shared RankSum `getElementForRead` retries the
pre-realign covering CIGAR when the attached haplotype CIGAR does
not cover `vc.getStart()`. BaseQ and ReadPos share that retry.
At `2:92307403 C/A` REF lists become `[30,30]` / `[65,91]`;
BaseQRankSum=−1.834, ReadPosRankSum=1.282. MQRankSum stays 1.834.
INFO DP stays 6. FORMAT/QUAL unchanged. Next INFO split
`2:92316296 A/T` DP/MQ/SOR is recorded only.
`classification: ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE`.
`production_change: one shared getElementForRead covering-CIGAR retry`.

6R.204: proof-only. At `2:92316296 A/T` live INFO DP is Java **2** vs
Rust **3** (FORMAT/QUAL already match). The two-read hom-alt
early-template leaves `annotation_likelihoods` empty, so Coverage
falls back to region-wide stored unique n=3 including
`H06HDADXX130110:2:1101:10046:78083` FLAG=147. Java Coverage is
loc-loop `retainEvidence` n=2. MQ 47.00 vs 40.25 and SOR 2.303 vs
1.179 are the same object split (record only). 6R.201 Coverage
cardinality is not reopened.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: NONE`.

6R.205: production. The two-read hom-alt early-template attaches
`annotation_likelihoods_from_stored_haplotypes` (marginalize then
`retainEvidence(±2)`). At `2:92316296 A/T` the attached object is
n=2; Coverage consumes it (INFO DP=2); MQ=47.00 and SOR=2.303 follow
the object. FORMAT/QUAL unchanged. 6R.199 stored-unique and 6R.201
Coverage are not reused as the arrow.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: one annotation-object construction arrow`.

6R.206: proof-only. At `2:92316416 C/A` live INFO DP is Java **1** vs
Rust **3** (FORMAT `GT=1/1 AD=0,1 DP=1 GQ=3 PL=45,3,0` QUAL 35.48
already match). The one-read hom-alt early-template
(`is_mid_a_one_read_hom_alt_site`) leaves `annotation_likelihoods`
empty, so Coverage falls back to region-wide stored unique n=3
(two MAPQ=47 mates that miss ±2 plus FLAG=147 MAPQ=21). Java Coverage
is loc-loop `retainEvidence` n=1 (that MAPQ=21 row only). MQ 21.00 vs
40.25 and SOR 1.609 vs 1.179 are the same object split (record only).
6R.205's two-read arm is not this site.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: NONE`.

6R.207: production. On **only** `is_mid_a_one_read_hom_alt_site`, attach
`annotation_likelihoods_from_stored_haplotypes` (marginalize then
`retainEvidence(±2)`). At `2:92316416 C/A` the attached object is
n=1; Coverage consumes it (INFO DP=1); MQ=21.00 and SOR=1.609 follow
the object. FORMAT/QUAL unchanged. 6R.199 stored-unique, 6R.201
Coverage, and 6R.205's two-read arm are not reused as the arrow.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: one annotation-object construction arrow`.

6R.208: proof-only. At `2:92317399 C/A` live INFO DP is Java **2** vs
Rust **1** and SOR **0.693** vs **1.609** (FORMAT `GT=1/1 AD=0,2 DP=2
GQ=6 PL=90,6,0` QUAL 78.32 already match; MQ print-close 27.00).
`try_shaped` is `None` (not 6R.205/207/199). SiteScore 6R.189
retainEvidence(±2) is n=2 (same-QNAME mates FLAG=99/147 MAPQ=27).
`per_variant_annotation_likelihoods` then collapses that pair to
attached n=1. Coverage and `calculateSOR` consume that object.
6R.205/207/199 are not this site.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: NONE`.

6R.209: production. On the Java loc-loop helper
`annotation_likelihoods_from_stored_haplotypes`, do not apply the
6R.180 same-QNAME collapse. At `2:92317399 C/A` retainEvidence n=2
(same-QNAME mates) is preserved; Coverage consumes n=2 (INFO DP=2);
SOR `[0,0;1,1]`=0.693. FORMAT/QUAL/MQ unchanged. The 6R.180 collapse
stays on FORMAT-subset and stored-unique callers; 6R.186 cluster-TG
stays n=1 because retainEvidence already drops the non-overlapping
row.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: one lifecycle-scoped removal of same-QNAME collapse from the Java loc-loop annotation object`.

6R.210: production. On `gap_sparse_shaped_early` (FORMAT-shaped gap path),
attach `annotation_likelihoods_from_stored_haplotypes`. At `2:92318199 C/T`
the attached object is loc-loop retainEvidence n=1 (the overlapping FLAG=99
MAPQ=24 read); Coverage consumes it (INFO DP=1); MQ=24.00 and SOR=1.609
follow the object; MQRankSum is omitted (empty REF fillQuals). FORMAT/QUAL
unchanged. Region-wide stored unique remains n=4. 6R.209 loc-loop helper
and 6R.199/205/207 construction are not reused as the arrow.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: one annotation-object construction arrow on gap_sparse_shaped_early`.

6R.211: production. On `is_p12_phase_e_gap_het_event` (FORMAT-shaped gap-tail
het path), attach `annotation_likelihoods_from_stored_haplotypes`. At
`2:92325193 C/T` the attached object is loc-loop retainEvidence n=3
(FLAG 163/81/83, MAPQ 34/24/25); Coverage consumes it (INFO DP=3);
MQ=28.03 follows the MAPQ list; SOR `[0,1;1,1]`=0.223 was already
print-close. The same-QNAME mate FLAG=83 MAPQ=34 misses overlap and is
dropped by retainEvidence, not by 6R.180 collapse. FORMAT/QUAL
unchanged. Region-wide stored unique remains n=4. Sibling
`2:92325205 G/A` closed on the same class. 6R.209 loc-loop helper and
6R.210 gap-sparse attach are not reused as the arrow.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: one annotation-object construction arrow on is_p12_phase_e_gap_het_event`.

6R.212: production. On `event_weak_sparse_het_pl` (FORMAT-shaped weak-sparse
het path, PL 55,0,21), attach `annotation_likelihoods_from_stored_haplotypes`.
At `2:92325268 C/T` the attached object is loc-loop retainEvidence n=3
(FLAG 81/83/83, MAPQ 24/25/34); Coverage consumes it (INFO DP=3);
MQ=28.03 follows that MAPQ list. The extra stored-unique row is FLAG=163
MAPQ=34, which misses overlap (ends 92325230). This is not the 6R.211
membership: at `92325193` FLAG=163 was kept and the FLAG=83 mate was
dropped. FORMAT/QUAL (`0/1 1,2 3 21 55,0,21` QUAL 47.64) unchanged.
SOR `[0,1;0,2]`=1.179 and MQRankSum=0.000 were already print-close.
6R.209 helper and 6R.210/211 attaches are not reused as the arrow.
`classification: ANNOTATION_SOURCE_OBJECT_DIVERGENCE`.
`production_change: one annotation-object construction arrow on event_weak_sparse_het_pl`.

6R.213: proof-only. Earliest remaining common-site split is `20:29455015 G/T`.
QD 0.96 vs 1.79 is not causal: FORMAT PL already differs (`69,0,2140` vs
`122,0,2304`). GT/AD/DP and INFO MQ/FS/SOR/RankSum/DP match. EventMap at
the loc is biallelic `G/T` (neighbor `20:29455019 G/A` already matches
Java QUAL 637.64 / PL `645,0,1147`). AF calc on Java PL GLs reproduces
QUAL 61.64 and QD 0.96; AF calc on Rust GLs reproduces QUAL 114.64 and
QD 1.79. GQ 69 vs 99 is the PL second-best cap. No production change:
annotation/Coverage/PairHMM/read-filter/haplotype-construction stay as
closed in 6R.209–6R.212. Next arrow is the biallelic genotype-likelihood
object itself (allele-likelihood matrix / marginalize), not QD.
`classification: GENOTYPE_LIKELIHOOD_DIVERGENCE`.
`production_change: NONE`.

6R.214: proof-only. First divergent object is the **pre-marginalization
haplotype likelihood matrix**, not `marginalize` arithmetic. Java
`assignGenotypeLikelihoods` / PairHMM input has **30** trimmed 130 bp
haplotypes (`20:29454995-29455124`, mapper REF=18 ALT=12, evidence 88
then retainEvidence 65) and still emits PL `69,0,2140` QUAL 61.64. Rust
`call_region` has **84** haplotypes (82×181 bp + 2×212 bp, span
`29454944-29455124`, mapper REF=66 ALT=18). All 30 Java FNV window
hashes are present in Rust; Rust has 13 additional unique sequences in
that same window. Stop: Case A `ALLELE_LIKELIHOOD_INPUT_DIVERGENCE`.
Do not compensate in max-marginalize / QUAL / GQ / QD. Next arrow is
haplotype trim / extra assembled haplotypes vs Java's 30.
`classification: ALLELE_LIKELIHOOD_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.215: proof-only. First extra haplotypes at `20:29455015 G/T` enter
**after** Java-equivalent SeqGraph k=25 `findBestPaths` (78 hashes
match). Java skips k=10 (`cycles_before_dangling`) and returns K=128
n=78. Rust `merge_rt_kbest_pre_remove_paths` then extracts
RT-before-remove k=10 and inserts **44** unique `kmer_size=10`
haplotypes (after_assemble 122; trim later 122→89; `call_region` hap_n
84). Case A: not trim, not SeqGraph k-best K, not dangling-merge, not
P12 cluster supplement as the first insert. Stop:
`HAPLOTYPE_MATERIALIZATION_DIVERGENCE`. Do not disable RT merge
globally. Next arrow is a Java-equivalent restriction of that
post-SeqGraph RT merge.
`classification: HAPLOTYPE_MATERIALIZATION_DIVERGENCE`.
`production_change: NONE`.

6R.216: proof-only. Java `assembleKmerGraphsAndHaplotypeCall` never
k-bests the ReadThreadingGraph after SeqGraph `findBestPaths`.
`createGraph(k=10)` returns null (`generateSeqGraph && hasCycles`)
before dangling recovery; `assemble()` still tries every configured
kmer (not a higher-k fallback). Rust SeqGraph k=25 remains hash-
identical to Java's 78; the 44 extras originate only at
`merge_rt_kbest_pre_remove_paths` / RT-before-remove k=10 (extract
skips the cycle abort). Global RT disable is not Java-equivalent
(L2 `g2-subset-live` p11 / P12 TTC still need merge_rt). Stop:
`RT_GRAPH_LIFECYCLE_DIVERGENCE`. Production unchanged; wait for
authorization before any cycle-abort-in-merge_rt change.
`classification: RT_GRAPH_LIFECYCLE_DIVERGENCE`.
`production_change: NONE`.

6R.217: proof-only. Diagnostic `createGraph` cycle abort at RT extract
(`abort_cyclic_before_dangling` matching Java prune-then-`hasCycles`
before dangling) makes k=10 contribute 0, drops all 44 extras, and
leaves unique merge = 78; k=25 stays hash-identical. Production extract
still uses abort=false (assemble n=122). p11/indel4 acyclic k=10 ALTs
survive; P12 TTC stays n=6. Scope is extract when `use_seq_graph`, not
dump/RT-fallback. Stop: `GRAPH_STATE_DIVERGENCE`. Production unchanged.
`classification: GRAPH_STATE_DIVERGENCE`.
`production_change: NONE`.

6R.218: ONE production change. SeqGraph-path RT extract now passes
`abort_cyclic_before_dangling = use_seq_graph && abort_seq_graph_on_cycles`
(Java `createGraph`). Cyclic k=10 contributes 0; extras 0; assemble 78;
trimmed hap_n 30. Measured PL/GQ/QUAL/QD match Java (`69,0,2140` / 69 /
61.64 / 0.96). p11/indel4 acyclic k=10 retained; P12 TTC n=6;
`use_seq_graph=false` still extracts 54. `RtExtractKey` includes the
abort flag. Dump builder and SeqGraph createGraph unchanged.
`classification: GRAPH_STATE_DIVERGENCE`.
`production_change: ONE`.

6R.219: proof-only. At `20:29455379 G/A`, Java FORMAT AD is `42,5`
(FORMAT DP 47) vs Rust `44,5` (FORMAT DP 49). Java
`DepthPerAlleleBySample` is informative `bestAllelesBreakingTies` on
the retainEvidence `AlleleLikelihoods`. Rust retainEvidence remarg is
already `47,5` (n=52, UNINF=0); annotation/INFO DP use that n=52
object; FORMAT AD/PL are a different 49-row subset. PairHMM residuals
are noncausal. Allele mapping is biallelic G/A. 6R.218 stays closed.
`classification: ALLELE_LIKELIHOOD_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.220: proof-only. The 52-row retainEvidence remarg (`AD 47,5` /
`PL 68,0,1937`) **is** isolated `try_genotype_variation_event`
(`hap_events=None`). Production `call_region` FORMAT is a different
object (`AD 44,5` / `PL 78,0,1811`). No BAM/overlap/QNAME/`keep_qnames`
predicate on the 52 reproduces FORMAT. Java
`calculateGLsForThisEvent` and `DepthPerAlleleBySample` share the same
retainEvidence object; default HC has no post-retainEvidence FORMAT
subset. The three REF rows cannot be named as a filter of the 52.
Case A. 6R.218 stays closed.
`classification: ALLELE_LIKELIHOOD_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.221: proof-only. Isolated and production-arg `try_genotype` stay
`AD 47,5` / `PL 68,0,1937`. `assign_genotype_likelihoods_for_region`
loc-loop is the FORMAT object `AD 44,5` / `PL 78,0,1811` (B = C).
Colocated merge does not fire. Event identity is unchanged. Annotation
unique n=52 on both (not a 52→49 row filter). Assign visits EventMap
loc 29455375 before the target with no emitted call (next inner).
`classification: GENOTYPE_LIKELIHOOD_LIFECYCLE_DIVERGENCE`.
`production_change: NONE`.

6R.222: proof-only. Loc-loop order through the target is
`[29455375, 29455379]`. Predecessor `20:29455375 T/A` REJECTS
(`VariantNotConfident`) but writes TLS region-likelihood rows
`(len=2808, n_haps=54, n_rows=52)`. Target `SiteScore` then hits that
cache and emits FORMAT `AD 44,5` / `PL 78,0,1811`. Skip-pred and
target-first stay A' (`47,5` / `68,0,1937`); filling the full PairHMM
table after the predecessor restores A'. Annotation 52×54 cells are
unchanged (`changed=0`). No second PairHMM. Java has no pointer-keyed
row cache between neighboring events.
`classification: SHARED_STATE_LIFECYCLE_DIVERGENCE`.
`production_change: NONE`.

6R.223: proof-only. First inner operation is
`with_region_likelihood_rows` cache HIT `(ptr, len=2808, n_haps=54)`
inside `SiteScore::from_allele_mapping`. Dense 52×54 hash matches the
predecessor T/A table; 52×2 membership differs by 3 rows (intersection
cells unchanged). AD `47,5`→`44,5` and PL `68,0,1937`→`78,0,1811`
change together in `genotype_from_marginalized_rows`. No second
PairHMM. Java has no equivalent pointer-keyed row cache.
`classification: GENOTYPE_LIKELIHOOD_CACHE_DIVERGENCE`.
`production_change: NONE`.

6R.224: proof-only. Predecessor `20:29455375 T/A` and target
`20:29455379 G/A` produce the same TLS key `(ptr, len=2808, n_haps=54)`
because successive owned retainEvidence subset Vecs of length 2808 reuse
the same allocator address after the predecessor Vec is dropped. Sparse
cells and read-index membership differ; the HIT returns predecessor
dense rows. Java has no pointer-keyed row cache.
`classification: CACHE_STORAGE_ALIASING`.
`production_change: NONE`.

6R.225: proof-only. `with_region_likelihood_rows` stores an owned dense
`Vec<ReadLikelihoodRow>` copy; the cache value outlives the source subset.
The key is `(as_ptr, len, n_haps)` — storage identity of a short-lived
retainEvidence Vec, not sparse-cell identity. Isolated predecessor has a
legitimate intra-event HIT (reshape then SiteScore on the same live slice).
Production `assign` on this region records 34 MISS + 1 invalid alias HIT
and 0 legitimate HITs. Diagnostic disable and a content-hash key both
restore P2 dense `0x67d17c7b8ef58d07`. Extra independent storage cannot
prevent the alias because the value is already owned. Java has no
equivalent cache. Decision Case B: smallest semantic identity is the
sparse-cell population plus `n_haps`. Uncached rebuild is semantically
sufficient (G1 KEEP is performance).
`classification: CACHE_KEY_IDENTITY_DIVERGENCE`.
`production_change: NONE`.

6R.226: one production change. `with_region_likelihood_rows` keys the
exact sparse-cell sequence `(read_index, haplotype_index, log10.to_bits())`
plus `n_haps`. Allocator reuse of a dropped 2808-cell subset Vec no
longer returns predecessor dense rows. Canonical `20:29455379 G/A`
dense hash is P2 `0x67d17c7b8ef58d07`; FORMAT `AD=47,5` `PL=68,0,1937`
`GQ=68` `DP=52` (P2 remarg). Java FORMAT remains `AD=42,5` `PL=84,0,1738`
`GQ=84` `DP=47` — that remaining membership/PL split is not patched here.
Intra-event same-slice HIT is preserved. Cache disabling was not chosen.
`classification: CACHE_KEY_IDENTITY_DIVERGENCE`.
`production_change: ONE`.

6R.227: proof-only. After 6R.226, isolated `try_genotype` and production
`assign_genotype_likelihoods_for_region` both consume P2 dense
`0x67d17c7b8ef58d07` and the identical 52×2 allele-likelihood matrix
(`AD=47,5` `PL=68,0,1937`). Predecessor `20:29455375 T/A` does not mutate
that matrix. Remaining vs Java is the 52-row retainEvidence remarg versus
Java's 47-read `AlleleLikelihoods` (five extra REF votes). Cache key
unchanged.
`classification: GENOTYPE_LIKELIHOOD_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.228: proof-only. Live Java `genotype-emit-at-loc` captures the actual
post-`retainEvidence` set at `20:29455379 G/A`: 47 reads, interval
`20:29455377-29455381`, votes AD `42,5`, PL `84,0,1738`. Rust event-local
overlap remarg is 52. Intersection 47, Rust-only 5, Java-only 0. The five
extras are informative REF and absent from Java hap_ll (n=236); original-BAM
mate-contig passes (6R.176 eliminated). First drop is before retainEvidence.
Cache key and `genotype_from_marginalized_rows` unchanged.
`classification: READ_FILTER_MEMBERSHIP_DIVERGENCE`.
`production_change: NONE`.

6R.229: proof-only. Live Java `hap-ll-membership-at-loc` tracks the five
Rust-only REF reads through hap_ll construction at `20:29455379 G/A`.
They are present through `normalizeLikelihoods` (245, bitmap `11111`)
and absent after `filterPoorlyModeledEvidence` (236, `00000`). Region,
clip, stub, `filterNonPassingReads`, mate-contig, PairHMM input, and
`AlleleLikelihoods` construction all KEEP. Java `max_ll < -8` DROP;
Rust same thresh KEEP (`max_ll ≥ -8`, `extra_retain=false`).
Cache key, `genotype_from_marginalized_rows`, and retainEvidence
unchanged.
`classification: POORLY_MODELED_FILTER_DIVERGENCE`.
`production_change: NONE`.

6R.230: proof-only. Live Java `five-ll-at-loc` vs Rust PairHMM evidence at
`20:29455379 G/A`. All five reads enter both PairHMM engines. Java raw
`max_ll` already `< -8` (`prim_max == norm_max`). Rust PairHMM sequences
are clipped to start `29455355`; Java keeps padded-window `148M` /
`33H113M2H`. Haplotype FNV sets are disjoint (Java 60 vs Rust 54).
300 Java prim cells have no matching Rust hap identity. Not the −8
predicate and not the ~1e-6 GKL-float residual.
Cache key, `genotype_from_marginalized_rows`, retainEvidence, and the
poorly-modeled formula unchanged.
`classification: READ_SEQUENCE_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.231: proof-only. The five Rust PairHMM sequences at `20:29455379 G/A`
are physical hard-clips to `region.extended_start=29455355`, produced by
`clip_finalized_reads_in_place` / `hard_clip_to_region` after
`AssemblyRegionTrimmer::apply_trim`. Same Java clip primitive
(`ReadClipper.hardClipToRegion`); Java interval is padded variant
`20:29455294-29455584` (`29455314 − 20`). Rust `29455355 = 29455375 − 20`
because EventMap `variation_events` lacks Java's leftmost SNP `29455314`.
PairHMM is downstream. Haplotype-population arrow deferred at that
interval dependency.
Cache key, `genotype_from_marginalized_rows`, retainEvidence, and the
poorly-modeled formula unchanged.
`classification: READ_INTERVAL_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.232: proof-only. Java variation event `20:29455314` is SNP `G>C` on a
single assembled haplotype (`c7acc50dfb9f9ecc`, CIGAR `460M`). Rust's
128 untrimmed haplotypes are all `G` at that coordinate, so EventMap
union never sees the event. Diagnostic: adding the SNP to Rust trim
restores left pad `29455294`; removing only that Java event does **not**
make Java match Rust (`29455328` remains). Clip/padding/PairHMM
downstream. Do not patch EventMap.
`classification: HAPLOTYPE_CONSTRUCTION_DIVERGENCE`.
`production_change: NONE`.

6R.233: proof-only. The missing carrier is already a k=25 SeqGraph
source→sink path in both implementations (C-edge multiplicity **3**;
post-prune RT graphs match 1064/1079). Java k-best `K=128` selects
`c7acc50dfb9f9ecc` at rank **116**. Rust production `K=128` returns
G at `29455314` on 128/128 paths; diagnostic `K=256` (not a product
change) returns the **same FNV at rank 129**. One of four pileup-C
reads does not thread C k-mers. EventMap/clip/PairHMM stay
downstream. Do not raise `K`.
`classification: ASSEMBLY_KBEST_DIVERGENCE`.
`production_change: NONE`.

6R.234: proof-only. Same carrier FNV `c7acc50dfb9f9ecc` is **not** the
same k-best path state: Java 18 edges vs Rust 19. Prefix 0–16 matches
(`log10(mult/out)` on total multiplicity). Rust splits the 38 bp Java
sink into 30+8 bp and adds `log10(63/78)=−0.09275`, moving rank
116 / −2.64786702 → 129 / −2.74062107. Production `K=128` is a
downstream cutoff. Do not raise `K`.
`classification: ASSEMBLY_KBEST_STATE_IDENTITY_DIVERGENCE`.
`production_change: NONE`.

6R.235: proof-only. The 38 vs 30+8 split is **graph topology**, not
k-best. Rust's 1 bp SeqGraph already joins at the first base of
`TGTTTCTT` (`in_degree=2`, plus a `GT` bubble). Initial
`zip_linear_chains` emits the independent 8 bp sink; Java's sink is one
38 bp vertex 36 with 4 incoming. The 30 bp vertex appears later at
`merge_common_suffices`. `63/78` is a consequence. `K=128` stays
downstream. Do not raise `K`.
`classification: GRAPH_TOPOLOGY_DIVERGENCE`.
`production_change: NONE`.

6R.236: proof-only. The 1 bp join is an **edge insertion**, not a
different vertex-identity rule. `971 → 428` (support 1, non-ref) is
absent after raw threading and prune; dangling-tail recovery
(`addEdge` weight 1, `dangling_java_exact=true`) splices it onto
reference k-mer 428, making `in_degree=2`. Java's post-cleanup
SeqGraph has no `TGTTTCTT` / `GT` vertices (one 38 bp sink). `K=128`
stays downstream. Do not raise `K`.
`classification: EDGE_INSERTION_DIVERGENCE`.
`production_change: NONE`.

6R.237: proof-only. Same dangling-tail candidate and SW CIGAR `3I9M`.
Java `mergeDanglingTail` computes `refIndexToMerge = 8 - 9 + 1 = 0`
and returns 0 (LCA-cycle guard). Rust `saturating_sub` yields 1 and
`addEdge(971,428,weight=1)`. `K=128` stays downstream. Do not raise
`K`.
`classification: DANGLING_TAIL_MERGE_PREDICATE_DIVERGENCE`.
`production_change: NONE`.

6R.238: proof-only. Java `refIndexToMerge == 0` is a **path-index
sentinel** (LCA / no `addEdge`), not a graph vertex id. Signed
`lastRefIndex - matchingSuffix + 1` can be 0; Rust
`saturating_sub + 1` maps that underflow onto path index 1. Negative
raw is unreachable from `longestSuffixMatch`. Graph vertex 0 can
exist independently (carrier ref source). `K=128` stays downstream.
Do not raise `K`.
`classification: INDEX_SENTINEL_SEMANTICS_DIVERGENCE`.
`production_change: NONE`.

6R.239: production. Replace saturating unsigned subtraction in
`plan_dangling_tail_merge` with checked
`(last_ref_idx + 1).checked_sub(matching_suffix)` so Java path-index
`0` remains the LCA / no-splice sentinel (`ref_index_zero_cycle`).
Canonical `8 - 9 + 1 = 0`: splice `971→428` absent, 38 bp sink,
carrier `c7acc50dfb9f9ecc` score `-2.64786702`, K256/K128 rank 116.
`K=128` was not changed.
`classification: INDEX_SENTINEL_SEMANTICS_DIVERGENCE`.
`production_change: ONE` (`dangling_tail_ref_index_to_merge`).

6R.240: proof-only. After 6R.239, Java and Rust EventMaps both contain
`20:29455314 G>C` on the canonical carrier (score `-2.64786702`,
K=128 unchanged). Pinned Java covering VCF still omits the locus.
Live Rust GLs `PL=21,0,1461` pass `passesEmitThreshold` at
`stand_emit=10` (QUAL 13.63) and fail Java `stand-call-conf=30`
(AF-monomorphic). First operation: emit predicate, not assembly.
`classification: EMISSION_PREDICATE_DIVERGENCE`.
`production_change: NONE`.

6R.241: production. `HcGenotypingConfig::strict_java` supplies Java 4.4
`standardConfidenceForCalling=30` to the existing `passesEmitThreshold`
equivalent (`java_emit_would_pass` / `filter_genotyped_calls_for_strict_java_emit`).
Non-strict / legacy keep `stand_emit=10`. Covering G>C GLs stay
`GT 0/1 PL 21,0,1461` QUAL ~13.63 and are no longer VCF-emitted.
`20:29455379 G/A` FORMAT is unchanged. K=128 unchanged.
`classification: EMISSION_PREDICATE_DIVERGENCE`.
`production_change: ONE` (`strict_java` calling confidence).

6R.242: proof-only. Fresh post-6R.241 covering VCF inventory. Covering
`20:29455314 G>C` is omitted by both engines. Totals Java 126 / Rust 126,
java-only 2, rust-only 2. First genuine remaining split is a **common**
record `20:29455649 T/TGTTTG` INFO DP 123 vs 230 (PL ±1 representation;
GT/AD/QUAL match; both emit at stand=30). Live `annotation_likelihoods`
is empty (`ann_n=0`); live emit INFO DP=230 is the region-wide fallback.
Do not jump to later rust-only `20:29456196 A/T`.
`classification: GENOTYPE_LIKELIHOOD_INPUT_DIVERGENCE` (annotation
AlleleLikelihoods, not genotyping PL).
`production_change: NONE`.

6R.243: proof-only. At `20:29455649 T/TGTTTG` colocated merge already
builds retainEvidence unique n=123 (Java INFO DP) then constructs
`GenotypedSiteCall` with `annotation_likelihoods=Vec::new()`.
`merged_handled_locs` skips SiteScore attach. Emit falls back to
stored-hap unique 230. PL ±1 is non-causal. Not emission.
`classification: ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE`.
`production_change: NONE`.

6R.244: proof-only. Why merge discards the 123-read object. Local
`subset` (retainEvidence n=123) is in scope at the Call constructor and
is not copied onto `annotation_likelihoods`. Java 4.4
`prepareReadAlleleLikelihoodsForAnnotation` reuses the genotyping
AlleleLikelihoods (contamination off); it does not build a new object.
`merged_handled_locs` then skips SiteScore attach (correct one-loc
merged VC; does not wipe a populated field). Empty emit sentinel →
INFO DP 230. PL ±1 independent. Candidates A+C (E as interaction).
`classification: ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE`.
`production_change: NONE`.

6R.245: proof-only. Minimal Java-equivalent attach is
`subset.into_owned()` onto `GenotypedSiteCall.annotation_likelihoods`
after `hap_rows` (Cow::Owned move; no second PairHMM). Empty sentinel
suppression is sufficient: Coverage of the subset identity is 123.
`production_change: NONE`. Proposed 6R.246: that single field assignment.

6R.246: production. `try_genotype_colocated_snp_indel_merge` Call sets
`annotation_likelihoods: subset.into_owned()` after `hap_rows`. Unique
read_index 123 / 2952 cells; INFO DP 230 → 123. FORMAT GT/AD/DP/GQ and
QUAL unchanged. PL ±1 retained. Fallback sentinel implementation
unchanged.
`classification: ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE`.

6R.247: proof-only. Remaining PL `570,0,3518` vs Java `570,0,3517` at
`20:29455649 T/TGTTTG` is not rounding. Rust 1/1 log10 GL
`-351.75162674083281900` has continuous PL `3517.516…`; both
`f64::round` and Java `Math.round`/`GLsToPLs` of those GLs emit 3518.
Java VCF 3517. First unequal value is the calculator 1/1 (merged 2/2)
genotype likelihood.
`classification: GENOTYPE_LIKELIHOOD_DIVERGENCE`.
`production_change: NONE`.

6R.248: proof-only. Merged allele 2 is `TGTTTG`; 2/2 is emitted 1/1.
Homozygous calculator is Java-equivalent: `GL(2/2)=Σ L(read|allele2)`
(Δ vs column sum 4.55e-13). Evidence 123×3. Allele floor and mapper
pools do not move PL. The 3517.5 crossing is already in the allele-row
`L(read|TGTTTG)` column.
`classification: GENOTYPE_LIKELIHOOD_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.249: proof-only. The TGTTTG allele-row is `max` over five EventMap
`T/TGTTTG` assembly haplotypes (idx 19–23, FNV `55012fcf3b430591` …
`79451c576721a729`). Rust `pool_max_log10` matches Java 4.4
`AlleleLikelihoods.marginalize` (max, not log-sum-exp). Reconstruction
matches all 123 cells. Allele-level −4.5 floor is a no-op. Java hap-list
dump for this ActiveFull is absent. Remaining causal inputs are the five
PairHMM columns.
`classification: ALLELE_LIKELIHOOD_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.250: proof-only. The five TGTTTG PairHMM columns share one Java-equivalent
read/quality plane (BQ 18, PCR Conservative, GOP/GCP). Kernel is NEON f64
vs Java GKL float. Known 7.41e-6 residual stacked over 123 reads cannot
cross 3517.5. Remaining candidate is the five haplotype sequences vs Java
(executable dump absent).
`classification: HAPLOTYPE_INPUT_DIVERGENCE`.
`production_change: NONE`.

6R.251: proof-only. The five Rust PairHMM 161-mers are EventMap-exact,
share CIGAR insertion `GTTTG`, and are trim subsequences of 388-bp k-best
parents. They differ from each other only by flank SNPs, not the insertion.
Java 4.4 construction/trim operations match; Java's concrete 161-mers are
unverified (executable absent).
`classification: HAPLOTYPE_SEQUENCE_INTERNALLY_CONSISTENT`.
`production_change: NONE`.

6R.252: recovered official GATK 4.4.0.0 JAR and ran `hap-trim-at-loc` on
the parity BAM/REF. Untrimmed five TGTTTG parents are byte-identical.
Trimmed haplotypes are not: Java `20:29455560-29455728` 174 bp
`90M5I79M` vs Rust `20:29455569-29455724` 161 bp `81M5I75M`; Rust is
Java[9:-4].
`classification: HAPLOTYPE_TRIM_DIVERGENCE`.
`production_change: NONE`.

6R.253: Java executable dump of `AssemblyRegionTrimmer.trim`. Variant span
matches (`29455590–29455703`). Java pads `A/AT` @ 29455644 with STR 84
(`75+longestSTR=9`); Rust keeps indel pad 75 (insertion STR skipped).
That is the 9-left / 4-right trim-span split.
`classification: TRIM_PADDING_SEMANTICS_DIVERGENCE`.
`production_change: NONE`.

6R.254: measurement-only. Frozen 123-read retainEvidence + actual Java
174-mer TGTTTG haplotypes through Rust NEON-f64 PairHMM. Control strip
of `CAAAGAGTA`/`TAAA` reproduces the 161-mer baseline (PL 3518). Java
174-mers apply a uniform ~−0.0337 log10 penalty, move emitted 2/2 GL
from −351.7516 to −351.8965, and integer PL **3518 → 3519**, away from
3517.5. Winner sets unchanged (12 unique-max rows keep H3/H1/H4).
`classification: TRIM_SPAN_MOVES_GL_BUT_NOT_PL_BOUNDARY`.
`production_change: NONE`.

6R.255: measurement-only. Frozen 123 reads and Rust 161-mer TGTTTG
haplotypes; only `hardClipToRegion` span varies. Rust clip reproduces
PL 3518. Java clip `20:29455560-29455728` adds 739 read bases, moves
emitted 2/2 from −351.7516 to −428.0464 (integer PL **4280**), away
from 3517.5. Left-only and right-only both hurt (3942 / 3711). 22/123
rows shift; membership unchanged.
`classification: READ_CLIP_INTERVAL_MOVES_GL_BUT_NOT_PL_BOUNDARY`.
`production_change: NONE`.

6R.256: measurement-only. Frozen 123 reads with Java 174-mer TGTTTG
haplotypes **and** Java clip `20:29455560-29455728` together. Three
controls reproduce (PL 3518 / 3519 / 4280). Joint emitted 2/2 is
−351.9157 (continuous PL 3519.157, integer **3519**), still above
3517.5. The two isolated effects are non-additive
(Δ_interaction = +76.2755 log10): 174-mer flanks cancel almost all of
the 6R.255 clip overhang penalty.
`classification: JOINT_JAVA_INPUTS_MOVE_GL_BUT_NOT_PL_BOUNDARY`.
`production_change: NONE`.

6R.257: measurement-only. Frozen Java 174-mers + Java clip; 123/123
clipped read arrays, five haplotype FNV, BQ/PCR/GCP match Java 4.4
`modifyReadQualities`. IQ/DQ Q6 floor never fires (0 bytes `< 6`).
Joint diagnostic remains PL 3519. Remaining difference is kernel
configuration (GKL vs NEON f64), not trim/clip.
`classification: PAIRHMM_CONFIGURATION_DIVERGENCE`.
`production_change: NONE`.

6R.258: measurement-only. Frozen 6R.257 inputs. First backend value is
GKL float `1.0f32/174` vs Rust `1.0f64/174` (bits `0x3f778a4c80000000` vs
`0x3f778a4c8178a4c8`). Substituting only that haplen quotient leaves
integer PL 3519. NEON is bit-identical to scalar logless on this matrix.
`classification: PAIRHMM_RECURRENCE_DIVERGENCE_NOT_SUFFICIENT`.
`production_change: NONE`.

6R.259: measurement-only. Frozen 6R.257 inputs. First remaining ph2pr-chain
value is GKL `-((float)32)/10.f` vs Rust `-(32 as f64)/10.0` (bits
`0xc04ccccd` vs `0xc00999999999999a`). Substituting only that exponent
moves continuous PL by −1.00e-6; integer PL stays 3519.
`classification: PAIRHMM_NUMERICAL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.260: measurement-only. Frozen 6R.257 inputs. First isolated `powf`
split is Q=20 (exponent −2.0 exact both sides): GKL `powf(10.f, −2.f)`
f32 `0x3c23d70a` vs Rust `10f64.powf(−2.0)` `0x3f847ae147ae147b`.
Injecting only the stored GKL ph2pr moves continuous PL by −8.26e-7;
integer PL stays 3519. `1-ph2pr` diagnostic moves *away* from Java.
`classification: PAIRHMM_POWF_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.261: measurement-only. Frozen 6R.257 inputs. GKL match is AVX
`VEC_SUB(1.0, distm)` in float (`_1_distm`), not a Context table.
Q=20: f32 `0x3f7d70a4` vs Rust `1.0-0.01` `0x3fefae147ae147ae`.
Injecting only that match moves continuous PL by +3.77e-6 (away from
Java); integer PL stays 3519.
`classification: PAIRHMM_MATCH_PRIOR_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.262: measurement-only. Frozen 6R.257 inputs. GKL mismatch is AVX
`VEC_DIV(distm, 3.0)` in float. Q=20: f32 `0x3b5a740d` vs Rust
`err/3` `0x3f6b4e81b4e81b4f` (`new-f32-division-rounding`). Injecting
only that mismatch moves continuous PL by −7.94e-7; integer PL stays
3519. `classification: PAIRHMM_MISMATCH_PRIOR_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.263: measurement-only. Frozen 6R.257 inputs. GKL AVX distm is
`VEC_BLENDV` / `_mm256_blendv_ps(mismatch, match, mask)` (AVX-512:
`_mm512_mask_blend_ps`). Bit-preserving: output bits equal the selected
match `0x3f7d70a4` or mismatch `0x3b5a740d` at Q=20. Injecting only
those selected distm values moves continuous PL by +2.97e-6 (away from
Java); integer PL stays 3519.
`classification: PAIRHMM_DISTM_BLEND_DIVERGENCE_NOT_CAUSAL`.
First new primitive: `NONE`. `production_change: NONE`.

6R.264: measurement-only. Frozen 6R.257 inputs. GKL AVX M-update is
`computeMXY` outer `VEC_MUL(sum, distmSel)` = `_mm256_mul_ps` (float).
First-cell-only injection of the INITIAL-normalized f32 product moves
continuous PL by −4.60e-8; all-M-multiply GKL-f32 CF moves +1.96e-6;
integer PL stays 3519.
`classification: PAIRHMM_M_MULTIPLY_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.265: measurement-only. Frozen 6R.257 inputs. GKL AVX X-update first
primitive is `computeMXY` `VEC_MUL(M_t_1, pMX)` = `_mm256_mul_ps`
(float). At `(i=1,j=1)` the product is 0. First-meaningful-cell
injection of the INITIAL-normalized f32 product leaves continuous PL
bit-identical (integer PL 3519).
`classification: PAIRHMM_X_UPDATE_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.266: measurement-only. Frozen 6R.257 inputs. GKL AVX X-update second
primitive is `computeMXY` `VEC_MUL(X_t_1, pXX)` = `_mm256_mul_ps`
(float, `pXX=ph2pr[gcp]`). At `(i=1,j=1)`/`(i=2,j=1)` the product is 0.
First-meaningful-cell injection of the INITIAL-normalized f32 product
leaves continuous PL bit-identical (integer PL 3519).
`classification: PAIRHMM_X_XX_MUL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.267: measurement-only. Frozen 6R.257 inputs. GKL AVX X-update third
primitive is `computeMXY` `VEC_ADD` = `_mm256_add_ps` (float). At
`(i=3,j=1)` one addend is 0. First both-nonzero-cell injection of the
INITIAL-normalized f32 add leaves continuous PL bit-identical (integer
PL 3519).
`classification: PAIRHMM_X_ADD_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.268: measurement-only. Frozen 6R.257 inputs. GKL AVX Y-update first
primitive is `computeMXY` `VEC_MUL(M_t_1_y, pMY)` = `_mm256_mul_ps`
(float, `pMY=ph2pr[del]`). At `(i=1,j=1)` the product is 0.
First-meaningful-cell injection of the INITIAL-normalized f32 product
leaves continuous PL bit-identical (integer PL 3519).
`classification: PAIRHMM_Y_MY_MUL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.269: measurement-only. Frozen 6R.257 inputs. GKL AVX Y-update second
primitive is `computeMXY` `VEC_MUL(Y_t_1, pYY)` = `_mm256_mul_ps`
(float, `pYY=ph2pr[gcp]`). At `(i=1,j=1)`/`(i=1,j=2)` the product is 0.
First-meaningful-cell injection of the INITIAL-normalized f32 product
leaves continuous PL bit-identical (integer PL 3519).
`classification: PAIRHMM_Y_YY_MUL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.270: measurement-only. Frozen 6R.257 inputs. GKL AVX Y-update third
primitive is `computeMXY` `VEC_ADD` of the two Y products = `_mm256_add_ps`
(float, not FMA). At `(i=1,j=1)`/`(i=1,j=2)` the add is trivial.
First-meaningful-cell injection of the INITIAL-normalized f32 add
leaves continuous PL bit-identical (integer PL 3519).
`classification: PAIRHMM_Y_ADD_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.271: measurement-only. Frozen 6R.257 inputs. GKL AVX M-update first
inner product is `computeMXY` `VEC_MUL(M_t_2, pMM)` = `_mm256_mul_ps`
(float, `pMM=set_mm_prob(ins,del)`). At `(i=1,j=1)`/`(i=2,j=1)` the product
is 0. First-meaningful-cell injection of the INITIAL-normalized f32
product moves continuous PL by 9.15e-8 toward Java (integer PL 3519).
`classification: PAIRHMM_M_UPDATE_MM_MUL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.272: measurement-only. Frozen 6R.257 inputs. GKL AVX M-update second
inner product is `computeMXY` `VEC_MUL(X_t_2, pGAPM)` = `_mm256_mul_ps`
(float, `pGAPM=1.0f-ph2pr[gcp]`). At `(i=1,j=1)`/`(i=2,j=1)`/`(i=2,j=2)`
the product is 0. First-meaningful-cell injection of the
INITIAL-normalized f32 product leaves continuous PL bit-identical
(integer PL 3519).
`classification: PAIRHMM_M_UPDATE_XGAPM_MUL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.273: measurement-only. Frozen 6R.257 inputs. GKL AVX M-update third
inner product is `computeMXY` `VEC_MUL(Y_t_2, pGAPM)` = `_mm256_mul_ps`
(float, `pGAPM=ctx._(1.0)-ph2pr[gcp]`; stripe-0 `Y_t_2=VEC_SET_LSE(init_Y)`).
First meaningful cell is `(i=1,j=1)`. First-meaningful-cell injection of
the INITIAL-normalized f32 product moves continuous PL 3.23e-8 away from
Java (integer PL 3519).
`classification: PAIRHMM_M_UPDATE_YGAPM_MUL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.274: measurement-only. Frozen 6R.257 inputs. GKL AVX M-update first
inner add is `computeMXY` `VEC_ADD(M*pMM, X*pGAPM)` = `_mm256_add_ps`
(float, not FMA). At `(i=1,j=1)`/`(i=2,j=2)`/`(i=3,j=2)` the add is
trivial. First both-nonzero-cell injection of the INITIAL-normalized
f32 add moves continuous PL 2.75e-7 away from Java (integer PL 3519).
`classification: PAIRHMM_M_UPDATE_MM_XGAPM_ADD_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.275: measurement-only. Frozen 6R.257 inputs. GKL AVX M-update second
inner add is `computeMXY` `VEC_ADD(partial_M, Y*pGAPM)` = `_mm256_add_ps`
(float, not FMA; left-assoc `(MM+XGAPM)+YGAPM`). At `(i=1,j=1)`/
`(i=2,j=2)`/`(i=3,j=3)` the add is trivial. First both-nonzero cell is
`(i=2,j=3)`. Y is not sub-ULP of the M partial. First-both-nonzero-cell
injection of the INITIAL-normalized f32 add moves continuous PL
1.36e-12 toward Java (integer PL 3519).
`classification: PAIRHMM_M_UPDATE_MXGAP_YGAPM_ADD_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.276: measurement-only. Frozen 6R.257 inputs. GKL AVX M-update outer
multiply is `computeMXY` `VEC_MUL(closed_6R275_f32_inner_sum, distmSel)`
= `_mm256_mul_ps` (float, not FMA). The left operand is the exact 6R.275
f32 inner sum. First nonzero cell is `(i=1,j=1)` (Q=32 match).
First-nonzero-cell injection moves continuous PL 4.60e-8 toward Java
(integer PL 3519).
`classification: PAIRHMM_M_UPDATE_SUM_DISTM_MUL_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.277: measurement-only. Frozen 6R.257 inputs. GKL AVX last-stripe
accumulation is `compute_full_prob` `sumM = VEC_ADD(sumM, M_t.d)` =
`_mm256_add_ps` (float, not FMA). The first both-nonzero result-lane add
is read 0, haplotype 0, column 116 (prior nonzero column 115); both
operands are subnormal and the f32 add is exact. First-both-nonzero
injection moves continuous PL 0 toward Java (integer PL 3519).
`classification: PAIRHMM_LAST_STRIPE_SUMM_ADD_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.278: measurement-only. Frozen 6R.257 inputs. GKL AVX last-stripe
accumulation is `compute_full_prob` `sumX = VEC_ADD(sumX, X_t.d)` =
`_mm256_add_ps` (float, not FMA). The first both-nonzero result-lane add
is read 0, haplotype 0, column 116 (prior nonzero column 115); both
operands are subnormal and the f32 add is exact. First-both-nonzero
injection moves continuous PL 0 toward Java (integer PL 3519).
`classification: PAIRHMM_LAST_STRIPE_SUMX_ADD_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.279: measurement-only. Frozen 6R.257 inputs. GKL AVX last-stripe
final add is `compute_full_prob` `sumMX.d = VEC_ADD(sumM, sumX)` =
`_mm256_add_ps` (float, not FMA). The first both-nonzero result lane is
read 0, haplotype 0 (read length 122, `remainingRows = 2`, lane 1). Both
operands are normal f32. The isolated f32 add differs from the f64 add of
the same widened operands by 41484288 ULP and moves continuous PL
5.06e-7 away from Java (integer PL 3519).
`classification: PAIRHMM_LAST_STRIPE_SUMMX_ADD_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.280: measurement-only. Frozen 6R.257 inputs. GKL AVX last-stripe
lane extract is `compute_full_prob` `result_avx2 = sumMX.f[remainingRows-1]`.
For read length 122, `remainingRows = 2` and the extracted lane is 1.
The read copies f32 bits `0x0cea27f6` with no arithmetic, conversion,
rounding, or narrowing. Substituting that exact lane moves continuous PL
5.06e-7 away from Java (integer PL 3519).
`classification: PAIRHMM_LAST_STRIPE_SUMMX_LANE_EXTRACT_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.281: measurement-only. Frozen 6R.257 inputs. GKL compares the
extracted f32 `0x0cea27f6` with `MIN_ACCEPTED` `1e-28f` (`0x10fd87b6`)
by `if (result_float < MIN_ACCEPTED)`. Both operands stay f32. The
boolean is true, and the same predicate on the Rust f64 ratio is true
on all 610 matrices. Forcing that decision without replacing the
likelihood moves continuous PL 0 (integer PL 3519).
`classification: PAIRHMM_RESULT_FLOAT_MIN_ACCEPTED_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.282: measurement-only. Frozen 6R.257 inputs. The first accepted
`result_float` is retain 0, haplotype 1, bits `0x12d79844`
(`result_float < 1e-28f` is false). GKL calls `log10f` from `<math.h>`.
On the GATK 4.4 Ubuntu 18.04 oracle that is glibc 2.27
`__ieee754_log10f`, not this host's libm. For this input the SSE2 and
FMA-contracted forms both return `0xc1d6ee20`. Injecting only that
logarithm's difference from `f64::log10` of the same bits moves
continuous PL `4.72e-7` toward Java (integer PL 3519).
`classification: PAIRHMM_RESULT_FLOAT_LOG10F_DIVERGENCE_NOT_CAUSAL`.
`production_change: NONE`.

6R.283: measurement-only. Frozen 6R.257 inputs. Coarse checkpoints
compare scale-matched GKL `log10f(result_float)` with the Rust log10
likelihood for all five TGTTTG haplotypes, then the floor-max aggregate,
the biallelic genotype likelihoods, and continuous PL. Injecting the
GKL vector into the existing Rust downstream path moves continuous PL
`3.40e-8` away from Java (integer PL 3519). `log10(2^120)` is not
subtracted: the reconstructed f32 is already probability-scale.
`classification: PAIRHMM_COARSE_CHECKPOINT_PAIRHMM_SCORES_NOT_MATERIAL`.
`production_change: NONE`.

6R.284: measurement-only. Frozen 6R.257 inputs were scored by the stock
`libgkl_pairhmm.so` inside `broadinstitute/gatk:4.4.0.0` on
`--platform linux/amd64` (glibc 2.27, AVX, `useDoublePrecision=false`,
1 thread). All 610 scored matrices took the float branch
(`result_float < 1e-28f` was false). Injecting that matrix into the
existing Rust downstream moves continuous PL `2.46e-6` away from the
pinned Java VCF boundary (integer PL 3519).
`classification: PAIRHMM_LIVE_GKL_LIKELIHOOD_BOUNDARY_NOT_MATERIAL`.
`production_change: NONE`.

6R.285: measurement-only. One GATK 4.4.0.0 genotyping pass on the frozen
locus. The five TGTTTG haplotypes are indices 19–23. The event is
`T/TTTG/TGTTTG`. Continuous hom-alt PL before `Math.round` is
`3517.46536813194416027`, and the integer vector is
`570,148,3485,0,2762,3517`. Injecting the live GKL matrix still moves
the Rust diagnostic by `2.46e-6`. Injecting Java's per-read `T` and
`TGTTTG` allele likelihoods into the Rust genotype calculator reproduces
Java's continuous hom-alt PL within `7e-5`.
`classification: DOWNSTREAM_PER_READ_ALLELE_LIKELIHOOD_MATERIAL`.
`production_change: NONE`.

6R.286: measurement-only. `createAlleleMapper` assigns haplotypes
0–13 to `T` (index 4 is the reference haplotype), 14–18 to `TTTG`,
and 19–23 to `TGTTTG`. Each allele likelihood is the max in that
group after the per-read `best - 4.5` floor. Replaying that map in
the Rust genotype calculator moves continuous hom-alt PL from
`3519.15718129777269496` to `3517.46544109771957665`.
`classification: DOWNSTREAM_HAPLOTYPE_TO_ALLELE_NORMALIZATION_CAUSAL`.
`production_change: NONE`.

6R.287: the same floor-then-max, using Rust's live `createAlleleMapper`
pools, was applied inside `try_genotype_colocated_snp_indel_merge`.
Zero haplotype values changed. Emitted PL stayed `570,0,3518`.
Continuous hom-alt at the emitter is `3517.51626740832807627`.
The allele columns still differ from the Java matrix by up to
`1.11985680133244614` per read.
`classification: DOWNSTREAM_HAPLOTYPE_TO_ALLELE_PRODUCTION_PATH_NOT_JAVA_PL`.
`production_change: NONE`.

6R.288: measurement-only. For read
`HISEQ1:11:H8GV6ADXX:1:2116:18670:99941` flags 99, Rust H0 and H1 are
already `-2.58310416668513199` at the PairHMM return. Normalize and
the merge copy that value. Java's corresponding best is
`-3.70296096801757812`. Substituting Java's 24 haplotype values for
this read leaves continuous hom-alt PL at `3517.51626740832716678`
and integer PL at 3518.
`classification: HAPLOTYPE_LIKELIHOOD_PAIRHMM_RETURN_NOT_CAUSAL`.
`production_change: NONE`.

6R.289: measurement-only. The same read enters PairHMM with different
bases. Java's evidence read is `51H97M` at 29455560 (97 bp). Rust's
scored read is `60H88M` at 29455569 (88 bp). Both are suffixes of the
raw `148M` record at 29455509. Java bases `[9..]` equal the Rust
bases. The 9 bp Java includes and Rust omits are `CAAAGAGGA`.
Base qualities, indel qualities, gap continuation, and haplotype 0
were not the first difference. No counterfactual: the base vectors
have different lengths.
`classification: PAIRHMM_INPUT_READ_BASES`.
`production_change: NONE`.

6R.290: measurement-only. The same read is still `148M` at 29455509 after
`finalizeRegion` on padded span `20:29455460-29455844`. The next clip,
`hardClipToRegion` on the trimmed padded span, produces Java `51H97M`
at 29455560 (span `20:29455560-29455728`) and Rust `60H88M` at
29455569 (span `20:29455569-29455724`). Java keeps raw offsets 51–59
(`CAAAGAGGA`). Rust removes them. No PL counterfactual.
`classification: READ_HARDCLIP_TRIMMED_SPAN_DIVERGENCE`.
`production_change: NONE`.

6R.291: measurement-only. The 13 overlapping trim events match, including
the indel `A>AT` at 29455644. Java `TandemRepeat` sets that event's
padding to `75 + 9 = 84`, so the running span becomes
`29455560-29455728`. Rust keeps indel padding 75 because the one-base
event span does not enter the STR helper, so that step yields
`29455569-29455719`. The shared indel at 29455649 then sets Rust's
right edge to `29455724`. No PL counterfactual.
`classification: TRIM_STR_PADDING_DIVERGENCE`.
`production_change: NONE`.

6R.292: measurement-only. For `A>AT` at 29455644 both sides use reference
window `20:29455460-29455844`. Java drops the anchor and searches from
29455645 (`TTTTTTTT…`) with alt suffix `T`, so the counts are `[8, 9]`
and the STR addition is 9. Rust's helper inspects only
`[29455644, 29455644]` (`A`) and returns `None` because that slice is
shorter than 2. The reference run is `T`, not `A`. No PL
counterfactual.
`classification: STR_REPEAT_SLICE_DIVERGENCE`.
`production_change: NONE`.

6R.293: coordinate-only counterfactual. Replaying the live trim loop with
padding 84 on `A>AT` at 29455644, and every other event unchanged,
yields `29455560-29455728` after that event, after `T>TGTTTG` at
29455649, and at the end of the loop. Production trim is unchanged.
`classification: STR_PADDING_COUNTERFACTUAL_REPRODUCES_JAVA_SPAN`.
`production_change: NONE`.

6R.294: the live path is run twice. The second run forces only the
trimmed padded span to `29455560-29455728`. `trim_modern` still returns
`29455569-29455724`. The diagnostic read becomes `51H97M` at 29455560.
Reference haplotypes become 169 bp on that span. Live continuous
hom-alt PL moves from `3517.51626740832807627` to
`3517.46541826365864836` (Java `3517.46536813194416027`). Integer
hom-alt PL moves from 3518 to 3517. Emitted PL moves from `570,0,3518`
to `570,0,3517`. Production remains `570,0,3518`.
`classification: TRIM_SPAN_CAUSAL_FOR_REMAINING_PL`.
`production_change: NONE`.

6R.295: the Java trim span stays forced. Against the AVX Java capture,
unfloored haplotype cells differ by at most `2.12153372559e-6`
(hap 20 of `HISEQ1:9:H8962ADXX:1:2216:7791:64678`). The largest allele
cell is haplotype 5, allele `T`, delta `2.38641348460e-6`. The sum of
the allele deltas is `5.005406882219e-5`, which is the remaining
continuous-PL residual. Production remains `570,0,3518`.
`classification: RESIDUAL_AT_PAIRHMM_OUTPUT`.
`production_change: NONE`.

6R.296: provenance only. Rust `post_kernel` is 251 reads × 24 raw
haplotype columns. Java `6r285_java_genotyping.tsv` `hap_ll` is 248
reads × columns 19–23, post-normalization, from
`Downstream285.dumpFrozenColumns`. Columns 0–18 are not in that file.
The 73-cell and 610-cell comparisons use only those five columns on
the 122 `allele_ll` reads. The cited haplotype 5 is the last of four
tied Rust T-pool maxima (haplotypes 2–5); Java has no haplotype index
for that allele. The allele sum `0.00005005406882219` is 366
reconstructed Rust pool maxima minus stored `allele_ll` cells.
`classification: CAPTURE_SCOPE_EVIDENCE_GAP`.
`production_change: NONE`.

6R.297: `Capture297` dumps all 24 post-normalization haplotype columns
at `assignGenotypeLikelihoods`. The kernel is `LOGLESS_CACHING`
because AVX does not run on this host. All 122 genotyping reads
match. All 2928 floored Rust cells equal the Java cells. Haplotypes
2, 3, 4, and 5 tie for the diagnostic read's T maximum, and the
same-run allele T equals that maximum. The stored AVX allele T does
not. Columns 19–23 versus the AVX file still have maximum absolute
delta `2.38641348460e-6`. AVX columns 0–18 remain uncaptured.
`classification: LOGLESS_24_COLUMN_MATRIX_IDENTICAL_AVX_0_18_UNCAPTURED`.
`production_change: NONE`.

6R.298: `6r285_java_genotyping.tsv` is produced by `Downstream285`.
`allele_ll` is `dumpEvent` at lines 311–320. The backend is
`AVX_LOGLESS_CACHING` with `useDoublePrecision=false`, not
`FASTEST_AVAILABLE`. Every stored likelihood cell is an exact
binary32 value. Haplotypes 0–13 are absent. On this arm64 host,
constructing that backend throws `Machine does not support AVX PairHMM`
because `libgkl_utils.dylib` is x86_64. No same-input AVX matrix was
produced.
`classification: JAVA_AVX_CAPTURE_NOT_REPRODUCIBLE`.
`production_change: NONE`.

6R.299: the LOGLESS `allele_ll` in `6r297_java_hap_ll.tsv` is 122×3
and equals the per-allele maximum of the 24 haplotype columns.
The span-forced Rust allele matrix matches those 366 cells. The
six genotype likelihoods differ. Integer PL is `570,148,3485,0,2762,3517`
on both sides. Rust continuous hom-alt remains
`3517.46541826365864836`; Java LOGLESS continuous hom-alt is
`3517.46534871093535912`.
`classification: CANONICAL_LOGLESS_DOWNSTREAM_DIVERGENCE`.
`production_change: NONE`.

6R.300: both sides use the 122×3 LOGLESS `allele_ll` matrix.
Homozygote per-read contributions match. Heterozygote per-read
contributions differ. The largest T/TTTG per-read delta is
`-0.00002069939267679` on
`HWI-D00360:8:H88U0ADXX:1:2101:19764:71575` flags 83.
Java uses `MathUtils.approximateLog10SumLog10`; Rust uses exact
`log10_sum_log10`.
`classification: GENOTYPE_LIKELIHOOD_PER_READ_OPERATION_DIVERGENCE`.
`production_change: NONE`.

6R.301: Java `approximateLog10SumLog10(double, double)` adds
`JacobianLogTable.get(larger - smaller)` when that difference is below 8.
Rust production `log10_sum_log10` adds the analytic
`log10(1 + 10^(other - max))`. The existing
`approximate_log10_sum_log10_pair` helper matches the Java result bits
on every heterozygote pair. Replacing only that combine on the
canonical allele matrix reproduces the Java continuous GL and PL
vectors. Production code is unchanged.
`classification: JAVA_JACOBIAN_TABLE_VS_RUST_ANALYTIC_COMBINE`.
`production_change: NONE`.

6R.302: the existing `approximate_log10_sum_log10_pair` helper is
already the production activity-profile implementation of Java
`approximateLog10SumLog10(double, double)`. Calling it for the
heterozygote arm reproduces the canonical Java GL and PL vectors.
Eleven canonical table lookups differ by one ulp from a Rust
`powf`/`log10` reconstruction; those ulps do not change the combined
result. Production code is unchanged.
`classification: JAVA_COMPATIBLE_JACOBIAN_DESIGN_CONFIRMED`.
`production_change: NONE`.

6R.303: the diploid and biallelic heterozygote genotype combines call
`approximate_log10_sum_log10_pair`. On the canonical LOGLESS allele
matrix the six continuous GLs and PLs match Java, and the production
PL conversion emits `570,0,3517`. The live call on the production trim
span still emits `570,0,3518`. `log10_sum_log10` is unchanged.
`classification: JAVA_COMPATIBLE_JACOBIAN_PRODUCTION_FIX_VALIDATED`.

6R.304: the Java trim span `29455560-29455728` composed with the
production Jacobian heterozygote combine reproduces the Java trim
geometry and the shared 122×3 allele cells. Rust also genotypes
`HISEQ1:11:H8GV6ADXX:2:1103:14252:55237` flags 97, whose three allele
likelihoods are 0. Dropping that row and replaying the 122 shared
rows matches the Java GL vector. The live 123-row continuous GL and
PL still differ. Integer PL on this counterfactual is
`570,148,3485,0,2762,3517`. The unforced call remains `570,0,3518`.
`classification: COMPOSED_COUNTERFACTUAL_MATRIX_DIVERGENCE`.
`production_change: NONE`.

6R.305: `HISEQ1:11:H8GV6ADXX:2:1103:14252:55237` flags 97 is paired,
mapped on contig 20, with its mate mapped on contig 7. Java
`filterNonPassingReads` removes it because
`MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE` fails, and the
`calculateGLs` matrix does not contain it. Rust's assembly-region
filter fails the same predicate, then
`score_pairhmm_from_records_java_mate_contig` puts the read back
into the genotyping matrix with likelihood 0 on every haplotype.
`classification: READ_MEMBERSHIP_JAVA_EXCLUDES_RUST_RETAINS`.
`production_change: NONE`.

6R.310: production. Indel trim padding uses
`tandem_repeat_at_event`, the Java anchored tandem-repeat count.
The scan starts one base after the event anchor, and the alternate
allele participates. For `20:29455644 A>AT` the count is 9 and the
padded span is `29455560-29455728`. SNP padding and non-repeat indels
stay on their previous widths. Failed-mate membership and clipped-read
order are unchanged, so the emitted PL is not yet the Java vector.
`classification: STR_PRODUCTION_PATCH_REPRODUCES_JAVA_GEOMETRY`.
`production_change: YES`.

6R.311: production. A paired read whose mapped mate is on another
contig is omitted from PairHMM and genotyping evidence. The existing
`passes_mate_on_same_contig_or_no_mapped_mate` predicate is unchanged.
Length, mapping-quality, and read-group failures still receive a zero
likelihood row. On `20:29455649` the failed-mate read is absent, the
genotyping matrix is the Java 122×3 allele matrix, and clipped-read
order is still the pre-6R.308 order, so the live GL is not yet
bit-identical to Java.
`classification: FAILED_MATE_PRODUCTION_PATCH_REPRODUCES_JAVA_MEMBERSHIP`.
`production_change: YES`.

6R.312: production. Clipped reads that enter PairHMM are ordered with
`java_read_coordinate_compare`, Java's `ReadCoordinateComparator`.
At one clipped start a forward read precedes a reverse read. On
`20:29455649` the 122 genotyping reads match Java's order, the
122×3 allele matrix is unchanged, and the six genotype likelihoods
and continuous PL values match the Java bits. The pre-6R.312
name sort remains only as a diagnostic counterfactual: it still
moves GL index 3 by 2 ULP and index 5 by 1 ULP.
`classification: JAVA_READ_ORDERING_PRODUCTION_PATCH_REPRODUCES_JAVA`.
`production_change: YES`.

6R.314: holdout. The four production semantics are checked on synthetic
inputs other than `20:29455649 T/TGTTTG` and `20:29455644 A>AT`.
Tandem-repeat trim, failed-mate membership, clipped-read order, and
the Jacobian heterozygote combine each pass, and one composite on
`chrHold:4500 G>GCAG` passes those checkpoints in that order.
`classification: STR_HOLDOUT_PASS`, `FAILED_MATE_HOLDOUT_PASS`,
`READ_ORDER_HOLDOUT_PASS`, `JACOBIAN_HOLDOUT_PASS`,
`COMPOSITE_HOLDOUT_PASS`.
`production_change: NONE`.
`HOLDOUT_6R243` was not run.

```text
HOLDOUT_6R314=1 cargo test -p gatk-haplotypecaller --test holdout_6r314_post_milestone_generalization -- --test-threads=1
cargo test -p gatk-haplotypecaller --lib holdout_6r314_failed_mate_zero_row -- --test-threads=1
```

```text
HOLDOUT_6R158=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
  cargo test -p gatk-haplotypecaller --test holdout_6r158_class_a3_preserve -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r158_class_a3_preserves_calculator_genotype -- --test-threads=1
HOLDOUT_6R159=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
  cargo test -p gatk-haplotypecaller --test holdout_6r159_homref_emission -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r159_homref_emission_boundary_contract -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r160_next_vcf_divergence_inventory -- --test-threads=1
HOLDOUT_6R160=1 cargo test -p gatk-haplotypecaller --test holdout_6r160_next_vcf_divergence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r161_informative_votes -- --test-threads=1
HOLDOUT_6R161=1 cargo test -p gatk-haplotypecaller --test holdout_6r161_informative_votes -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r162_extra_deletion_haplotype -- --test-threads=1
HOLDOUT_6R162=1 cargo test -p gatk-haplotypecaller --test holdout_6r162_extra_deletion_haplotype -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r163_dangling_merge_parity -- --test-threads=1
HOLDOUT_6R163=1 cargo test -p gatk-haplotypecaller --test holdout_6r163_dangling_merge -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r164_no_dangling_merge_haplotype_materialization -- --test-threads=1
HOLDOUT_6R164=1 cargo test -p gatk-haplotypecaller --test holdout_6r164_no_dangling_merge -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r165_info_annotation_divergence -- --test-threads=1
HOLDOUT_6R165=1 cargo test -p gatk-haplotypecaller --test holdout_6r165_info_annotation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r166_fs_sor_java_evidence_source -- --test-threads=1
HOLDOUT_6R166=1 cargo test -p gatk-haplotypecaller --test holdout_6r166_info_annotation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r167_mq_java_evidence_source -- --test-threads=1
HOLDOUT_6R167=1 cargo test -p gatk-haplotypecaller --test holdout_6r167_mq_evidence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r168_mq_rms_formula -- --test-threads=1
HOLDOUT_6R168=1 cargo test -p gatk-haplotypecaller --test holdout_6r168_mq_rms -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r169_live_info_inventory -- --test-threads=1
HOLDOUT_6R169=1 cargo test -p gatk-haplotypecaller --test holdout_6r169_live_info -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r170_read_pos_rank_sum_java_evidence_source -- --test-threads=1
HOLDOUT_6R170=1 cargo test -p gatk-haplotypecaller --test holdout_6r170_read_pos_rank_sum -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r171_read_pos_rank_sum_empty_result_contract -- --test-threads=1
HOLDOUT_6R171=1 cargo test -p gatk-haplotypecaller --test holdout_6r171_read_pos_rank_sum_empty -- --test-threads=1
cargo test -p gatk-haplotypecaller --lib forensic_6r172_semantic_matrix_production_fn -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r172_read_pos_rank_sum_undefined_vs_zero -- --test-threads=1
HOLDOUT_6R172=1 cargo test -p gatk-haplotypecaller --test holdout_6r172_read_pos_rank_sum_undefined -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r173_live_vcf_reconnaissance -- --test-threads=1
HOLDOUT_6R173=1 cargo test -p gatk-haplotypecaller --test holdout_6r173_live_vcf_reconnaissance -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r174_info_dp_java_coverage_source -- --test-threads=1
HOLDOUT_6R174=1 cargo test -p gatk-haplotypecaller --test holdout_6r174_info_dp -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r175_sor_java_evidence_source -- --test-threads=1
HOLDOUT_6R175=1 cargo test -p gatk-haplotypecaller --test holdout_6r175_sor -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r176_pairhmm_mate_contig_membership -- --test-threads=1
HOLDOUT_6R176=1 cargo test -p gatk-haplotypecaller --test holdout_6r176_sor -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r177_inbreeding_coeff_emit_predicate -- --test-threads=1
HOLDOUT_6R177=1 cargo test -p gatk-haplotypecaller --test holdout_6r177_inbreeding -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r178_inbreeding_coeff_min_samples -- --test-threads=1
HOLDOUT_6R178=1 cargo test -p gatk-haplotypecaller --test holdout_6r178_inbreeding -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r179_indel_info_membership -- --test-threads=1
HOLDOUT_6R179=1 cargo test -p gatk-haplotypecaller --test holdout_6r179_indel_info -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r180_indel_info_annotation_evidence -- --test-threads=1
HOLDOUT_6R180=1 cargo test -p gatk-haplotypecaller --test holdout_6r180_indel_info -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r181_cluster_tg_annotation_membership -- --test-threads=1
HOLDOUT_6R181=1 cargo test -p gatk-haplotypecaller --test holdout_6r181_cluster_tg -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r182_per_variant_annotation_object -- --test-threads=1
HOLDOUT_6R182=1 cargo test -p gatk-haplotypecaller --test holdout_6r182_per_variant_object -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r183_poorly_modeled_stored_membership -- --test-threads=1
HOLDOUT_6R183=1 cargo test -p gatk-haplotypecaller --test holdout_6r183_poorly_modeled -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r184_post_refresh_refilter_lifecycle -- --test-threads=1
HOLDOUT_6R184=1 cargo test -p gatk-haplotypecaller --test holdout_6r184_post_refresh_refilter -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r185_post_refresh_java_order_filter -- --test-threads=1
HOLDOUT_6R185=1 cargo test -p gatk-haplotypecaller --test holdout_6r185_post_refresh_filter -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r186_cluster_tg_annotation_object -- --test-threads=1
HOLDOUT_6R186=1 cargo test -p gatk-haplotypecaller --test holdout_6r186_cluster_tg_annotation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r187_fresh_info_reconnaissance -- --test-threads=1
HOLDOUT_6R187=1 cargo test -p gatk-haplotypecaller --test holdout_6r187_fresh_info -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r188_wrong_annotation_source_object -- --test-threads=1
HOLDOUT_6R188=1 cargo test -p gatk-haplotypecaller --test holdout_6r188_wrong_annotation_source -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r189_site_annotation_uses_stored_haplotype_object -- --test-threads=1
HOLDOUT_6R189=1 cargo test -p gatk-haplotypecaller --test holdout_6r189_site_annotation_stored_hap -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r190_fresh_info_reconnaissance -- --test-threads=1
HOLDOUT_6R190=1 cargo test -p gatk-haplotypecaller --test holdout_6r190_fresh_info -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r191_read_pos_rank_sum_uses_annotation_likelihoods -- --test-threads=1
HOLDOUT_6R191=1 cargo test -p gatk-haplotypecaller --test holdout_6r191_read_pos_rank_sum -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r192_fresh_info_reconnaissance -- --test-threads=1
HOLDOUT_6R192=1 cargo test -p gatk-haplotypecaller --test holdout_6r192_fresh_info -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r193_baseq_rank_sum_emission -- --test-threads=1
HOLDOUT_6R193=1 cargo test -p gatk-haplotypecaller --test holdout_6r193_baseq_rank_sum -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r194_mq_rank_sum_emission -- --test-threads=1
HOLDOUT_6R194=1 cargo test -p gatk-haplotypecaller --test holdout_6r194_mq_rank_sum -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r195_qd_follows_qual -- --test-threads=1
HOLDOUT_6R195=1 cargo test -p gatk-haplotypecaller --test holdout_6r195_qd_follows_qual -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r196_qual_indel_af_prior -- --test-threads=1
HOLDOUT_6R196=1 cargo test -p gatk-haplotypecaller --test holdout_6r196_qual_indel_af_prior -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r197_fresh_info_reconnaissance -- --test-threads=1
HOLDOUT_6R197=1 cargo test -p gatk-haplotypecaller --test holdout_6r197_fresh_info -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r198_per_variant_annotation_object -- --test-threads=1
HOLDOUT_6R198=1 cargo test -p gatk-haplotypecaller --test holdout_6r198_per_variant_object -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r199_cluster_downstream_annotation_object -- --test-threads=1
HOLDOUT_6R199=1 cargo test -p gatk-haplotypecaller --test holdout_6r199_cluster_downstream_annotation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r200_coverage_no_secondary_overlap_filter -- --test-threads=1
HOLDOUT_6R200=1 cargo test -p gatk-haplotypecaller --test holdout_6r200_coverage_secondary_overlap -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r201_coverage_consumes_annotation_cardinality -- --test-threads=1
HOLDOUT_6R201=1 cargo test -p gatk-haplotypecaller --test holdout_6r201_coverage_cardinality -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r202_ranksum_annotation_boundary -- --test-threads=1
HOLDOUT_6R202=1 cargo test -p gatk-haplotypecaller --test holdout_6r202_ranksum_boundary -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r203_ranksum_get_element_java_cigar_semantics -- --test-threads=1
HOLDOUT_6R203=1 cargo test -p gatk-haplotypecaller --test holdout_6r203_ranksum_get_element -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r204_coverage_evidence_boundary -- --test-threads=1
HOLDOUT_6R204=1 cargo test -p gatk-haplotypecaller --test holdout_6r204_coverage_evidence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r205_two_read_hom_alt_annotation_object -- --test-threads=1
HOLDOUT_6R205=1 cargo test -p gatk-haplotypecaller --test holdout_6r205_two_read_hom_alt_annotation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r206_annotation_object_membership -- --test-threads=1
HOLDOUT_6R206=1 cargo test -p gatk-haplotypecaller --test holdout_6r206_annotation_object_membership -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r207_mid_a_one_read_hom_alt_annotation_object -- --test-threads=1
HOLDOUT_6R207=1 cargo test -p gatk-haplotypecaller --test holdout_6r207_mid_a_one_read_hom_alt_annotation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r208_mid_b_annotation_evidence -- --test-threads=1
HOLDOUT_6R208=1 cargo test -p gatk-haplotypecaller --test holdout_6r208_mid_b_annotation_evidence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r209_preserve_loc_loop_mate_evidence -- --test-threads=1
HOLDOUT_6R209=1 cargo test -p gatk-haplotypecaller --test holdout_6r209_preserve_loc_loop_mate_evidence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r210_annotation_source_boundary -- --test-threads=1
HOLDOUT_6R210=1 cargo test -p gatk-haplotypecaller --test holdout_6r210_annotation_source_boundary -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r211_annotation_membership_boundary -- --test-threads=1
HOLDOUT_6R211=1 cargo test -p gatk-haplotypecaller --test holdout_6r211_annotation_membership_boundary -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r212_event_weak_sparse_annotation_boundary -- --test-threads=1
HOLDOUT_6R212=1 cargo test -p gatk-haplotypecaller --test holdout_6r212_event_weak_sparse_annotation_boundary -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r213_first_genotype_qual_boundary -- --test-threads=1
HOLDOUT_6R213=1 cargo test -p gatk-haplotypecaller --test holdout_6r213_first_genotype_qual_boundary -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r214_marginalized_allele_likelihood_boundary -- --test-threads=1
HOLDOUT_6R214=1 cargo test -p gatk-haplotypecaller --test holdout_6r214_marginalized_allele_likelihood -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r215_haplotype_population_boundary -- --test-threads=1
HOLDOUT_6R215=1 cargo test -p gatk-haplotypecaller --test holdout_6r215_haplotype_population -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r216_rt_kbest_java_restriction -- --test-threads=1
HOLDOUT_6R216=1 cargo test -p gatk-haplotypecaller --test holdout_6r216_rt_kbest_java_restriction -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r217_rt_cycle_abort_scope -- --test-threads=1
HOLDOUT_6R217=1 cargo test -p gatk-haplotypecaller --test holdout_6r217_rt_cycle_abort_scope -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r218_java_seqgraph_cycle_gate -- --test-threads=1
HOLDOUT_6R218=1 cargo test -p gatk-haplotypecaller --test holdout_6r218_java_seqgraph_cycle_gate -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r219_ad_membership_boundary -- --test-threads=1
HOLDOUT_6R219=1 cargo test -p gatk-haplotypecaller --test holdout_6r219_ad_membership -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r220_format_ad_subset_boundary -- --test-threads=1
HOLDOUT_6R220=1 cargo test -p gatk-haplotypecaller --test holdout_6r220_format_ad_subset -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r221_call_region_assign_boundary -- --test-threads=1
HOLDOUT_6R221=1 cargo test -p gatk-haplotypecaller --test holdout_6r221_call_region_assign -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r222_loc_loop_likelihood_lifecycle -- --test-threads=1
HOLDOUT_6R222=1 cargo test -p gatk-haplotypecaller --test holdout_6r222_loc_loop_likelihood -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r223_first_inner_likelihood_transition -- --test-threads=1
HOLDOUT_6R223=1 cargo test -p gatk-haplotypecaller --test holdout_6r223_first_inner_likelihood -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r224_tls_cache_key_identity -- --test-threads=1
HOLDOUT_6R224=1 cargo test -p gatk-haplotypecaller --test holdout_6r224_tls_cache_key -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r225_cache_semantic_contract -- --test-threads=1
HOLDOUT_6R225=1 cargo test -p gatk-haplotypecaller --test holdout_6r225_cache_semantic_contract -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r226_cache_key_uses_sparse_population_identity -- --test-threads=1
HOLDOUT_6R226=1 cargo test -p gatk-haplotypecaller --test holdout_6r226_cache_key_sparse_identity -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r227_post_cache_genotype_likelihood_boundary -- --test-threads=1
HOLDOUT_6R227=1 cargo test -p gatk-haplotypecaller --test holdout_6r227_post_cache_genotype_likelihood -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r228_java_rust_retain_evidence_membership -- --test-threads=1
HOLDOUT_6R228=1 cargo test -p gatk-haplotypecaller --test holdout_6r228_java_rust_retain_evidence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r229_java_hap_ll_membership_pipeline -- --test-threads=1
HOLDOUT_6R229=1 cargo test -p gatk-haplotypecaller --test holdout_6r229_java_hap_ll_membership -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r230_five_read_likelihood_boundary -- --test-threads=1
HOLDOUT_6R230=1 cargo test -p gatk-haplotypecaller --test holdout_6r230_five_read_likelihood -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r231_read_clipping_boundary -- --test-threads=1
HOLDOUT_6R231=1 cargo test -p gatk-haplotypecaller --test holdout_6r231_read_clipping -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r232_missing_variation_event_boundary -- --test-threads=1
HOLDOUT_6R232=1 cargo test -p gatk-haplotypecaller --test holdout_6r232_missing_variation_event -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r233_c_carrier_haplotype_boundary -- --test-threads=1
HOLDOUT_6R233=1 cargo test -p gatk-haplotypecaller --test holdout_6r233_c_carrier_haplotype -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r234_c_carrier_kbest_rank_boundary -- --test-threads=1
HOLDOUT_6R234=1 cargo test -p gatk-haplotypecaller --test holdout_6r234_c_carrier_kbest_rank -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r235_carrier_sink_edge_split -- --test-threads=1
HOLDOUT_6R235=1 cargo test -p gatk-haplotypecaller --test holdout_6r235_carrier_sink_edge -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r236_seqgraph_sink_topology -- --test-threads=1
HOLDOUT_6R236=1 cargo test -p gatk-haplotypecaller --test holdout_6r236_seqgraph_sink_topology -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r237_dangling_tail_splice_decision -- --test-threads=1
HOLDOUT_6R237=1 cargo test -p gatk-haplotypecaller --test holdout_6r237_dangling_tail_splice -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r238_dangling_index_sentinel_contract -- --test-threads=1
HOLDOUT_6R238=1 cargo test -p gatk-haplotypecaller --test holdout_6r238_dangling_index_sentinel -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r239_dangling_index_sentinel_java_semantics -- --test-threads=1
HOLDOUT_6R239=1 cargo test -p gatk-haplotypecaller --test holdout_6r239_dangling_index_sentinel -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r240_covering_gc_emit -- --test-threads=1
HOLDOUT_6R240=1 cargo test -p gatk-haplotypecaller --test holdout_6r240_covering_gc_emit -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r241_strict_java_emit_confidence -- --test-threads=1
HOLDOUT_6R241=1 cargo test -p gatk-haplotypecaller --test holdout_6r241_strict_java_emit_confidence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r242_fresh_parity_inventory -- --test-threads=1
HOLDOUT_6R242=1 cargo test -p gatk-haplotypecaller --test holdout_6r242_fresh_parity_inventory -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r243_annotation_likelihood_lifecycle -- --test-threads=1
HOLDOUT_6R243=1 cargo test -p gatk-haplotypecaller --test holdout_6r243_annotation_likelihood_lifecycle -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r244_colocated_merge_annotation_lifecycle -- --test-threads=1
HOLDOUT_6R244=1 cargo test -p gatk-haplotypecaller --test holdout_6r244_colocated_merge_annotation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r245_colocated_merge_annotation_propagation -- --test-threads=1
HOLDOUT_6R245=1 cargo test -p gatk-haplotypecaller --test holdout_6r245_colocated_merge_annotation_propagation -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r246_colocated_merge_annotation_attach -- --test-threads=1
HOLDOUT_6R246=1 cargo test -p gatk-haplotypecaller --test holdout_6r246_colocated_merge_annotation_attach -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r247_first_pl_divergence -- --test-threads=1
HOLDOUT_6R247=1 cargo test -p gatk-haplotypecaller --test holdout_6r247_first_pl_divergence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r248_merged_hom_alt_gl -- --test-threads=1
HOLDOUT_6R248=1 cargo test -p gatk-haplotypecaller --test holdout_6r248_merged_hom_alt_gl -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r249_tgtttg_allele_row -- --test-threads=1
HOLDOUT_6R249=1 cargo test -p gatk-haplotypecaller --test holdout_6r249_tgtttg_allele_row -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r250_tgtttg_pairhmm_columns -- --test-threads=1
HOLDOUT_6R250=1 cargo test -p gatk-haplotypecaller --test holdout_6r250_tgtttg_pairhmm_columns -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r251_haplotype_sequence_identity -- --test-threads=1
HOLDOUT_6R251=1 cargo test -p gatk-haplotypecaller --test holdout_6r251_haplotype_sequence_identity -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r252_java_actual_haplotype_comparison -- --test-threads=1
HOLDOUT_6R252=1 cargo test -p gatk-haplotypecaller --test holdout_6r252_java_actual_haplotype_comparison -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r253_java_vs_rust_trim_span -- --test-threads=1
HOLDOUT_6R253=1 cargo test -p gatk-haplotypecaller --test holdout_6r253_java_vs_rust_trim_span -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r254_java_trim_bases_pairhmm_consequence -- --test-threads=1
HOLDOUT_6R254=1 cargo test -p gatk-haplotypecaller --test holdout_6r254_java_trim_bases_pairhmm_consequence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r255_read_clip_interval_pairhmm_consequence -- --test-threads=1
HOLDOUT_6R255=1 cargo test -p gatk-haplotypecaller --test holdout_6r255_read_clip_interval_pairhmm_consequence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r256_joint_java_input_plane -- --test-threads=1
HOLDOUT_6R256=1 cargo test -p gatk-haplotypecaller --test holdout_6r256_joint_java_input_plane -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r257_pairhmm_input_plane -- --test-threads=1
HOLDOUT_6R257=1 cargo test -p gatk-haplotypecaller --test holdout_6r257_pairhmm_input_plane -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r258_pairhmm_recurrence -- --test-threads=1
HOLDOUT_6R258=1 cargo test -p gatk-haplotypecaller --test holdout_6r258_pairhmm_recurrence -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r259_gkl_powf_ph2pr -- --test-threads=1
HOLDOUT_6R259=1 cargo test -p gatk-haplotypecaller --test holdout_6r259_gkl_powf_ph2pr -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r260_gkl_powf_ph2pr -- --test-threads=1
HOLDOUT_6R260=1 cargo test -p gatk-haplotypecaller --test holdout_6r260_gkl_powf_ph2pr -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r261_gkl_float_one_minus_ph2pr -- --test-threads=1
HOLDOUT_6R261=1 cargo test -p gatk-haplotypecaller --test holdout_6r261_gkl_float_one_minus_ph2pr -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r262_gkl_float_ph2pr_div3 -- --test-threads=1
HOLDOUT_6R262=1 cargo test -p gatk-haplotypecaller --test holdout_6r262_gkl_float_ph2pr_div3 -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r263_gkl_avx_distm_blend -- --test-threads=1
HOLDOUT_6R263=1 cargo test -p gatk-haplotypecaller --test holdout_6r263_gkl_avx_distm_blend -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r264_gkl_avx_m_update_distm_mul -- --test-threads=1
HOLDOUT_6R264=1 cargo test -p gatk-haplotypecaller --test holdout_6r264_gkl_avx_m_update_distm_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r265_gkl_avx_x_update -- --test-threads=1
HOLDOUT_6R265=1 cargo test -p gatk-haplotypecaller --test holdout_6r265_gkl_avx_x_update -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r266_gkl_avx_x_update_xx_mul -- --test-threads=1
HOLDOUT_6R266=1 cargo test -p gatk-haplotypecaller --test holdout_6r266_gkl_avx_x_update_xx_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r267_gkl_avx_x_update_add -- --test-threads=1
HOLDOUT_6R267=1 cargo test -p gatk-haplotypecaller --test holdout_6r267_gkl_avx_x_update_add -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r268_gkl_avx_y_update_my_mul -- --test-threads=1
HOLDOUT_6R268=1 cargo test -p gatk-haplotypecaller --test holdout_6r268_gkl_avx_y_update_my_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r269_gkl_avx_y_update_yy_mul -- --test-threads=1
HOLDOUT_6R269=1 cargo test -p gatk-haplotypecaller --test holdout_6r269_gkl_avx_y_update_yy_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r270_gkl_avx_y_update_add -- --test-threads=1
HOLDOUT_6R270=1 cargo test -p gatk-haplotypecaller --test holdout_6r270_gkl_avx_y_update_add -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r271_gkl_avx_m_update_mm_mul -- --test-threads=1
HOLDOUT_6R271=1 cargo test -p gatk-haplotypecaller --test holdout_6r271_gkl_avx_m_update_mm_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r272_gkl_avx_m_update_xgapm_mul -- --test-threads=1
HOLDOUT_6R272=1 cargo test -p gatk-haplotypecaller --test holdout_6r272_gkl_avx_m_update_xgapm_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r273_gkl_avx_m_update_ygapm_mul -- --test-threads=1
HOLDOUT_6R273=1 cargo test -p gatk-haplotypecaller --test holdout_6r273_gkl_avx_m_update_ygapm_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r274_gkl_avx_m_update_mm_xgapm_add -- --test-threads=1
HOLDOUT_6R274=1 cargo test -p gatk-haplotypecaller --test holdout_6r274_gkl_avx_m_update_mm_xgapm_add -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r275_gkl_avx_m_update_mxgap_ygapm_add -- --test-threads=1
HOLDOUT_6R275=1 cargo test -p gatk-haplotypecaller --test holdout_6r275_gkl_avx_m_update_mxgap_ygapm_add -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r276_gkl_avx_m_update_sum_distm_mul -- --test-threads=1
HOLDOUT_6R276=1 cargo test -p gatk-haplotypecaller --test holdout_6r276_gkl_avx_m_update_sum_distm_mul -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r277_gkl_avx_last_stripe_summ_add -- --test-threads=1
HOLDOUT_6R277=1 cargo test -p gatk-haplotypecaller --test holdout_6r277_gkl_avx_last_stripe_summ_add -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r278_gkl_avx_last_stripe_sumx_add -- --test-threads=1
HOLDOUT_6R278=1 cargo test -p gatk-haplotypecaller --test holdout_6r278_gkl_avx_last_stripe_sumx_add -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r279_gkl_avx_last_stripe_summx_add -- --test-threads=1
HOLDOUT_6R279=1 cargo test -p gatk-haplotypecaller --test holdout_6r279_gkl_avx_last_stripe_summx_add -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r280_gkl_avx_last_stripe_summx_lane_extract -- --test-threads=1
HOLDOUT_6R280=1 cargo test -p gatk-haplotypecaller --test holdout_6r280_gkl_avx_last_stripe_summx_lane_extract -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r281_gkl_result_float_min_accepted -- --test-threads=1
HOLDOUT_6R281=1 cargo test -p gatk-haplotypecaller --test holdout_6r281_gkl_result_float_min_accepted -- --test-threads=1
cargo test -p gatk-haplotypecaller --test forensic_6r282_gkl_result_float_log10f -- --test-threads=1
HOLDOUT_6R282=1 cargo test -p gatk-haplotypecaller --test holdout_6r282_gkl_result_float_log10f -- --test-threads=1
```

## How parity is established

Algorithmic equivalence is **not** “the Rust file looks like the Java class.” It is:

1. Identify the **Java 4.4 contract** from the pinned source (and live Docker when the
   executable can show the observable).
2. Build an **observable test** (graph dump, EventMap, VCF field, RNG prefix).
3. Locate the **first proven divergence** — stop stacking later symptoms.
4. Apply the **smallest general fix** (no locus hard-codes, no widening of P12 bands).
5. Re-run independent regression / holdout tests (`six_r*`, lib suite,
   `p12_call_none_mid_b_test`).

Rust-native modules are preferred over cloning Java type trees. Details:
[`PARITY_MILESTONE_6R.md`](PARITY_MILESTONE_6R.md).

## Classes of Java contracts found on mid-B

These are **general** contracts, demonstrated on mid-B, not coordinate special cases:

| Contract | Java (4.4) | What was wrong |
|----------|------------|----------------|
| Graph uniqueness / k-mer ladder | Uniqueness on `refHaplotype.getBases()` (extended region), not the ±500 padded assembly REF | Rust skipped k=25 on a padded window that is non-unique |
| Dangling-head mismatch cap | `best_prefix_match` hard-aborts when the mismatch budget overflows | Rust recovered heads Java would reject |
| `getBasesForPath` | Expand every `inDegree==0` source when `expandSource` is true | Rust dropped some sources |
| Allele keep | Well-supported SNPs kept even if two haplotypes share them (default HC does not run the unique-supporter collapse) | Shared SNPs dropped |
| Trimmer `maxEnd` | `max(maxEnd, vc.getEnd()+padding)` | Rust accumulated padding per event (94M vs Java 54M) |
| AF EM loop | Dirichlet update **every** iteration, including the last, then P(no variant) | Rust broke before the last update (QUAL 78.583 vs 78.32) |
| MLEAC / MLEAF | `round(EM expected alt count)` / `MLEAC / AN`; **not** the called GT | Rust copied 1/1 → MLEAC=2 |
| QualByDepth | If raw QD ≥ 35: `30 + Random(47382911).nextGaussian()*3` | Rust capped at 30 with no jitter |
| QualByDepth RNG lifetime | One JVM-static `Utils.randomGenerator`; draw iff raw QD ≥ 35; order = sequential walker emit | Thread-local RNG + jitter inside Rayon region emit reseeds per worker (716-cluster restarted at 25.36) |
| `findBestPaths` retention | Keep if SW CIGAR ref-span equals the reference haplotype CIGAR ref-span (≥ 30, no `N`). Sequence length is not a gate. | Rust dropped alts with `len < 75%` of padded ref (`28M171D160M` / 188 bp vs 359) |
| EventMap CIGAR alleles | `processCigarForInitialEvents` emits D/I with no allele-length cap (regular bases; skip unresolved edge I) | Rust dropped CIGAR events with `REF/ALT.len() > 40` (`171D` REF=172) |
| EventMap vs padded REF | EventMap uses the haplotype CIGAR only; `trimTo` does not re-SW equal-length SNP haps against the untrimmed pad | Supplemental Indel SW vs pad invented a spanning D; `prefer_dominant_spanning_indels` then dropped SNPs Java still emits |
| EventMap union | `getAllVariantContexts` keeps every per-haplotype EventMap allele; no nested-SNP drop inside another hap’s spanning indel | `collect_variation_events` / EventMap regen applied `prefer_dominant_spanning_indels` (Rust-only) |

## What remains unknown / out of scope

- Genome-wide, autosome, or multi-sample HC equivalence.
- Sharing one process-global `Random` with reservoir downsampling on intervals that
  overflow `max-reads-per-alignment-start` **before** the first high-QD site (mid-B
  has two reads; streams coincide).
- VCF QUAL print rounding of 78.323 → 78.32.
- Waivers still in force: **W-H1**, **W-H3**, and other claim-matrix scoped rows.

## Tests

Reusable gates (not a forensic diary):

```text
cargo test -p gatk-haplotypecaller --lib -- --test-threads=1
cargo test -p gatk-haplotypecaller --test p12_call_none_mid_b_test
HOLDOUT_6R43=1 cargo test -p gatk-haplotypecaller --test holdout_6r43_test
```

`six_r*` tests under `gatk-haplotypecaller` pin the mid-B contracts above without
requiring the 6R markdown reports. Chr20_tiny genotype-entry holdouts are env-gated
(`HOLDOUT_6R130`…`HOLDOUT_6R262`); production k-best is unchanged unless that env is set.

Independent-region discovery (not whole-codebase parity; 6R.43 snapshot):
[`parity/6R.43_HOLDOUT_MATRIX.md`](parity/6R.43_HOLDOUT_MATRIX.md).
Retention contract vs `p12_snp_cluster`: [`parity/6R.45_RETENTION.md`](parity/6R.45_RETENTION.md).
Trim vs missing 171D EventMap allele: [`parity/6R.46_TRIM.md`](parity/6R.46_TRIM.md).
EventMap CIGAR allele-length (no 40 bp cap): [`parity/6R.47_EVENTMAP.md`](parity/6R.47_EVENTMAP.md).
QualByDepth process-global RNG stream: [`parity/6R.48_QD_RNG.md`](parity/6R.48_QD_RNG.md).
6R.49: skip EventMap supplemental indel SW when hap length equals the trimmed reference haplotype.
6R.50: EventMap union does not apply `prefer_dominant_spanning_indels`.

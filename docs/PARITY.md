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

## Independent chr20_tiny genotype-boundary holdouts (6R.130–6R.221)

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
(`HOLDOUT_6R130`…`HOLDOUT_6R221`); production k-best is unchanged unless that env is set.

Independent-region discovery (not whole-codebase parity; 6R.43 snapshot):
[`parity/6R.43_HOLDOUT_MATRIX.md`](parity/6R.43_HOLDOUT_MATRIX.md).
Retention contract vs `p12_snp_cluster`: [`parity/6R.45_RETENTION.md`](parity/6R.45_RETENTION.md).
Trim vs missing 171D EventMap allele: [`parity/6R.46_TRIM.md`](parity/6R.46_TRIM.md).
EventMap CIGAR allele-length (no 40 bp cap): [`parity/6R.47_EVENTMAP.md`](parity/6R.47_EVENTMAP.md).
QualByDepth process-global RNG stream: [`parity/6R.48_QD_RNG.md`](parity/6R.48_QD_RNG.md).
6R.49: skip EventMap supplemental indel SW when hap length equals the trimmed reference haplotype.
6R.50: EventMap union does not apply `prefer_dominant_spanning_indels`.

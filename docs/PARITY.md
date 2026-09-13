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

## Independent chr20_tiny genotype-boundary holdouts (6R.130–6R.185)

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
(`HOLDOUT_6R130`…`HOLDOUT_6R185`); production k-best is unchanged unless that env is set.

Independent-region discovery (not whole-codebase parity; 6R.43 snapshot):
[`parity/6R.43_HOLDOUT_MATRIX.md`](parity/6R.43_HOLDOUT_MATRIX.md).
Retention contract vs `p12_snp_cluster`: [`parity/6R.45_RETENTION.md`](parity/6R.45_RETENTION.md).
Trim vs missing 171D EventMap allele: [`parity/6R.46_TRIM.md`](parity/6R.46_TRIM.md).
EventMap CIGAR allele-length (no 40 bp cap): [`parity/6R.47_EVENTMAP.md`](parity/6R.47_EVENTMAP.md).
QualByDepth process-global RNG stream: [`parity/6R.48_QD_RNG.md`](parity/6R.48_QD_RNG.md).
6R.49: skip EventMap supplemental indel SW when hap length equals the trimmed reference haplotype.
6R.50: EventMap union does not apply `prefer_dominant_spanning_indels`.

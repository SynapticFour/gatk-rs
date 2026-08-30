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

## Test status (6R.41 / 6R.42 hygiene)

After this documentation cleanup: `gatk-haplotypecaller --lib` **576** passed
(`--test-threads=1`); `p12_call_none_mid_b_test` **1** passed;
`cargo fmt --all -- --check` ok. Dangling-recovery unit tests were moved out of
the production module so the N-3 size gate still holds; algorithm unchanged.

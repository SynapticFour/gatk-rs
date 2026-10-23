# 6R milestone — canonical mid-B HaplotypeCaller path

**Date:** 2026-08-30
**Pinned Java:** GATK **4.4.0.0**, SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`
**Image:** `us.gcr.io/broad-gatk/gatk:4.4.0.0` ([`GATK_PINNED.env`](GATK_PINNED.env))
**Interval:** `2:92317000-92319000` (ActiveFull `2:92317262–92317491`)

This is the short public note for that milestone. The site-by-site log is
the engineering-record half of [`PARITY.md`](PARITY.md). Later 6R and 7R
notes are indexed from [`parity/README.md`](parity/README.md). This page is
not a claim-matrix Yes row beyond the scoped mid-B row already in
[`CLAIM_MATRIX.md`](CLAIM_MATRIX.md).

## What was proven on this region

Rust matched GATK 4.4 on the three oracle SNPs:

`92317399 C/A`, `92317407 T/C`, `92317412 G/C`

with FORMAT `GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0`, MLEAC=1, MLEAF=0.500,
QUAL 78.32 (printed), QD 25.36 / 28.73 / 30.97.

**Canonical mid-B: CONVERGED on that region. Whole-codebase GATK 4.4
parity: NOT ESTABLISHED.**

## What was fixed (classes, not a diary)

Java-faithful contracts, each with an observable test: reference-haplotype
k-mer uniqueness boundary; dangling-head mismatch-cap abort;
`getBasesForPath` any-source expansion; allele-keep (no unique-supporter
collapse on default HC); trimmer `maxEnd`; AF EM last Dirichlet update;
MLEAC from EM not called GT; QualByDepth `Random(47382911)` jitter.

The method was: Java contract, observable test, first divergence, smallest
general fix, regression. The narrative of later sites (6R.43 through the
`20:29455649 T/TGTTTG` witness and 6R.314 holdouts) stays in
[`PARITY.md`](PARITY.md), [`parity/6R.312_FINAL_CLOSURE.md`](parity/6R.312_FINAL_CLOSURE.md),
and [`parity/6R.314_HOLDOUT.md`](parity/6R.314_HOLDOUT.md). Those files are
the evidence. They are not restated here.

One correction against this milestone’s original closing line: production
SeqGraph k-best is no longer `legacy_1024` when the policy variable is
unset. See [`parity/KBEST_RESOURCE_POLICY.md`](parity/KBEST_RESOURCE_POLICY.md)
(7R.22). That change is not an additional Yes row.

## What this milestone left unknown

Genome-wide or autosome HC equivalence, other samples, intervals where the
QD random stream interleaves differently, and the claim-matrix waivers
(W-H1, W-H3, and the rest named there).

## Recorded test note from the 6R.41 / 6R.42 cleanup

These counts were written with the 2026-08-30 documentation pass. They are
not a current workspace result:

- `gatk-haplotypecaller --lib`: **576** passed (`--test-threads=1`)
- `p12_call_none_mid_b_test`: **1** passed
- `cargo fmt --all -- --check`: ok

Dangling-recovery unit tests were moved out of the production module so the
N-3 size gate still held. The algorithm was unchanged by that move.

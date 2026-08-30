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
```

`six_r*` tests under `gatk-haplotypecaller` pin the contracts above without requiring
the 6R markdown reports.

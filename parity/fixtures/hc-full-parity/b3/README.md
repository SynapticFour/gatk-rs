# B.3 — `apply` count vs `callRegion` fast path

## Semantics (GATK 4)

- **`AssemblyRegionWalker`**: one `apply(AssemblyRegion, …)` per region emitted by `AssemblyRegionIterator` (active **and** inactive).
- **`HaplotypeCallerEngine.callRegion`**: if `!region.isActive()`, returns immediately via `referenceModelForNoVariation` — **no** `assembleReads` / graph work (**inactive fast path**).

Rust models this as [`WalkerApplyStats`](../../../../gatk-haplotypecaller/src/walker_apply.rs): `total_apply` equals region count; `inactive_fast_path` counts regions where [`call_disposition`](../../../../gatk-haplotypecaller/src/walker_apply.rs) is `InactiveReferenceFastPath`.

## Artifact (`*.tsv`)

Two lines:

1. Header: `total_apply`, `inactive_fast_path`, `active_full` (tab-separated)
2. One data row; invariant `total_apply == inactive_fast_path + active_full`

## Cases manifest

[`cases.tsv`](./cases.tsv) adds **`padding`** vs B.2: shard padding as `u64`, or `-` for GATK default (**100**). Use a smaller padding on tiny references when disjoint `-L` segments must stay in **separate** padded spans (multiple `apply` steps).

## Regenerating expected

Same pattern as [`../b2/README.md`](../b2/README.md); use `PARITY_CARGO_TARGET_DIR` when needed:

```bash
export PARITY_CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-$PWD/target}"
cargo run -p gatk-haplotypecaller --example hc_full_parity_gate_dump -- \
  apply-summary parity/fixtures/reference.fa parity/fixtures/sample.bam chr1:5-15
# Non-default padding (see `chr1_disjoint_pad5_two_inactive` in `cases.tsv`):
cargo run -p gatk-haplotypecaller --example hc_full_parity_gate_dump -- \
  apply-summary parity/fixtures/reference.fa parity/fixtures/sample.bam "chr1:1-5;chr1:20-25" 5
```

## Java L2 gate (optional)

Count `apply` / match `callRegion` branch in a pinned GATK run only after you add a small Java counter or aspect hook; **L1** is Rust golden vs this TSV.

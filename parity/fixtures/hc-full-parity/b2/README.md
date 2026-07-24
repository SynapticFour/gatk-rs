# B.2 — Assembly region iterator golden files

## Artifact (`*.tsv`)

Produced by `hc_full_parity_gate_dump assembly-regions` (see `gatk-haplotypecaller/examples/hc_full_parity_gate_dump.rs`).

- Line 1: **header** — `contig\tstart\tend\tis_active\textended_start\textended_end\textension`
- Following lines: one row per [`AssemblyRegion`](../../../../gatk-haplotypecaller/src/assembly_region_iterator.rs) in **full walker traversal** order (all shards; same path as B.4).

## Cases manifest

[`cases.tsv`](./cases.tsv) columns:

- `case_id`
- `ref`, `bam` — paths relative to repo root
- `interval_cli`
- `padding` — `-` = default **100**
- `expected_tsv`

## Regenerating expected (Rust)

Use `PARITY_CARGO_TARGET_DIR` (or repo `target/`) like the B.1 README. Example:

```bash
export PARITY_CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-$PWD/target}"
exp=parity/fixtures/hc-full-parity/b2/expected/chr1_5_15.sample_bam.tsv
mkdir -p "$(dirname "$exp")"
cargo run -p gatk-haplotypecaller --example hc_full_parity_gate_dump -- \
  assembly-regions parity/fixtures/reference.fa parity/fixtures/sample.bam chr1:5-15 > "$exp"
# Optional 5th arg: shard padding `u64` (default 100), aligned with B.1/B.3.
```

## Edge fixtures (B.2.4 / B.2.5)

- **`b2_read_prefix5.sam`** — one 5 bp read at `chr1:1`; intervals past the aligned span yield inactive / empty pileup rows.
- **`b2_no_reads.sam`** — header-only BAM (RG present for Java `HcContext`); exercises zero-coverage shard traversal.
- **`b2_softclip_2s3m.sam`** — `2S3M` read; pileup uses LIBS / `AlignmentContext` semantics (not “all reads on contig”).

**B.2.2** zero-depth locus rows are additionally gated under **`b2-empty-locus/`** (same `locus-pileup` subcommand).

## Java parity (L2)

Strict **`run_hc_full_parity_l2.sh`** compares Rust dumps to frozen `java_dumps/b2/<case_id>_<PIN>.tsv` from `HcFullParityGateDump assembly-regions`. Refresh: `scripts/parity/run_hc_full_parity_java_refresh.sh`.

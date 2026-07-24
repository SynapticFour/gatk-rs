# B.1 — Read shard golden files

## Artifact (`*.tsv`)

UTF-8, no BOM. One line per **padded span** (after merge + pad), in shard order (sequence-dictionary order of contigs that appear in `-L`).

Columns (tab-separated):

1. `contig`
2. `span_start` — 1-based inclusive
3. `span_end` — 1-based inclusive

## Cases manifest

[`cases.tsv`](./cases.tsv) columns:

- `case_id`
- `ref` — path relative to repo root
- `interval_cli` — same grammar as `parse_intervals_cli_string` (`;`-separated includes)
- `padding` — decimal `u64` (GATK default = **100**)
- `expected_tsv` — path relative to repo root

## Regenerating expected (Rust)

From repo root (recommended: `export PARITY_CARGO_TARGET_DIR="$PWD/target"` so `hts-sys` bindings stay under the workspace target):

```bash
while IFS=$'\t' read -r cid ref iv pad exp; do
  [[ -z "$cid" || "$cid" == \#* ]] && continue
  mkdir -p "$(dirname "$exp")"
  cargo run -p gatk-haplotypecaller --example hc_full_parity_gate_dump -- \
    read-shards "$ref" "$iv" "$pad" > "$exp"
done < parity/fixtures/hc-full-parity/b1/cases.tsv
```

Review the diff, commit only intentional semantic changes.

## Java parity (optional L2)

GATK does not emit read-shard bounds as a stable CLI file. To compare against Java **L2**:

1. Add a tiny **Java main** (or Picard-style tool) under `scripts/parity/java/` that constructs the same `MultiIntervalLocalReadShard` / padding as `AssemblyRegionWalker` and prints this TSV format; or
2. Record shard span bounds from a one-off `DEBUG` log line and normalize to this schema.

Until L2 exists, **L1** = byte match of Rust dump to `expected/*.tsv` (CI + `run_hc_full_parity_b1_read_shards.sh`).

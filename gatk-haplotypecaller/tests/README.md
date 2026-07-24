# HaplotypeCaller test harness layers (Sprint L-2)

Three layers — keep new tests in the right bucket.

| Layer | Purpose | Where |
|-------|---------|--------|
| **Unit** | Algorithm predicates, GLs, EventMap sync — no BAM/FASTA fixtures | `tests/genotyping/`, `tests/discovery/` (pulled into the crate via `#[path]` for private access); also `src/**/#[cfg(test)]` |
| **Parity oracle** | Compare Rust vs Java dumps / TSV oracles — must not define production emit | `parity/reports/`, `parity/fixtures/`, scripts under `scripts/parity/`; gated by `--features parity_harness` |
| **P12 fixture** | NA12878 chr2:92.3M integration / sign-off | Top-level `tests/p12_*.rs`, `tests/post_shadow_*`, phase probes |

## Unit extraction (L-1)

| Path | Included from |
|------|----------------|
| `genotyping/engine_unit.rs` | `hc_genotyping_engine` (`#[path]`) |
| `discovery/event_discovery_unit.rs` | `read_event_discovery` (`#[path]`) |

Run: `cargo test -p gatk-haplotypecaller --lib`

## Parity bridges (L-4)

`CallRegionArgs::parity_aligned()` / `legacy_read_bridges()` exist only under `cfg(test)` or `--features parity_harness`.

```bash
cargo test -p gatk-haplotypecaller --features parity_harness
```

## Oracle ≠ emit (L-3)

```bash
python3 scripts/parity/oracle_emit_audit.py
```

See [`docs/CLAIM_MATRIX.md`](../../docs/CLAIM_MATRIX.md) for asserted gates.

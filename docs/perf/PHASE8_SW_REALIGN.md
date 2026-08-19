# Phase 8 — SW / realign leaf

Observable Java contract unchanged: same SoftClip / Indel SW scores, same
`createReadAlignedToRef` CIGAR+pos updates, same clipped bases for SW.

## Cuts

| Area | Change | Why |
|------|--------|-----|
| SW TLS | Single zero of `sw`/`btrack` (ensure no longer double-fills) | Cut fill work on every DP align |
| `last_index_of` | Bookend gate + interior slice eq | SoftClip exact-hit path (common realign) |
| Realign | `record_cigar_differs` iterates BAM CIGAR (no Vec) | Hot per-read |
| Clip-for-SW | Drain in place instead of `to_vec` window | Soft-clipped reads |

## Also in this PR

- wall-losers Finalize concat fix — see [`WALL_LOSERS_F1_DIG.md`](WALL_LOSERS_F1_DIG.md)

## Prove

- `cargo test -p gatk-haplotypecaller smith_waterman`
- Existing realign / unclip unit tests
- Optional: TRACE rematch w09 `prep_realign` share vs phase7
- Follow-up rolling scores / HC length benches: [`SW_HC_OPTIMIZE.md`](SW_HC_OPTIMIZE.md)

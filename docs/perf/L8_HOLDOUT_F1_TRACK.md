# L8 holdout F1 track (separate from wall)

**Do not conflate** with beat-Java wall. **Do not widen P12 bands.**

## Status (phase6 / FastestAvailable default → emit AD fix)

| Slice | Interval | Java sites | Rust sites | Rust F1 | Gate |
|-------|----------|----------:|----------:|--------:|------|
| Holdout | `20:15000000-15050000` | 55 | **51** | **~0.904** | **pass** (local) |

Prior fail: rust sites 41 / F1≈0.844. Dig + fix: [`CALLRATE_EMIT_AD_DIG.md`](CALLRATE_EMIT_AD_DIG.md).

Regen: `P12_SKIP_JAVA=1 ./scripts/parity/run_hc_full_parity_j6_dense_holdout.sh`  
Eval: `scripts/parity/p13_truth_eval.py` + `thresholds_dense_holdout.json`

Precision stays high; **recall** is the gap (shared≈40 of 55 Java sites).

## Working hypothesis (evidence-class)

Missing sites cluster as assembly/EventMap under-call vs Java on this offset window —
same class as historical Peak-cut k-best path-edge collapse, **not** PairHMM SIMD
numerics (unit SIMD↔scalar green). Treat as L8 generalize / call-rate, not wall TLS.

## Allowed next moves

1. Diff missing loci: Java−Rust set; classify SNP vs indel; EventMap / k-best / filter.
2. Holdout TRACE on one missing locus region (small −L) vs Java dump.
3. Keep L9 FORMAT+F1 battery green on primary dense before claiming holdout fixed.

## Rejected

- Widening holdout F1 thresholds or P12 bands to “pass”
- Flipping PairHMM default away from FastestAvailable to chase F1
- Blocking wall PRs on holdout alone (track in parallel)

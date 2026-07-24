# Fair HaplotypeCaller comparison (dedicated host)

**Status:** no dedicated-host fair comparison published yet.

Run on label `gatk-rs-benchmark` via [`.github/workflows/benchmark.yml`](../../.github/workflows/benchmark.yml)
or locally:

```bash
./scripts/perf/run_fair_hc_comparison.sh
```

Contract:

| Item | Value |
|------|--------|
| Configs | `rust_logless_scalar`, `rust_simd`, `java_fastest_available`, `java_logless_caching` |
| Repeats | ≥5 → **median ± sample stdev** |
| Regions | small / medium / large (nested GIAB windows) |
| Metrics | wall, user, sys, Peak-RSS, optional RAPL energy |
| Primary Java baseline | **`FASTEST_AVAILABLE`** (native AVX verified) |

Host doctrine: [`docs/ci/PERF_BENCHMARK_HOST.md`](../ci/PERF_BENCHMARK_HOST.md).

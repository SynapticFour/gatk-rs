# Performance host specifications

**Status:** not yet captured on the dedicated `gatk-rs-benchmark` host.

Run [`scripts/perf/capture_host_specs.sh`](../../scripts/perf/capture_host_specs.sh)
on that machine (or wait for [`.github/workflows/benchmark.yml`](../../.github/workflows/benchmark.yml))
to replace this stub with CPU model, RAM, kernel, governor, SMT, and AVX2/AVX-512 flags.

Doctrine: [`docs/ci/PERF_BENCHMARK_HOST.md`](../ci/PERF_BENCHMARK_HOST.md).

> Do not cite timing / Peak-RSS numbers as dedicated-host results until this
> file (or a stamped copy under `docs/perf/runs/dedicated_*/`) lists real hardware.

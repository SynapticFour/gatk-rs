# Cargo bench profile audit

## Current settings (workspace `Cargo.toml`)

```toml
[profile.release]
lto = true
codegen-units = 1
panic = "abort"
opt-level = 3

[profile.bench]
inherits = "release"
debug = false
lto = false
codegen-units = 16
```

So **`cargo bench` is still not full release-equivalent** (no fat LTO, cgu=16), but
**`debug = false` is now default** (ledger **M1** / focus matrix: ~23% PairHMM Criterion
win vs debuginfo). LTO stays off for CI/Docker memory. Overrides via
`CARGO_PROFILE_BENCH_*` / `RUSTFLAGS`.

## Reproducible matrix

| Tool | Path |
|------|------|
| Matrix runner | [`scripts/perf/run_bench_profile_matrix.sh`](../../scripts/perf/run_bench_profile_matrix.sh) |
| Focused signed run | gitignored `docs/perf/runs/bench_profile_focus_*/SUMMARY.md` (local matrix) |
| PGO workflow (opt-in) | [`scripts/perf/run_pgo_hc.sh`](../../scripts/perf/run_pgo_hc.sh) |

```bash
# Full matrix (slow; separate target dirs)
BENCH_MATRIX_QUICK=1 ./scripts/perf/run_bench_profile_matrix.sh

# Focused IDs used for signed numbers below (shared target/, 3s Criterion):
#   smith_waterman_align/soft_clip/128x96
#   assembly_graph_depth/threading_build_medium_k10
#   pairhmm simd_r100_h/32
```

| Config id | LTO | CGU | debug | `RUSTFLAGS` |
|-----------|-----|-----|-------|-------------|
| `current_bench` | false | 16 | true | — |
| `no_lto_cgu16_nodebug` | false | 16 | false | — |
| `no_lto_cgu1` | false | 1 | false | — |
| `release_equiv` | fat | 1 | false | — |
| `lto_thin` | thin | 1 | false | — |
| `release_native` | fat | 1 | false | `-C target-cpu=native` |

## Measured impact (aarch64-apple-darwin, mid Criterion time)

Ratios are **vs `current_bench`** (lower is faster). Source: focus run above.

| Kernel | current | −debug (cgu16) | cgu1 no-LTO | fat LTO | thin LTO | +native |
|--------|--------:|---------------:|------------:|--------:|---------:|--------:|
| PairHMM SIMD `simd_r100_h/32` | 373 µs | **0.77×** | 0.78× | 0.89× | **0.78×** | 0.88× |
| Smith–Waterman `128×96` SoftClip | 39 µs | 1.05×† | **0.78×** | 0.96× | 0.83× | 0.85× |
| Assembly threading `medium_k10` | 888 µs | **0.84×** | 0.91× | 0.98× | 1.03× | **0.76×** |

†SW −debug delta within noise on this host (±5%).

### Isolation takeaways

| Contrast | Finding |
|----------|---------|
| **`debug = true`** | **Largest, clear hit on PairHMM (~23% slower).** Assembly ~16% slower. |
| **`codegen-units` 16→1** (no LTO) | Small / mixed; SW liked cgu1 here; PairHMM ≈ flat. |
| **Fat LTO** | **No win** on these microkernels vs optimized no-LTO; PairHMM **regressed** vs `no_lto_cgu1`. Build cost high. |
| **Thin LTO** | ≈ no-LTO cgu1 for PairHMM; cheaper than fat. |
| **`target-cpu=native`** | Helps assembly (~24% vs current); little PairHMM gain (NEON already on by default on this host). |

Genotype AD multipass was not in the focus run (PairHMM filter previously wrong). Expect similar **debuginfo sensitivity** (CIGAR/hash loops), not LTO-dominated; re-run with:

`BENCH_MATRIX_FILTERS='genotype_dense_ad/multipass_ad/R256' BENCH_MATRIX_CONFIGS='current_bench no_lto_cgu16_nodebug release_equiv' ./scripts/perf/run_bench_profile_matrix.sh`

(requires `--features parity_harness`).

## Recommendations (still no default `Cargo.toml` edit)

### Keep three artifact classes

| Purpose | How to build | Use for |
|---------|--------------|---------|
| **Portable production** | `cargo build --release` (current: fat LTO, cgu=1, no native) | Ship, CI wall, Java compare |
| **Host-optimized production** | release + `RUSTFLAGS='-C target-cpu=native'` | Dedicated bench host only |
| **Microbenchmark (dev/CI)** | `cargo bench` | Fast iteration |
| **Publishable microbench** | matrix `no_lto_cgu16_nodebug` or `release_equiv` | Docs claiming kernel × speedups |

### What to do about `[profile.bench]`

**Do not flip default bench to fat LTO + cgu=1.** Measured: no kernel win, long link times, CI memory risk (original reason for the override).

**Preferred policy (Option B — hybrid):**

1. **Leave** `lto = false` and `codegen-units = 16` for default `cargo bench`.
2. **Optionally** set `debug = false` in a follow-up PR (only change justified by PairHMM ~1.3× on this host). Still keep LTO off.
3. **Document** that published Criterion numbers must name the config (`current_bench` vs `no_lto_cgu16_nodebug` vs `release_equiv` / `release_native`).
4. For release-parity claims, run:

```bash
CARGO_PROFILE_BENCH_LTO=fat CARGO_PROFILE_BENCH_CODEGEN_UNITS=1 CARGO_PROFILE_BENCH_DEBUG=0 \
  cargo bench -p gatk-haplotypecaller --bench pairhmm --locked -- simd_r100_h/32
```

**Do not** put `target-cpu=native` in committed `Cargo.toml`.

### Architecture builds

| Build | Flags |
|-------|-------|
| Portable aarch64 / x86_64 | default release (runtime SIMD dispatch) |
| Host-opt | `-C target-cpu=native` or pinned `apple-m1` / `znver4` |
| Avoid | baking `+avx2` into portable Linux binaries if old CPUs matter |

## PGO — promising for HC wall, not for default yet

Tight DP kernels (PairHMM/SW) already look LTO-insensitive; **branchy genotype/AD** is the better PGO target. Design only:

[`scripts/perf/run_pgo_hc.sh`](../../scripts/perf/run_pgo_hc.sh) — `prepare` → `train` (`PGO_TRAIN_CMD`) → `optimize` → compare via fair HC wall.

**Gates before any default enablement:** ≥3% wall-losers median improvement, no F1/P12 regression, opt-in CI job only.

## Status

- [x] Audit current profiles  
- [x] Matrix runner + PGO script  
- [x] Quantify PairHMM / SW / assembly on focus matrix  
- [x] Recommend hybrid bench config (**no Cargo.toml change yet**)  
- [ ] Optional follow-up PR: `profile.bench.debug = false` only  
- [ ] Genotype AD in matrix  
- [ ] PGO fair-wall Δ before default  

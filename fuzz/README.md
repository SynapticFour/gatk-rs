# Fuzz targets (Rust-native R3)

Portable smoke (stable Rust, no libFuzzer):

```bash
cargo test -p gatk-core --test fuzz_smoke_io_edges
```

Full libFuzzer targets (requires nightly + `cargo-fuzz`):

```bash
cargo install cargo-fuzz
cargo +nightly fuzz run allele_from_string   # needs --features core-targets (default)
cargo +nightly fuzz run parse_cigar_str
cargo +nightly fuzz run hc_differential      # lean: no gatk-core / htslib
```

Lean compile check (no nightly / cargo-fuzz):

```bash
cargo check --manifest-path fuzz/Cargo.toml --no-default-features --bin hc_differential
```

## HC differential fuzzer (`hc_differential`)

Finds **Java GATK4 vs gatk-rs** HaplotypeCaller divergences on synthetic BAMs
(edge CIGARs, soft-clips, mate overlap, low MAPQ, BQ jitter) — cases truth-sets
do not cover.

| Layer | Role |
|-------|------|
| `fuzz_targets/hc_differential.rs` | libFuzzer generative surface (`#[path]` → shared `scenario.rs`) |
| `gatk-rs-equiv differential-fuzz` | Full campaign: materialize BAM → both callers → shrink → fixture + optional `gh` issue |
| `fuzz/run_hc_differential.sh` | M4-safe wrapper (default 8 iterations, `RAYON_NUM_THREADS=2`) |

Campaign (requires `samtools`, Java GATK 4.4 jar or `gatk` on PATH, gatk-rs binary):

```bash
# Optional: export JAVA_GATK_JAR=/path/to/gatk-package-4.4.0.0-local.jar
./fuzz/run_hc_differential.sh --iterations 8

# On divergence: write fixture under gatk-haplotypecaller/tests/fixtures/regressions/
# and optionally open a GitHub issue (label parity-divergence):
./fuzz/run_hc_differential.sh --open-github-issue
```

Replay a minimized fixture:

```bash
./fuzz/run_hc_differential.sh \
  --replay-fixture gatk-haplotypecaller/tests/fixtures/regressions/<id>
```

LibFuzzer generative-only:

```bash
./fuzz/run_hc_differential.sh --libfuzzer
```

This crate is intentionally **not** a workspace member (cargo-fuzz convention).

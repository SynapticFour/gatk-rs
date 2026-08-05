# Rust-native HC performance showcase (parity-safe)

Living roadmap: turn `gatk-rs` into a clear Rust-over-Java showcase for HaplotypeCaller
**without** inventing non-Java downsampling, widening P12 bands, or marketing genome-wide
equivalence from smoke RSS.

Doctrine: workspace rule `rust-native-algorithm-parity` (under `.cursor/rules/`),
[`AIR_M4_GIAB_RECIPE.md`](AIR_M4_GIAB_RECIPE.md), [`PAIRHMM_SPEEDUP.md`](PAIRHMM_SPEEDUP.md),
[`CLAIM_MATRIX.md`](../CLAIM_MATRIX.md).

## Mission constraints

| Allowed | Forbidden |
|---------|-----------|
| Ownership / `Arc` / COW / `&[u8]` layouts | Non-Java force-downsample |
| TLS scratch bounds + release | Widening P12 bands (N-7) |
| Feature-gated SIMD with differential proof | Silent Log10 → SIMD default flip |
| jemalloc / arena tuning (measured) | Locus-overfit “pins” as design |
| Rayon hygiene (no nested oversubscription) | Public Peak-RSS claims from smoke only |

Observable contract: pinned GATK 4.4 HaplotypeCaller math + emit gates (L2–L5; GIAB when signed).

## Why Rust wins here

| Lever | Advantage | HC surface |
|-------|-----------|------------|
| Ownership / `Arc` / COW | Share BAM evidence; explicit mutation cost | `shared_bam`, progressive release, realign |
| Zero-copy slices | DP / threading without UTF-8 `String` tax | PairHMM; assembly spine |
| TLS + drop | Bound per-thread DP planes | PairHMM / SW scratch |
| Feature-gated SIMD | Safe NEON/AVX2 without C UB soup | `pairhmm_simd/` (opt-in until GIAB) |
| Global allocator | Lower Peak-RSS fragmentation vs glibc | `jemalloc` feature on CLI |
| Rayon | Parallel regions without hand-rolled pools | Sequential hap score under pressure |

## Phases

### A — Memory ownership (Peak-RSS)

1. Byte-native `AssemblyRead` / kmers (`Vec<u8>` / `&[u8]`); UTF-8 only at VCF emit.
2. Single `finalizeRegion` buffer shared by assemble + PairHMM.
3. COW-aware trim / realign (unique ownership before `make_mut`).
4. Sequential haplotype scoring when `GATK_RS_HC_SEQUENTIAL=1`.
5. Feature-gated jemalloc + MALLOC docs for GIAB / Air.

### B — CPU (parity-gated)

- Keep production default **Log10** PairHMM until signed GIAB/hap.py with `FASTEST_AVAILABLE`.
- Harden SIMD packing + Criterion phenotype benches (read len / hap count classes).
- Index-based / borrowed kmers; avoid unnecessary graph clones in k-best.

### C — Scalability / I/O

- Safer shard streaming within GATK iterator semantics.
- Share BAM header across workers; avoid per-region reopen where possible.
- Public claims only from dedicated benchmark host + realistic windows.

## Claim rules

- No marketing “X% less memory” from `chr1:1-32` smoke.
- GIAB `ci-subset` stays unsigned in `CLAIM_MATRIX.md` until Peak-RSS + hap.py gates.
- Never widen P12 bands. Prefer “Java does X on this evidence class” over locus pins.
- Measure: bomb / 50 kb / 100 kb / 500 kb with `GATK_RS_HC_SEQUENTIAL=1`, threads 1–2.

## Local recipe

See [`AIR_M4_GIAB_RECIPE.md`](AIR_M4_GIAB_RECIPE.md). With jemalloc builds:

```bash
cargo build -p gatk-cli --release --locked --features jemalloc
export GATK_RS_HC_SEQUENTIAL=1 RAYON_NUM_THREADS=1
# optional: MALLOC_CONF=background_thread:true,dirty_decay_ms:1000,narenas:2
```

## Measurement status (Phase A land)

| Check | Result |
|-------|--------|
| `scripts/ci/check_hc_rss_regression.sh` | OK (unit + dict Peak-RSS); HC bomb window skipped when realistic BAM/ref not staged |
| Excellence N-2 / N-5 / N-7 | PASS (env allowlist; unwrap ratchet; P12 band freeze untouched) |
| Dense bomb / 50 kb / 100 kb / 500 kb Peak-RSS | **Pending** — requires staged `parity/realworld/...` BAM + `hs37d5.simple.fa`; record numbers in [`HC_MEMORY_PROFILE.md`](HC_MEMORY_PROFILE.md) only after measurement |
| GIAB `ci-subset` | Still **unsigned** in [`CLAIM_MATRIX.md`](../CLAIM_MATRIX.md) |

Holdout windows (chr21 + chr20 offset) stay under L8/L9 signoff scripts; do not delete P12-scoped pins from Phase A ownership work.

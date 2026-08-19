# HaplotypeCaller performance ledger

Doctrine: optimize by `wall × workload fraction × achievable speedup`, not
microbench-only. No P12 widening, no semantic shortcuts, no nested Rayon.

## Methodology (step 1)

| Item | Policy |
|------|--------|
| Product wall | thr=2, unset `GATK_RS_HC_SEQUENTIAL`, `FASTEST_AVAILABLE` |
| Java baseline | GATK 4.4 `FASTEST_AVAILABLE` (native PairHMM verified) |
| Signed product lane | `GIAB_MODE=wall-losers` |
| Microbench default | `[profile.bench]`: release opts, **`debug=false`**, `lto=false`, cgu=16 |
| Publishable microbench | Name config (`current_bench` / `no_lto_cgu16_nodebug` / `release_equiv`) |
| Production profile | `scripts/perf/run_hc_profile.sh` → `hc_profile.json` + TRACE |

## Top-3 wall contributors (evidence class)

Sources: production profiles (gitignored local `hc_profile` under `docs/perf/runs/`);
local TRACE mega ([BEAT_JAVA_WALL_NEXT.md](BEAT_JAVA_WALL_NEXT.md));
phase7 Instruments ([PHASE7_WALL_PROFILE.md](PHASE7_WALL_PROFILE.md)).

### Product window w09 `21:9.5–9.7Mb` (thr=2, FASTEST_AVAILABLE) — tip profile

Run wall **57.6 s**. Stage walls (nested guards; can overlap TRACE):

| Rank | Stage | Wall s | Notes |
|-----:|-------|-------:|-------|
| **1** | **PairHMM** | **22.3** | Pack occupancy **4.4%**; prefix reuse **94.1%** |
| **2** | **Smith–Waterman** | **21.9** | Realign |
| **3** | **Haplotype gen / event discovery** | **~12.8 each** | Assemble prep |
| — | Genotype assignment | 3.2 | Nested AD **0.12 s** (memo working) |

### Remaining CI losers (200 kb heads; thr=2) — tip profile

| Window | Run wall | #1 PairHMM | #2 SW | #3 Hap/event | Pack occ | Prefix |
|--------|--------:|----------:|------:|-------------:|---------:|-------:|
| w11 `21:11.0–11.2Mb` | 146.6 s | 118.3 | 92.2 | ~20.6 / 19.8 | 3.0% | 94.6% |
| w26 `20:26.0–26.2Mb` | 40.9 s | 26.7 | 17.0 | ~9.5 / 7.8 | 6.3% | 91.1% |
| w29 `20:29.4–29.6Mb` | 103.1 s | 68.0 | 53.4 | ~19.6 / 20.7 | 3.0% | 95.1% |

Claim lane remains full 1 Mb wall-losers. Artifacts gitignored under `docs/perf/runs/`.

### Mega dense TRACE (pre–G1 tip class → post-G1)

| Rank | Contributor | Pre-G1 | Post-G1 tip |
|-----:|-------------|--------|-------------|
| 1 | Genotype assign / AD | ~130–172 s | **~1.2 s** (no longer dominant) |
| 2 | PairHMM | ~9–20 s | ~13 s |
| 3 | Realign SW | ~4–18 s | ~4.7 s |

**Priority rule:** use the measured profile for the target claim window. Post-rematch
wall-losers median **~1.15×** → PairHMM pack occupancy + SW first on remaining
slower shards (`01_chr21_w11`, `00_chr20_w26/w29`).


## Signed Java rematch (CI) — **current tip**

`GIAB_MODE=wall-losers` @ `perf/hc-wall-ledger-g1-p1-p2` tip
([run 32244936304](https://github.com/SynapticFour/gatk-rs/actions/runs/32244936304)):

| Metric | Value |
|--------|------:|
| Median wall Rust/Java | **~1.15×** |
| Σ wall | **1.27×** (32m54 / 25m53) |
| Worst | `01_chr21_w11` **1.65×** |
| Best | `01_chr21_w29` **0.82×** (Rust faster) |
| Peak RSS | ~0.14–0.66× Java (wins) |
| Equivalence | PASS `max_\|ΔF1\|=0.0004` |

### vs prior signed baseline (`6a478ca`, run 32003572266)

| Metric | Prior | Tip | Δ |
|--------|------:|----:|---|
| Median | 1.61× | **1.15×** | −0.46 |
| Σ | 1.65× | **1.27×** | −0.38 |
| Worst | w09 2.01× | w11 1.65× | — |
| Rust Σ wall | 2654 s | 1974 s | **−680 s** |
| Equivalence | PASS 0.0004 | PASS 0.0004 | flat |

Per-shard Rust/Java: w09 **2.01→1.06×**; three shards within ±10% of Java; one
Rust win. Note: `01_chr21_w11` ratio **1.41→1.65×** (Java abs dropped more than
Rust; Rust still −57 s wall).

Prior baseline kept for history: median **~1.61×**, Σ **1.65×**, worst w09 **2.01×**.

## Experiment ledger

| ID | Date | Change | Workload | Baseline | New | Speedup | Peak RSS | Alloc | Parity | Keep? |
|----|------|--------|----------|----------|-----|---------|----------|-------|--------|-------|
| M1 | 2026-08-19 | `[profile.bench] debug=false` | Criterion PairHMM `simd_r100_h/32` (focus matrix) | ~373 µs (debug=true) | ~0.77× prior | **~1.30×** micro | n/a | n/a | n/a (build) | **KEEP** |
| G1 | 2026-08-19 | Likelihood TLS borrow cache + empty `read_id`/`read_index`; indel AD memo; SiteScore uses `with_region_likelihood_rows` | `genotype_likelihood_reshape` Criterion | uncached R512×H32×A16 **659 µs** | cached borrow **44 ns** | **~15 000×** reshape leaf | n/a | no String/format; 1 reshape/locus | marg + phase_c + ad_result_memo PASS | **KEEP** |
| P1 | 2026-08-19 | Drop dead f64 `prior` TLS (pack/AVX2/NEON/Logless) | PairHMM SIMD vs scalar | same numerics | same | RSS hygiene (~25% smaller f64 DP Peak) | expect ↓ | fewer Peak bytes | pairhmm_simd_vs_scalar + logless PASS | **KEEP** |
| P2 | 2026-08-19 | TLS-reuse PairHMM transitions (`logless_fill_transitions`) pack/NEON/AVX2 | `pairhmm_logless_simd/simd_r100_h/32` | ~1 alloc/read | 0 steady-state | small wall (amortized); **alloc↓** | n/a | 1 fewer Vec/read | pairhmm_simd_vs_scalar PASS | **KEEP** |
| W09 | 2026-08-19 | Production profile baseline (post G1/P1 tip) | w09 `21:9.5–9.7Mb` thr=2 | n/a | run wall **57.6 s**; PairHMM 22.3 / SW 21.9 / geno 3.2 (AD 0.12) | n/a | see TRACE | — | observe-only | baseline |
| MEGA | 2026-08-19 | Local rematch mega `21:9825–9828k` (same tip) | was assign **~130–172 s** | run **18.0 s**; geno assign Σ **1.19 s** (AD 0.075); PairHMM 13.1 / SW 4.7 | **~10×** region wall vs prior TRACE class; genotype no longer dominant | — | — | observe-only | **KEEP signal** |
| WL1 | 2026-08-19 | Tip rematch vs Java (G1/P1/P2 + RT/SW TLS) | `wall-losers` CI 32244936304 | median 1.61× / Σ 1.65× | **median 1.15× / Σ 1.27×** | **~1.40×** product median | still wins | — | PASS ΔF1=0.0004 | **KEEP** |
| P3a | 2026-08-19 | `PREFIX_REUSE_OVER_SIMD_FRAC` 0.35→0.50 | w11 200 kb profile | occ 3.0% / prefix 94.6% | **unchanged** | none | — | — | pairhmm_simd PASS | **REVERT** |
| P3b | 2026-08-19 | `PREFIX_REUSE_MIN_HAPS` NEON 3→6 / AVX2 5→8 | w11 200 kb profile | same pack unit counts | **identical occupancy** | none | — | — | pairhmm_simd PASS | **REVERT** (named consts @ 3/5) |
| W11+ | 2026-08-19 | Loser profiles (w11/w26/w29 heads) | local thr=2 | n/a | PairHMM≻SW≻hap/event on all three | n/a | — | — | observe-only | baseline |


### Production profile artifact

`docs/perf/runs/hc_profile_ledger_w09/` — `hc_profile.md` + JSON + TRACE (gitignored).


### G1 notes

- **Why faster:** multi-allelic sites previously rebuilt dense `R×H` rows (and `format!("read_{i}")`) once per allele. TLS cache rebuilds once; SiteScore borrows without cloning. Observable GLs unchanged.
- **Wall hypothesis:** cuts `marginalize_wall` on multi-allelic dense sites; AD pileup remains the larger mega-window leaf ([GENOTYPE_ASSIGN_COMPLEXITY.md](GENOTYPE_ASSIGN_COMPLEXITY.md)).
- **Indel AD memo:** identical `(reads,loc,ref,alt)` rescans return cached counts.
- Raw A/B: `docs/perf/runs/ledger_g1_reshape_ab.txt`

### Prove commands

```bash
cargo test -p gatk-haplotypecaller --lib ad_result_memo --locked
cargo test -p gatk-haplotypecaller --lib pairhmm_logless --locked
cargo test -p gatk-haplotypecaller --test pairhmm_simd_vs_scalar_test --locked
cargo test -p gatk-haplotypecaller --test hc_genotyping_marginalization_fixture_test --locked
cargo test -p gatk-haplotypecaller --test phase_c_genotyping_parity_test --locked
cargo bench -p gatk-haplotypecaller --bench genotype_assign --features parity_harness --locked -- genotype_likelihood_reshape
```

### Production profile (step 2)

```bash
./scripts/perf/run_hc_profile.sh --out-dir docs/perf/runs/hc_profile_ledger_w09 -- \
  -R parity/realworld/assets/hs37d5.simple.fa \
  -I parity/realworld/na12878_ci_loser_windows/01_chr21_w09.bam \
  -O /tmp/w09_ledger.vcf \
  -L 21:9500000-9700000
```

```bash
# Remaining slower CI shards (200 kb heads for local dig)
./scripts/perf/run_hc_profile.sh --out-dir docs/perf/runs/hc_profile_ledger_w11 -- \
  -R parity/realworld/assets/hs37d5.simple.fa \
  -I parity/realworld/na12878_ci_loser_windows/01_chr21_w11.bam \
  -O /tmp/w11_ledger.vcf -L 21:11000000-11200000
./scripts/perf/run_hc_profile.sh --out-dir docs/perf/runs/hc_profile_ledger_w26 -- \
  -R parity/realworld/assets/hs37d5.simple.fa \
  -I parity/realworld/na12878_ci_loser_windows/00_chr20_w26.bam \
  -O /tmp/w26_ledger.vcf -L 20:26000000-26200000
./scripts/perf/run_hc_profile.sh --out-dir docs/perf/runs/hc_profile_ledger_w29 -- \
  -R parity/realworld/assets/hs37d5.simple.fa \
  -I parity/realworld/na12878_ci_loser_windows/00_chr20_w29.bam \
  -O /tmp/w29_ledger.vcf -L 20:29000000-29200000
```

## Next candidates (ordered)

1. **PairHMM read-axis packs / wavefront** — hap-axis prefix path saturated (~3–6% packs); prefix-vs-pack knobs A/B **reverted** ([PHASE9_PAIRHMM_PACK.md](PHASE9_PAIRHMM_PACK.md))
2. **SW striped SIMD** — gate met (SW #2 on all loser heads); only with oracle + Criterion ([SW_HC_OPTIMIZE.md](SW_HC_OPTIMIZE.md))
3. Collapse `try_genotype` to single AD authority (mega windows; parity-sensitive)
4. QNAME dedupe without `Vec<u8>` ownership (Arc/intern)
5. RT arena node bases (preserve BTree neighbor order)
6. Host-native release builds for dedicated bench only; rematch wall-losers after each KEEP

## Rejected / do-not

- P12 band widening; nested Rayon; unsafe for bounds only; toy-only microbench claims

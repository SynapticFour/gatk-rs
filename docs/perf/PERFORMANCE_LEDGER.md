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

Sources: **this tip** production profile (gitignored local `hc_profile` under
`docs/perf/runs/`); local TRACE mega ([BEAT_JAVA_WALL_NEXT.md](BEAT_JAVA_WALL_NEXT.md));
phase7 Instruments ([PHASE7_WALL_PROFILE.md](PHASE7_WALL_PROFILE.md)).

### Product window w09 `21:9.5–9.7Mb` (thr=2, FASTEST_AVAILABLE) — **this tip**

Run wall **57.6 s**. Stage walls (nested guards; can overlap TRACE):

| Rank | Stage | Wall s | Notes |
|-----:|-------|-------:|-------|
| **1** | **PairHMM** | **22.3** | Pack occupancy **4.4%**; prefix reuse **94.1%** |
| **2** | **Smith–Waterman** | **21.9** | Realign |
| **3** | **Haplotype gen / event discovery** | **~12.8 each** | Assemble prep |
| — | Genotype assignment | 3.2 | Nested AD **0.12 s** (memo working) |

### Mega dense TRACE (pre–G1 tip class)

| Rank | Contributor | Wall |
|-----:|-------------|------|
| 1 | Genotype assign / AD | ~130–172 s |
| 2 | PairHMM | ~9–20 s |
| 3 | Realign SW | ~4–18 s |

**Priority rule:** use the measured profile for the target claim window. Wall-losers /
mega → genotype AD authority next. w09-class → PairHMM (transitions TLS, pack
occupancy) and SW first.


## Signed Java rematch baseline (CI, pre–local tip)

`GIAB_MODE=wall-losers` @ `6a478ca` (run 32003572266):

| Metric | Value |
|--------|------:|
| Median wall Rust/Java | **~1.61×** |
| Σ wall | **1.65×** |
| Worst | w09 **2.01×** |
| Peak RSS | ~0.16–0.62× Java (wins) |
| Equivalence | PASS `max_\|ΔF1\|=0.0004` |

**Note:** Local tip includes uncommitted AD memo / RT packing / SW rolling / etc.
CI 1.61× is **not** this tip until rematched.

## Experiment ledger

| ID | Date | Change | Workload | Baseline | New | Speedup | Peak RSS | Alloc | Parity | Keep? |
|----|------|--------|----------|----------|-----|---------|----------|-------|--------|-------|
| M1 | 2026-08-19 | `[profile.bench] debug=false` | Criterion PairHMM `simd_r100_h/32` (focus matrix) | ~373 µs (debug=true) | ~0.77× prior | **~1.30×** micro | n/a | n/a | n/a (build) | **KEEP** |
| G1 | 2026-08-19 | Likelihood TLS borrow cache + empty `read_id`/`read_index`; indel AD memo; SiteScore uses `with_region_likelihood_rows` | `genotype_likelihood_reshape` Criterion | uncached R512×H32×A16 **659 µs** | cached borrow **44 ns** | **~15 000×** reshape leaf | n/a | no String/format; 1 reshape/locus | marg + phase_c + ad_result_memo PASS | **KEEP** |
| P1 | 2026-08-19 | Drop dead f64 `prior` TLS (pack/AVX2/NEON/Logless) | PairHMM SIMD vs scalar | same numerics | same | RSS hygiene (~25% smaller f64 DP Peak) | expect ↓ | fewer Peak bytes | pairhmm_simd_vs_scalar + logless PASS | **KEEP** |
| P2 | 2026-08-19 | TLS-reuse PairHMM transitions (`logless_fill_transitions`) pack/NEON/AVX2 | `pairhmm_logless_simd/simd_r100_h/32` | ~1 alloc/read | 0 steady-state | small wall (amortized); **alloc↓** | n/a | 1 fewer Vec/read | pairhmm_simd_vs_scalar PASS | **KEEP** |
| W09 | 2026-08-19 | Production profile baseline (post G1/P1 tip) | w09 `21:9.5–9.7Mb` thr=2 | n/a | run wall **57.6 s**; PairHMM 22.3 / SW 21.9 / geno 3.2 (AD 0.12) | n/a | see TRACE | — | observe-only | baseline |
| MEGA | 2026-08-19 | Local rematch mega `21:9825–9828k` (same tip) | was assign **~130–172 s** | run **18.0 s**; geno assign Σ **1.19 s** (AD 0.075); PairHMM 13.1 / SW 4.7 | **~10×** region wall vs prior TRACE class; genotype no longer dominant | — | — | observe-only | **KEEP signal** |


### Production profile artifact

`docs/perf/runs/hc_profile_ledger_w09/` — `hc_profile.md` + JSON + TRACE.


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

Record stage walls into a new ledger row when available.

## Next candidates (ordered)

1. **PairHMM pack occupancy** on w09 (only **4.4%** SIMD packs; **94%** prefix reuse) — retune prefix-vs-pack threshold with TRACE proof
2. **SW realign** wall (~22 s on w09) — already rolling; next is striped/SIMD only with oracle
3. Collapse `try_genotype` to single AD authority (mega windows; parity-sensitive)
4. QNAME dedupe without `Vec<u8>` ownership (Arc/intern)
5. RT arena node bases (preserve BTree neighbor order)
6. Host-native release builds for dedicated bench only; wall-losers rematch vs Java

## Rejected / do-not

- P12 band widening; nested Rayon; unsafe for bounds only; toy-only microbench claims

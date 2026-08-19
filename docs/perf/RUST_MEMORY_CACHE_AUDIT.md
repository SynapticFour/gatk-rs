# Rust memory / cache audit — HaplotypeCaller

**Scope:** allocation, locality, Arc/RC, TLS, bounds checks, pointer-chasing — **not**
algorithmic redesign. No cosmetic Rustification.

**Method:** hot-path inspection of PairHMM / likelihood, read-threading / assembly,
genotype / AD. Each item has a class + measurable hypothesis.

**Already landed (do not re-propose):**
- PairHMM TLS DP planes; read-qual scratch in `likelihood_engine`
- RT `KmerKey` packing + `RollingKmer` (docs/perf/RT_KMER_PACKED_KEYS.md)
- SW rolling scores (docs/perf/SW_HC_OPTIMIZE.md)
- AD decode cache + AD result memo (docs/perf/GENOTYPE_ASSIGN_COMPLEXITY.md)

---

## Classification legend

| Tag | Meaning |
|-----|---------|
| **A** | Allocation reduction |
| **L** | Locality improvement |
| **B** | Branch reduction |
| **C** | Compiler optimization (unchecked / monomorphize / inline) |
| **K** | Cache / working-set improvement |
| **S** | Synchronization / RC / TLS-borrow reduction |

Priority: **P0** = likely wall-visible; **P1** = measurable microbench; **P2** = RSS/hygiene.

---

## 1. PairHMM / likelihood engine

### P0 — Dead TLS `prior` planes on Logless / AVX2 / NEON / pack-f64

| Evidence | Notes |
|----------|-------|
| `pairhmm_simd/avx2.rs` scratch `prior.resize` | Hot kernel uses stack `prior_arr`, never `scratch.prior` |
| `pairhmm_simd/neon.rs` same | |
| `pairhmm_simd/pack.rs` `F64Scratch.prior` | Comment: rolling path unused; f32 path still needs prior |
| `pairhmm_logless.rs` LoglessScratch | Emission via inline `logless_match_mismatch_prior` |

| Class | **K** / **A** (Peak sticky TLS) |
|-------|----------------------------------|
| **Hypothesis** | Drop unused f64 `prior` → ~20–25% smaller PairHMM TLS high-water on SIMD/Logless. Measure: Peak RSS / `hc_profile` alloc bytes on dense window; Criterion unchanged numerically. |

### P0 — Fresh `transitions` Vec every read (SIMD / pack / wavefront)

| Evidence | |
|----------|--|
| `pack.rs` `transitions = vec![[0.0;6]; rn+1]` | each `score_haps_logless_packed_*` |
| `avx2.rs` / `neon.rs` | `logless_build_transitions` → owned Vec per call |
| `wavefront/prep.rs` `ReadPrep::build` | 4 Vecs (f64/f32 transitions + match/mm) per read |

| Class | **A** / **K** |
|-------|----------------|
| **Hypothesis** | TLS-reuse transitions/`ReadPrep` → ~1–4 fewer heap allocs/read. Measure: `dhat` alloc count on `pairhmm` Criterion `simd_r100_h/32`; expect 5–15% fewer allocs when hap_count ≥ 32. |

### P1 — Per-haplotype TLS `RefCell` borrow (scalar backends)

| Evidence | |
|----------|--|
| `likelihood_engine.rs` | Outer `PAIRHMM_READ_SCRATCH` borrow OK |
| `pairhmm_logless` / `pairhmm_log10` | `.with(borrow_mut)` **inside** hap `.map` |

SIMD/pack/wavefront already borrow once around the hap loop.

| Class | **S** / **B** |
|-------|----------------|
| **Hypothesis** | Pass `&mut scratch` for all haps → remove N RefCell checks/read. Expect ~2–8% on LoglessScalar Criterion at hap_count ≥ 32. |

### P1 — Checked indexing on rolling / scalar leftover paths

| Evidence | |
|----------|--|
| `pack.rs` `score_one_f64_rolling` | `m_curr[j]` etc. checked |
| Contrast | `fill_prefix_row_*`, wavefront use `get_unchecked` after length gates |

| Class | **C** |
|-------|------|
| **Hypothesis** | Unchecked after `ensure` → ~5–15% on leftover/rolling microbench only (not main AVX pack). |

### P1 — Score `Vec<f64>` + explode to AoS triples

| Evidence | |
|----------|--|
| All score APIs return `Vec<f64>` | new alloc/read |
| `engine.rs` PairHMM region | scores → `RegionReadLikelihood` triples → later dense reshape |

| Class | **A** / **L** |
|-------|----------------|
| **Hypothesis** | Write into caller/TLS out buffer; keep dense `R×H` for genotyping → cut R allocs + better scan locality. Measure: PairHMM stage wall + genotype `marginalize_wall` on profiled mega window. |

### P2 — Length-group `Vec<&[u8]>` / `lengths` collect (AVX2/NEON)

| Evidence | `avx2.rs` / `neon.rs` by_len HashMap (TLS-reused) then `subset: Vec<&[u8]>` per group |

| Class | **A** |
|-------|------|
| **Hypothesis** | Score by indices into `haplotypes` → fewer tiny Vecs when many length groups. Small wall unless length diversity is high. |

### Already good

- Read-qual TLS scratch (`likelihood_engine`); no `to_vec` on HC score path  
- Single TLS DP borrow for packed/AVX2/NEON/wavefront across haps  
- Contiguous SoA DP planes (working set is the residual cache issue, not pointers)

---

## 2. Read-threading / assembly graph

### Already fixed

- Sliding-window map keys: `KmerKey` packed + `RollingKmer` (zero Arc on ACGT k≤64)  
- Unique maps: `HashMap`/`HashSet` (not BTree)  
- Pending sequences consumed via `take` (not clone)

### P0 — Node payloads remain `Arc<[u8]>` + `materialize_arc` per vertex

| Evidence | |
|----------|--|
| `kmer_key.rs` `materialize_arc` | packed → decode → `Arc::from` |
| `read_threading_graph.rs` `create_vertex` | nodes `Vec<Arc<[u8]>>` |
| `assembly.rs` `KmerNode.kmer` | Arc payload |

| Class | **A** / **S** / **L** |
|-------|------------------------|
| **Hypothesis** | Intern suffix byte or store packed key + arena bytes → Arc count ≈ segments only, not vertices. Measure: dhat Arc::from on `threading_build_high_*`; expect Arc count ≈ `#usable_segments` not `#windows+#vertices`. |

### P0 — Usable-segment still `Arc::from` slice

| Evidence | `read_threading_graph.rs` `sequences_from_read` |

| Class | **A** / **S** |
|-------|----------------|
| **Hypothesis** | Index into read buffer / region arena → 0 Arc/segment on ACGT. Measure: dhat at high depth. |

### P1 — Adjacency `HashMap<usize, BTreeSet<usize>>` + per-step `collect::<Vec<_>>`

| Evidence | |
|----------|--|
| `read_threading_graph` / `assembly` | BTreeSet **required** for first-suffix-match order |
| `extend_chain_by_one` | neighbors collected to Vec repeatedly |

| Class | **A** / **L** / **K** |
|-------|------------------------|
| **Hypothesis** | Keep ordered neighbors as sorted `Vec<usize>` (append+sort once) or CSR with sorted runs → fewer BTree node hops; iterate without collect. Measure: `threading_build` + LLC misses; **must** preserve iteration order (parity). |

### P1 — Unreserved HashMap growth during build

| Evidence | Builder maps start empty (`unique_kmers`, `edges`, `outgoing`) |

| Class | **A** / **K** |
|-------|----------------|
| **Hypothesis** | Reserve from `#pending windows` estimate → fewer rehashes. Measure: dhat realloc count on `threading_build_high_k10`. |

### P1 — Bench contamination: `edges_sorted()` inside threading_build

| Evidence | `benches/assembly_graph.rs` `threading_build_only` |

| Class | **A** (measurement) |
|-------|---------------------|
| **Hypothesis** | Removing sort from build bench → lower reported build time (fairer). Not a product change. |

### P2 — Junction-tree path still pre-packing (`Vec<u8>` / `BTreeMap` / snapshot clone)

| Evidence | `junction_tree_graph.rs` |

| Class | **A** / **L** |
|-------|----------------|
| **Hypothesis** | Port `KmerKey` + HashMap unique maps → same class of win as RT. Only if JT path is on production wall. |

### Do **not** replace

- `BTreeSet` adjacency without an equivalent deterministic tie-break (observable topology).

---

## 3. Genotype / AD / allele mapping

### Already fixed

- AdDecodeCache (CIGAR/seq once per Arc)  
- AdResultMemo (identical AD/softclip rescans)  
- Pad/slice short-circuits in `try_genotype`

### P0 — `format!("read_{i}")` + dense reshape per allele

| Evidence | |
|----------|--|
| `hc_genotyping_engine/mod.rs` `region_likelihoods_to_rows` | String + `Vec<f64>` per read |
| `genotype_site_score.rs` | calls reshape **every allele** |
| Also | region_summary, sparse paths, narrow helpers |

| Class | **A** / **L** |
|-------|----------------|
| **Hypothesis** | Store `ReadIndex` (no String); cache one dense `R×H` (or filtered view) per locus → cut allocs ≈ `R × (1 + A×passes)`. Measure: `dhat` + `marginalize_wall_s` on mega / `genotype_marginalize_pl` bench. Aligns with GENOTYPE_ASSIGN deferred item. |

### P0 — QNAME `HashSet<Vec<u8>>` on AD miss / LL dedupe

| Evidence | |
|----------|--|
| `read_allele_depths_*_dedupe` | `qname.to_owned()` into HashSet |
| Softclip | same (memoized) |
| `dedupe_likelihood_subset_by_qname` | **not** memoized — often per allele |

| Class | **A** / **K** |
|-------|----------------|
| **Hypothesis** | Key by `Arc` ptr (+ mate flag) or interned qname id → fewer String/Vec keys. Measure: alloc rate on memo-miss AD + dense multi-allelic LL dedupe. |

### P1 — Likelihood subset `.into_owned()` / `.cloned()` / `subset.clone()`

| Evidence | `genotype_site_pipeline.rs`, early_template gap paths |

| Class | **A** |
|-------|------|
| **Hypothesis** | Index bitset / in-place filter → O(R) fewer shallow clones per allele. Modest wall vs AD; visible in alloc counters. |

### P1 — Contiguous allele / event ownership

| Evidence | `try_genotype` takes owned `VariationEvent` (clones at call sites) |

| Class | **A** |
|-------|------|
| **Hypothesis** | `&VariationEvent` → remove S·A event clones. Tiny wall; cleaner API. |

### P2 — TLS memo double-entry on miss

| Evidence | `ad_result_memo` drop → compute → re-borrow to insert |

| Class | **S** |
|-------|------|
| **Hypothesis** | Single borrow with entry API → noise-level only vs O(R·C) AD. |

### Not a genotype-assign Arc problem

- Production AD takes `&[SharedBamRecord]`; RC bumps are upstream of assign.

---

## 4. Cross-cutting themes

| Theme | Where | Prefer |
|-------|-------|--------|
| Contiguous / index graphs | RT adjacency, hap×read scores | Arena + `u32` indices; CSR with **sorted** neighbor runs |
| Avoid pointer payload in hot maps | K-mer keys (done), node bases (open) | Packed keys / suffix bytes; Arc only for long-lived share |
| TLS scratch reuse | PairHMM transitions, wavefront prep | Same pattern as DP planes |
| One reshape, many consumers | Genotype rows | Dense matrix owned once per region/locus |
| Measure before “arena everything” | — | dhat + Criterion + `hc_profile` nested walls |

---

## 5. Suggested experiment order (no semantics change)

| # | Change | Prove |
|---|--------|-------|
| 1 | Drop dead f64 `prior` TLS planes | Peak RSS / PairHMM TLS; SIMD unit gate |
| 2 | TLS transitions / `ReadPrep` | `dhat` + `pairhmm` Criterion |
| 3 | Genotype: kill `format!(read_*)` + cache rows/locus | `marginalize_wall` + mega TRACE |
| 4 | Scalar PairHMM single scratch borrow | LoglessScalar Criterion |
| 5 | RT reserve HashMaps + avoid neighbor `collect` | `threading_build` + dhat |
| 6 | QNAME key intern / Arc-ptr for LL dedupe | dense multi-allelic alloc profile |
| 7 | Optional: index-based node bases (arena) | RT high-depth; **parity** indel/SNP topology |

---

## 6. Explicit non-goals (this audit)

- Farrar/striped SW or PairHMM algorithm swaps  
- Replacing BTree adjacency without order-preserving substitute  
- Cosmetic iterator/`into_iter` churn without alloc/cache evidence  
- PGO / LTO profile changes (see `BENCH_PROFILE_MATRIX.md`)

---

## Status

Audit only — **no code changes** in this pass. Implement from the experiment table with gates above.

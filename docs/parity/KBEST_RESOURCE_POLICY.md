# SeqGraph k-best resource policy (6R.140, production default 7R.22)

Unset `GATK_RS_EXPERIMENTAL_KBEST_POLICY` selects `unbounded_diagnostic`: no frontier-count refusal. That is the production default. Setting the variable to `unbounded_diagnostic` selects the same policy. `legacy_1024` is the explicit historical heap-length refusal, not the default. `byte_budget` is an explicit experimental byte ceiling. No new frontier cap was added. The frontier-default record is `docs/parity/7R.22_JAVA_EQUIVALENT_KBEST_DEFAULT.md`. That note still says the RSS abort was unchanged; the production RSS contract below is the current one.

**`docs/CLAIM_MATRIX.md` is unchanged.**

This is not a product claim, a public CLI, or an approved byte budget.

## What is separated

| axis | meaning | production default |
|------|---------|--------------------|
| **K** | requested completed paths / per-vertex visit budget | unchanged (`num_best_haplotypes_per_graph`, typically 128) |
| **Resource policy** | how much **frontier** may be live | `unbounded_diagnostic` (no frontier-count refusal) |

Do not equate K with the 1024 heap-object proxy. That proxy is the `legacy_1024` refusal, and it applies only when that policy is selected explicitly.

## Policies

| name | frontier refusal | how it is selected |
|------|------------------|--------------------|
| `unbounded_diagnostic` | none. The production search still skips expansion at path-edge cap 4096 and expansion cap 12_000. Those caps are separate from this frontier policy. | **production default** when the policy env is unset. Explicit `unbounded_diagnostic` is the same policy. |
| `legacy_1024` | `BinaryHeap.len() >= 1024` (same timing as pre-6R.140 SeqGraph k-best) | explicit `GATK_RS_EXPERIMENTAL_KBEST_POLICY=legacy_1024` only |
| `byte_budget` | experimental live-frontier **bytes** ≥ the caller ceiling | explicit `byte_budget` plus `GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET`. **No ceiling is a product budget.** |

## Env (internal / experimental)

Unset, empty, `0`, or `off` selects `unbounded_diagnostic`. Invalid values **fail closed** (configuration error). Not a public CLI.

| variable | role |
|----------|------|
| `GATK_RS_EXPERIMENTAL_KBEST_POLICY` | `legacy_1024` \| `unbounded_diagnostic` \| `byte_budget` |
| `GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET` | positive integer **bytes**; required iff `byte_budget`; forbidden otherwise |
| `GATK_RS_KBEST_DIAGNOSTICS` | `1` → stderr line `KBEST_RESOURCE_DIAG ...` |

## Accounting (experimental, copied `PathState`)

Counted in `peak_frontier_bytes` (and the byte-budget check):

- `heap.len() × size_of::<HeapItem>()` (includes `PathState` / `Vec` header)
- `Σ edges.capacity() × 16` (owned edge payload)

**Not** counted: malloc rounding, `BinaryHeap` capacity slack (`peak_heap_vec_capacity_bytes` is diagnostic only), completed result paths, SeqGraph, PairHMM, BAM, fragmentation.

Process RSS (`rss_before` / `rss_peak` / `rss_after`) is sampled around the run. **It is not k-best-attributable** merely because it was recorded (6R.138).

## RSS abort

`GATK_RS_HC_RSS_ABORT_MIB` is a separate guard from the frontier policy. Unset leaves the guard off.

When it is set and process RSS is at or above the limit, production k-best returns `GatkError::memory` with `hc_rss_limit_exceeded_message` and does not return the partial path vector as success. That is `find_best_haplotypes_seq_graph`, `find_best_haplotypes_seq_graph_with_policy`, and the read-threading k-best search. The assembler and RT-first checks that use the same message return that memory error as well.

`find_best_haplotypes_seq_graph_forensic` and `find_best_haplotypes_seq_graph_compact` still stop the search, set `rss_aborted`, and return `Ok` with the paths collected so far.

## Logical vs heap vs bytes

A future compact representation must **not** redefine the resource contract as `BinaryHeap.len()`.

| representation | logical frontier | heap entries |
|----------------|------------------|--------------|
| copied / compact eager | live unpopped partials | same |
| compact lazy heads | live unpopped **+ deferred siblings** | heap heads only |

6R.136: logical peak 1812 with lazy heap peak 1529.

## API

- Production: `find_best_haplotypes_seq_graph` → `SeqKbestResourcePolicy::from_runtime()` (unset → `unbounded_diagnostic`).
- Explicit policy, including `legacy_1024` and `byte_budget`: `find_best_haplotypes_seq_graph_with_policy`.
- Forensic RSS abort still returns `Ok` with `rss_aborted`: `find_best_haplotypes_seq_graph_forensic` and `find_best_haplotypes_seq_graph_compact`.
- `DIAGNOSTIC_GENEROUS_BYTE_BUDGET` (8 MiB) is **only** for 6R.140 canonical experiments. It is not a selected contract.

## Next

Curves across budgets/workloads are **6R.141**, only if this instrumentation stays trustworthy. Do not treat 8 MiB, 839 KiB, or 1812 states as policy.

//! GATK `ChainPruner`, `LowWeightChainPruner`, and `AdaptiveChainPruner` on [`crate::assembly::AssemblyGraph`].

use crate::assembly::{AssemblyGraph, AssemblyGraphPruningParams};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet, VecDeque};

#[derive(Debug, Clone, Copy)]
struct ChainEdge {
    from: usize,
    to: usize,
    support: u32,
}

#[derive(Debug, Clone)]
struct LinearChain {
    edges: Vec<ChainEdge>,
    first: usize,
    last: usize,
}

fn find_all_chains(graph: &AssemblyGraph) -> Vec<LinearChain> {
    let sources = graph.source_nodes();
    let mut chain_starts: VecDeque<usize> = sources.into();
    let mut chains = Vec::new();
    let mut seen: HashSet<usize> = chain_starts.iter().copied().collect();

    while let Some(start) = chain_starts.pop_front() {
        for to in graph.outgoing_nodes(start) {
            if let Some(edge) = graph.edge_pruning_support(start, to) {
                let chain = extend_chain(
                    graph,
                    ChainEdge {
                        from: start,
                        to,
                        support: edge,
                    },
                );
                let end = chain.last;
                chains.push(chain);
                if seen.insert(end) {
                    chain_starts.push_back(end);
                }
            }
        }
    }
    chains
}

fn extend_chain(graph: &AssemblyGraph, start: ChainEdge) -> LinearChain {
    let mut edges = vec![start];
    let first = start.from;
    let mut last = start.to;
    loop {
        let outs = graph.outgoing_nodes(last);
        if outs.len() != 1 {
            break;
        }
        let next_to = outs[0];
        if graph.incoming_count(last) > 1 || next_to == first {
            break;
        }
        let Some(support) = graph.edge_pruning_support(last, next_to) else {
            break;
        };
        edges.push(ChainEdge {
            from: last,
            to: next_to,
            support,
        });
        last = next_to;
    }
    LinearChain { edges, first, last }
}

fn low_weight_chains_to_remove(
    chains: &[LinearChain],
    graph: &AssemblyGraph,
    prune_factor: u32,
) -> Vec<LinearChain> {
    chains
        .iter()
        .filter(|c| {
            c.edges
                .iter()
                .all(|e| e.support < prune_factor && !graph.edge_is_ref(e.from, e.to))
        })
        .cloned()
        .collect()
}

fn qual_to_error_prob(qual: u8) -> f64 {
    10_f64.powf(-(qual as f64) / 10.0)
}

fn error_prob_to_qual(p: f64) -> u8 {
    let q = (-10.0 * p.max(1e-300).log10()).round() as i32;
    q.clamp(0, 254) as u8
}

fn qual_to_log_error_prob(qual: u8) -> f64 {
    let e = qual_to_error_prob(qual);
    e.ln()
}

fn qual_to_log_prob(qual: u8) -> f64 {
    let e = qual_to_error_prob(qual);
    (1.0 - e).ln()
}

fn fast_bernoulli_entropy(z: f64) -> f64 {
    if z <= 0.0 || z >= 1.0 {
        return 0.0;
    }
    let o = 1.0 - z;
    -z * z.ln() - o * o.ln()
}

fn log_binomial(n: usize, k: usize) -> f64 {
    crate::activity_scoring::log_binomial_coefficient_natural(n as u32, k as u32)
}

/// GATK `Mutect2Engine.logLikelihoodRatio(int refCount, int altCount, double errorProbability)`.
fn mutect_log_likelihood_ratio(n_ref: usize, n_alt: usize, error_probability: f64) -> f64 {
    if n_alt == 0 {
        return 0.0;
    }
    let qual = error_prob_to_qual(error_probability);
    let n = n_ref + n_alt;
    let digamma = statrs::function::gamma::digamma;
    let f_tilde_ratio = (digamma((n_ref + 1) as f64) - digamma((n_alt + 1) as f64)).exp();
    let epsilon = qual_to_error_prob(qual);
    let z_bar_alt = (1.0 - epsilon) / (1.0 - epsilon + epsilon * f_tilde_ratio);
    let log_epsilon = qual_to_log_error_prob(qual);
    let log_one_minus_epsilon = qual_to_log_prob(qual);
    let read_sum =
        z_bar_alt * (log_one_minus_epsilon - log_epsilon) + fast_bernoulli_entropy(z_bar_alt);
    let beta_entropy = -(n as f64 + 1.0).ln() - log_binomial(n, n_alt);
    beta_entropy + read_sum
}

fn chain_log_odds(chain: &LinearChain, graph: &AssemblyGraph, error_rate: f64) -> (f64, f64) {
    let left_total: usize = graph
        .outgoing_nodes(chain.first)
        .iter()
        .filter_map(|&t| graph.edge_pruning_support(chain.first, t))
        .map(|s| s as usize)
        .sum();
    let right_total: usize = graph
        .incoming_nodes(chain.last)
        .iter()
        .filter_map(|&f| graph.edge_pruning_support(f, chain.last))
        .map(|s| s as usize)
        .sum();
    let left_mult = chain.edges[0].support as usize;
    let right_mult = chain.edges.last().map(|e| e.support as usize).unwrap_or(0);
    let left = if graph.is_source(chain.first) {
        0.0
    } else {
        mutect_log_likelihood_ratio(left_total.saturating_sub(left_mult), left_mult, error_rate)
    };
    let right = if graph.is_sink(chain.last) {
        0.0
    } else {
        mutect_log_likelihood_ratio(
            right_total.saturating_sub(right_mult),
            right_mult,
            error_rate,
        )
    };
    (left, right)
}

fn max_weight_chain_index(chains: &[LinearChain], graph: &AssemblyGraph) -> usize {
    chains
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            let wa = a.edges.iter().map(|e| e.support).max().unwrap_or(0);
            let wb = b.edges.iter().map(|e| e.support).max().unwrap_or(0);
            wa.cmp(&wb)
                .then_with(|| a.edges.len().cmp(&b.edges.len()))
                .then_with(|| graph.kmer_at(a.first).cmp(graph.kmer_at(b.first)))
                .then_with(|| chain_bases(a, graph).cmp(&chain_bases(b, graph)))
        })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

fn chain_bases(chain: &LinearChain, graph: &AssemblyGraph) -> Vec<u8> {
    if chain.edges.is_empty() {
        return Vec::new();
    }
    graph.kmer_at(chain.first).to_vec()
}

fn adaptive_chains_to_remove(
    chains: &[LinearChain],
    graph: &AssemblyGraph,
    params: &AssemblyGraphPruningParams,
) -> Vec<LinearChain> {
    if chains.is_empty() {
        return Vec::new();
    }
    let pass1 = likely_error_chains(chains, graph, params.initial_error_rate_for_pruning, params);
    let error_count: usize = pass1
        .iter()
        .flat_map(|c| c.edges.iter().map(|e| e.support as usize))
        .sum();
    let total_bases: usize = chains
        .iter()
        .flat_map(|c| c.edges.iter().map(|e| e.support as usize))
        .sum();
    let error_rate = if total_bases == 0 {
        params.initial_error_rate_for_pruning
    } else {
        error_count as f64 / total_bases as f64
    };
    likely_error_chains(chains, graph, error_rate, params)
}

fn likely_error_chains(
    chains: &[LinearChain],
    graph: &AssemblyGraph,
    error_rate: f64,
    params: &AssemblyGraphPruningParams,
) -> Vec<LinearChain> {
    let mut chain_lods: Vec<(usize, f64, f64)> = Vec::with_capacity(chains.len());
    for (i, c) in chains.iter().enumerate() {
        let (left, right) = chain_log_odds(c, graph, error_rate);
        chain_lods.push((i, left, right));
    }

    let mut vertex_to_seedable: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    let mut vertex_to_good_in: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    let mut vertex_to_good_out: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();

    for &(idx, left, right) in &chain_lods {
        let chain = &chains[idx];
        if right >= params.pruning_log_odds_threshold {
            vertex_to_good_in.entry(chain.last).or_default().push(idx);
        }
        if left >= params.pruning_log_odds_threshold {
            vertex_to_good_out.entry(chain.first).or_default().push(idx);
        }
        if left >= params.pruning_seeding_log_odds_threshold
            && right >= params.pruning_seeding_log_odds_threshold
        {
            vertex_to_seedable.entry(chain.first).or_default().push(idx);
            vertex_to_seedable.entry(chain.last).or_default().push(idx);
        }
    }

    #[derive(Eq, PartialEq)]
    struct HeapItem {
        chain_idx: usize,
        score: OrderedFloat,
        tie: String,
    }
    impl Ord for HeapItem {
        fn cmp(&self, other: &Self) -> Ordering {
            other
                .score
                .cmp(&self.score)
                .then_with(|| self.tie.cmp(&other.tie))
        }
    }
    impl PartialOrd for HeapItem {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }
    #[derive(Eq, PartialEq, PartialOrd, Ord)]
    struct OrderedFloat(u64);
    impl OrderedFloat {
        fn new(v: f64) -> Self {
            OrderedFloat(v.to_bits())
        }
    }

    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    let max_idx = max_weight_chain_index(chains, graph);
    heap.push(HeapItem {
        chain_idx: max_idx,
        score: OrderedFloat::new(f64::INFINITY),
        tie: format!("{}:{}", chains[max_idx].first, chains[max_idx].last),
    });

    let mut processed_vertices: HashSet<usize> = HashSet::new();
    for (&v, seed_chains) in &vertex_to_seedable {
        if seed_chains.len() > 2 {
            if let Some(out) = vertex_to_good_out.get(&v) {
                for &idx in out {
                    let c = &chains[idx];
                    let (l, _) = chain_log_odds(c, graph, error_rate);
                    heap.push(HeapItem {
                        chain_idx: idx,
                        score: OrderedFloat::new(l),
                        tie: format!("{}:{}", c.first, c.last),
                    });
                }
            }
            if let Some(inc) = vertex_to_good_in.get(&v) {
                for &idx in inc {
                    let c = &chains[idx];
                    let (_, r) = chain_log_odds(c, graph, error_rate);
                    heap.push(HeapItem {
                        chain_idx: idx,
                        score: OrderedFloat::new(r),
                        tie: format!("{}:{}", c.first, c.last),
                    });
                }
            }
            processed_vertices.insert(v);
        }
    }

    let mut good: HashSet<usize> = HashSet::new();
    let mut vertices_with_outgoing_good: HashSet<usize> = HashSet::new();
    let mut variant_count = 0usize;

    while let Some(item) = heap.pop() {
        let chain = &chains[item.chain_idx];
        if !good.insert(item.chain_idx) {
            continue;
        }
        let new_variant = vertices_with_outgoing_good.insert(chain.first);
        if new_variant {
            variant_count += 1;
        }
        if new_variant && variant_count > params.max_unpruned_variants {
            continue;
        }
        for v in [chain.first, chain.last] {
            if processed_vertices.contains(&v) {
                continue;
            }
            if let Some(out) = vertex_to_good_out.get(&v) {
                for &idx in out {
                    let c = &chains[idx];
                    let (l, _) = chain_log_odds(c, graph, error_rate);
                    heap.push(HeapItem {
                        chain_idx: idx,
                        score: OrderedFloat::new(l),
                        tie: format!("{}:{}", c.first, c.last),
                    });
                }
            }
            if let Some(inc) = vertex_to_good_in.get(&v) {
                for &idx in inc {
                    let c = &chains[idx];
                    let (_, r) = chain_log_odds(c, graph, error_rate);
                    heap.push(HeapItem {
                        chain_idx: idx,
                        score: OrderedFloat::new(r),
                        tie: format!("{}:{}", c.first, c.last),
                    });
                }
            }
            processed_vertices.insert(v);
        }
    }

    chains
        .iter()
        .enumerate()
        .filter_map(|(i, c)| {
            if good.contains(&i) {
                None
            } else {
                Some(c.clone())
            }
        })
        .collect()
}

fn remove_chains(graph: &mut AssemblyGraph, to_remove: &[LinearChain]) {
    for c in to_remove {
        for e in &c.edges {
            graph.remove_edge(e.from, e.to);
        }
    }
    graph.remove_isolated_nodes();
}

/// GATK `LowWeightChainPruner` + optional `AdaptiveChainPruner`.
pub fn apply_gatk_pruning(graph: &mut AssemblyGraph, params: &AssemblyGraphPruningParams) -> u32 {
    let before = graph.edge_count();
    let chains = find_all_chains(graph);
    let to_remove = if params.use_adaptive_pruning {
        adaptive_chains_to_remove(&chains, graph, params)
    } else {
        low_weight_chains_to_remove(&chains, graph, params.min_prune_factor)
    };
    remove_chains(graph, &to_remove);
    before.saturating_sub(graph.edge_count()) as u32
}

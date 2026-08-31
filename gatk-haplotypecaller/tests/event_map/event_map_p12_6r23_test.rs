//! 6R.23 TEST-ONLY: why k=85 C/A threading is a REF-disconnected island.
//! Does not change production assembly, dangling, k, gates, EventMap, or W-H1.

#[cfg(test)]
mod traces {
    use crate::assembly::{
        AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead,
    };
    use crate::assembly_region_finalize::{
        assembly_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, padded_reference_loc, records_to_assembly_reads,
    };
    use crate::bio_ids::KmerSize;
    use crate::kmer_key::{key_from_window, KmerKey, MAX_PACKED128_K};
    use crate::read_projection::query_index_at_reference_position;
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use rust_htslib::bam::record::CigarString;
    use std::collections::HashSet;
    use std::path::Path;

    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const K: usize = 85;
    const MIN_Q: u8 = 10;

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn load_mid_b() -> Option<(
        AssemblyRead,
        Vec<AssemblyRead>,
        Vec<rust_htslib::bam::Record>,
        u64,
    )> {
        let (ref_fasta, bam) = fixture_paths()?;
        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).ok()?;
        let specs = parse_intervals_cli_string(&dict, "2:92317000-92319000").ok()?;
        let walk = crate::walker_traversal::traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            &bam,
            &crate::read_model::ReadFilterParams::gatk_standard_hc(),
            &crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .ok()?;
        let regions = crate::walker_traversal::flatten_assembly_regions(&walk);
        let region = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= SITE_CA
                && r.end.get() >= SITE_CA
        })?;
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).ok()?;
        let (pad_start, _) = padded_reference_loc(region, &dict);
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        Some((reference, assembly_reads, finalized, pad_start))
    }

    fn graph_params() -> AssemblyGraphParams {
        AssemblyGraphParams {
            kmer_size: KmerSize::try_from_usize(K).expect("k"),
            min_base_quality: MIN_Q,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: 128,
            max_haplotype_bases: 4096,
            start_threading_only_at_existing_vertex: false,
        }
    }

    fn hq_segments(read: &AssemblyRead) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut last_good: Option<usize> = None;
        let bytes = read.bases.as_slice();
        for end in 0..=bytes.len() {
            let unusable = end == bytes.len()
                || read.base_quals[end] < MIN_Q
                || !matches!(bytes[end], b'A' | b'C' | b'G' | b'T' | b'N');
            if unusable {
                if let Some(start) = last_good {
                    if end - start >= K {
                        out.push((start, end));
                    }
                }
                last_good = None;
            } else if last_good.is_none() {
                last_good = Some(end);
            }
        }
        out
    }

    fn ref_kmer_set(ref_bases: &[u8]) -> HashSet<Vec<u8>> {
        let mut s = HashSet::new();
        if ref_bases.len() < K {
            return s;
        }
        for i in 0..=ref_bases.len() - K {
            s.insert(ref_bases[i..i + K].to_vec());
        }
        s
    }

    fn non_unique_in_seq(bases: &[u8]) -> HashSet<KmerKey> {
        let mut seen = HashSet::new();
        let mut non = HashSet::new();
        if bases.len() < K {
            return non;
        }
        for i in 0..=bases.len() - K {
            let key = key_from_window(bases, i, K);
            if !seen.insert(key.clone()) {
                non.insert(key);
            }
        }
        non
    }

    fn find_start_java(bases: &[u8], non_unique: &HashSet<KmerKey>) -> Option<usize> {
        if bases.len() < K {
            return None;
        }
        let last = bases.len().saturating_sub(K);
        for i in 0..last {
            let key = key_from_window(bases, i, K);
            if !non_unique.contains(&key) {
                return Some(i);
            }
        }
        None
    }

    fn compact(k: &[u8]) -> String {
        let s = String::from_utf8_lossy(k);
        if s.len() <= 16 {
            s.into_owned()
        } else {
            format!("{}..{}", &s[..8], &s[s.len() - 8..])
        }
    }

    fn revcomp(bases: &[u8]) -> Vec<u8> {
        bases
            .iter()
            .rev()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                x => x,
            })
            .collect()
    }

    fn undirected_component(graph: &AssemblyGraph, seeds: &[usize]) -> HashSet<usize> {
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); graph.node_count()];
        for e in graph.edges_sorted() {
            adj[e.from].push(e.to);
            adj[e.to].push(e.from);
        }
        let mut seen = HashSet::new();
        let mut stack = Vec::new();
        for &s in seeds {
            if seen.insert(s) {
                stack.push(s);
            }
        }
        while let Some(v) = stack.pop() {
            for &w in &adj[v] {
                if seen.insert(w) {
                    stack.push(w);
                }
            }
        }
        seen
    }

    fn lookup_vertex(graph: &AssemblyGraph, kmer: &[u8]) -> Option<usize> {
        if let Some(id) = graph.vertex_id_for_kmer(kmer) {
            return Some(id);
        }
        graph
            .nodes()
            .iter()
            .find(|n| n.kmer.as_ref() == kmer)
            .map(|n| n.id)
    }

    fn synth_read(ref_bases: &[u8], subs: &[(usize, u8)]) -> AssemblyRead {
        let mut bases = ref_bases.to_vec();
        for &(off, b) in subs {
            if off < bases.len() {
                bases[off] = b;
            }
        }
        AssemblyRead {
            base_quals: vec![30; bases.len()],
            bases,
        }
    }

    struct ControlTopo {
        label: &'static str,
        n_read_kmers: usize,
        n_identical_to_ref: usize,
        n_alt_only: usize,
        first_alt: Option<usize>,
        last_alt: Option<usize>,
        island_vertices: usize,
        island_ref_nodes: usize,
        connected: bool,
    }

    fn classify_control(
        label: &'static str,
        reference: &AssemblyRead,
        reads: &[AssemblyRead],
        seed_kmers: &[Vec<u8>],
    ) -> ControlTopo {
        let (graph, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
            reference,
            reads,
            &graph_params(),
        )
        .expect("control graph");
        let ref_set = ref_kmer_set(&reference.bases);
        let mut n_kmers = 0usize;
        let mut n_id = 0usize;
        let mut n_alt = 0usize;
        let mut first_alt = None;
        let mut last_alt = None;
        for read in reads {
            if read.bases.len() < K {
                continue;
            }
            for i in 0..=read.bases.len() - K {
                n_kmers += 1;
                let km = &read.bases[i..i + K];
                if ref_set.contains(km) {
                    n_id += 1;
                } else {
                    n_alt += 1;
                    if first_alt.is_none() {
                        first_alt = Some(i);
                    }
                    last_alt = Some(i);
                }
            }
        }
        let mut seeds = Vec::new();
        for km in seed_kmers {
            if let Some(id) = lookup_vertex(&graph, km) {
                seeds.push(id);
            }
        }
        let island_vertices;
        let island_ref;
        let connected;
        if seeds.is_empty() {
            island_vertices = 0;
            island_ref = 0;
            connected = false;
        } else {
            let comp = undirected_component(&graph, &seeds);
            island_vertices = comp.len();
            island_ref = comp.iter().filter(|v| graph.ref_nodes.contains(v)).count();
            connected = island_ref > 0;
        }
        ControlTopo {
            label,
            n_read_kmers: n_kmers,
            n_identical_to_ref: n_id,
            n_alt_only: n_alt,
            first_alt,
            last_alt,
            island_vertices,
            island_ref_nodes: island_ref,
            connected,
        }
    }

    fn seed_kmers_overlapping(bases: &[u8], offs: &[usize]) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        if bases.len() < K {
            return out;
        }
        for &off in offs {
            if off >= bases.len() {
                continue;
            }
            let first = off.saturating_sub(K - 1);
            let last = off.min(bases.len() - K);
            if first > last {
                continue;
            }
            for st in first..=last {
                out.push(bases[st..st + K].to_vec());
            }
        }
        out.sort();
        out.dedup();
        out
    }

    #[test]
    fn six_r23_k85_alt_island_threading_audit() {
        let Some((reference, assembly_reads, finalized, pad_start)) = load_mid_b() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        assert_eq!(assembly_reads.len(), 2);
        let ref_bases = reference.bases.as_slice();
        let off_ca = (SITE_CA.saturating_sub(pad_start)) as usize;
        let off_tc = (SITE_TC.saturating_sub(pad_start)) as usize;
        let off_gc = (SITE_GC.saturating_sub(pad_start)) as usize;
        eprintln!(
            "PAD start={pad_start} len={} off_CA={off_ca} off_TC={off_tc} off_GC={off_gc} REF@CA={} @TC={} @GC={}",
            ref_bases.len(),
            ref_bases.get(off_ca).copied().unwrap_or(b'?') as char,
            ref_bases.get(off_tc).copied().unwrap_or(b'?') as char,
            ref_bases.get(off_gc).copied().unwrap_or(b'?') as char
        );

        let params = graph_params();
        let (mut graph, summary) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &reference,
            &assembly_reads,
            &params,
        )
        .expect("k=85 graph");
        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = 2;
        let pruned = graph.apply_pruning(&pruning);
        eprintln!(
            "RAW/PRUNE k={K} nodes={} edges={} low_complexity={} unique_kmers={} non_unique={} prune_removed={pruned}",
            graph.node_count(),
            graph.edge_count(),
            summary.is_low_complexity,
            summary.unique_kmer_count,
            summary.non_unique_kmer_count
        );
        assert!(graph.node_count() > 0);
        assert_eq!(pruned, 0);

        let ref_set = ref_kmer_set(ref_bases);
        eprintln!(
            "KMER_IDENTITY k={K} packed128_max={MAX_PACKED128_K} uses_bytes_keys={} ref_unique_85mers={}",
            K > MAX_PACKED128_K,
            ref_set.len()
        );

        let mut all_non_unique = HashSet::new();
        all_non_unique.extend(non_unique_in_seq(ref_bases));
        for ar in &assembly_reads {
            for (s, e) in hq_segments(ar) {
                all_non_unique.extend(non_unique_in_seq(&ar.bases[s..e]));
            }
        }
        eprintln!(
            "NON_UNIQUE_KEYS (preprocess-style, HQ slices + REF) n={}",
            all_non_unique.len()
        );

        let mut seed_ids = Vec::new();
        let mut n_read_kmers = 0usize;
        let mut n_identical = 0usize;
        let mut n_alt_only = 0usize;
        let mut n_rc_hits_ref = 0usize;
        let mut n_key_mismatch_same_bytes = 0usize;
        let mut first_alt_global: Option<(usize, usize)> = None;
        let mut last_alt_global: Option<(usize, usize)> = None;

        for (ri, ar) in assembly_reads.iter().enumerate() {
            let cigar = CigarString(finalized[ri].cigar().iter().copied().collect());
            let qname = String::from_utf8_lossy(finalized[ri].qname()).into_owned();
            let segs = hq_segments(ar);
            eprintln!(
                "READ[{ri}] qname={qname} aln_start1={} seq_len={} hq_segments={:?} cigar={cigar}",
                finalized[ri].pos() + 1,
                ar.bases.len(),
                segs
            );
            let mut mismatches = Vec::new();
            for qi in 0..ar.bases.len() {
                let Some(rp) = crate::read_projection::reference_position_at_query_index(
                    finalized[ri].pos(),
                    &cigar,
                    qi,
                ) else {
                    continue;
                };
                if rp < 0 {
                    continue;
                }
                let gpos = rp as u64 + 1;
                if gpos < pad_start {
                    continue;
                }
                let off = (gpos - pad_start) as usize;
                let Some(&rb) = ref_bases.get(off) else {
                    continue;
                };
                let b = ar.bases[qi];
                if b != rb && matches!(b, b'A' | b'C' | b'G' | b'T') {
                    mismatches.push((gpos, qi, rb as char, b as char));
                }
            }
            eprintln!(
                "  MISMATCHES_VS_PADDED_REF n={} max_match_run={} {:?}",
                mismatches.len(),
                {
                    let mut qis: Vec<usize> = mismatches.iter().map(|t| t.1).collect();
                    qis.sort_unstable();
                    let mut max_run = if qis.is_empty() {
                        ar.bases.len()
                    } else {
                        qis[0]
                    };
                    for w in qis.windows(2) {
                        max_run = max_run.max(w[1].saturating_sub(w[0]).saturating_sub(1));
                    }
                    if let Some(&last) = qis.last() {
                        max_run = max_run.max(ar.bases.len().saturating_sub(last + 1));
                    }
                    max_run
                },
                mismatches
            );
            if ar.bases.len() >= K {
                let aln0 = finalized[ri].pos();
                let ref_off = (aln0 as u64 + 1).saturating_sub(pad_start) as usize;
                if ref_off + K <= ref_bases.len() {
                    let rk = &ref_bases[ref_off..ref_off + K];
                    let ak = &ar.bases[..K];
                    let n_diff = rk.iter().zip(ak.iter()).filter(|(a, b)| a != b).count();
                    eprintln!(
                        "  FIRST_85MER vs REF@{} diffs={n_diff} read={} REF={} keys_eq={}",
                        pad_start + ref_off as u64,
                        compact(ak),
                        compact(rk),
                        key_from_window(ak, 0, K) == key_from_window(rk, 0, K)
                    );
                }
            }
            for site in [SITE_CA, SITE_TC, SITE_GC] {
                let qi = query_index_at_reference_position(
                    finalized[ri].pos(),
                    &cigar,
                    (site - 1) as i64,
                );
                match qi {
                    Some(i) if i < ar.bases.len() => {
                        let rb = {
                            let off = (site.saturating_sub(pad_start)) as usize;
                            ref_bases.get(off).copied().unwrap_or(b'?')
                        };
                        eprintln!(
                            "  SITE {site} qi={i} read={} q={} REF={}",
                            ar.bases[i] as char,
                            ar.base_quals.get(i).copied().unwrap_or(0),
                            rb as char
                        );
                    }
                    _ => eprintln!("  SITE {site} no query base"),
                }
            }

            for (s, e) in &segs {
                let bases = &ar.bases[*s..*e];
                let start = find_start_java(bases, &all_non_unique);
                let start_kmer = start.map(|i| bases[i..i + K].to_vec());
                let start_in_ref = start_kmer
                    .as_ref()
                    .map(|km| ref_set.contains(km))
                    .unwrap_or(false);
                let start_is_ref_node = start_kmer
                    .as_ref()
                    .and_then(|km| lookup_vertex(&graph, km))
                    .map(|id| graph.ref_nodes.contains(&id))
                    .unwrap_or(false);
                eprintln!(
                    "  FIND_START seg={s}..{e} start={start:?} start_in_ref_set={start_in_ref} start_is_ref_node={start_is_ref_node} kmer={}",
                    start_kmer.as_deref().map(compact).unwrap_or_else(|| "-".into())
                );

                if bases.len() < K {
                    continue;
                }
                let mut trace: Vec<(&'static str, usize)> = Vec::new();
                let mut prev_class = "";
                for i in 0..=bases.len() - K {
                    n_read_kmers += 1;
                    let km = &bases[i..i + K];
                    let key = key_from_window(bases, i, K);
                    let key2 = key_from_window(km, 0, K);
                    if key != key2 {
                        n_key_mismatch_same_bytes += 1;
                    }
                    if ref_set.contains(&revcomp(km)) {
                        n_rc_hits_ref += 1;
                    }
                    let in_ref = ref_set.contains(km);
                    if in_ref {
                        n_identical += 1;
                    } else {
                        n_alt_only += 1;
                        if first_alt_global.is_none() {
                            first_alt_global = Some((ri, i));
                        }
                        last_alt_global = Some((ri, i));
                    }
                    let vid = lookup_vertex(&graph, km);
                    let is_ref_node = vid.map(|id| graph.ref_nodes.contains(&id)).unwrap_or(false);
                    let class = if is_ref_node {
                        "REF"
                    } else if in_ref {
                        "REFSEQ_NOT_NODE"
                    } else if vid.is_some() {
                        "ALT"
                    } else {
                        "ABSENT"
                    };
                    if class != prev_class {
                        trace.push((class, i));
                        prev_class = class;
                    }
                    if !in_ref {
                        if let Some(id) = vid {
                            seed_ids.push(id);
                        }
                    }
                }
                eprintln!("  THREAD_TRACE_TRANSITIONS (class,kmer_start)={trace:?}");
                if let Some(st) = start {
                    let mut path = String::new();
                    let mut prev: Option<usize> = None;
                    let mut n_ref_v = 0usize;
                    let mut n_alt_v = 0usize;
                    let last = bases.len() - K;
                    for i in st..=last {
                        let km = &bases[i..i + K];
                        let vid = lookup_vertex(&graph, km);
                        let tag = match vid {
                            Some(id) if graph.ref_nodes.contains(&id) => {
                                n_ref_v += 1;
                                "R"
                            }
                            Some(_) => {
                                n_alt_v += 1;
                                "A"
                            }
                            None => {
                                n_alt_v += 1;
                                "?"
                            }
                        };
                        if i == st
                            || i == last
                            || (i > st && {
                                let prev_km = &bases[i - 1..i - 1 + K];
                                let prev_ref = lookup_vertex(&graph, prev_km)
                                    .map(|id| graph.ref_nodes.contains(&id))
                                    .unwrap_or(false);
                                let cur_ref = tag == "R";
                                prev_ref != cur_ref
                            })
                        {
                            if !path.is_empty() {
                                path.push_str(" -> ");
                            }
                            path.push_str(&format!("{tag}@{i}"));
                        }
                        if let (Some(p), Some(id)) = (prev, vid) {
                            let _ = graph.edge_support(p, id);
                        }
                        prev = vid;
                    }
                    eprintln!(
                        "  THREAD_PATH start={st}..{last} n_ref_vertices={n_ref_v} n_nonref_vertices={n_alt_v} spine={path}"
                    );
                }
            }
        }

        seed_ids.sort_unstable();
        seed_ids.dedup();
        eprintln!(
            "READ_KMER_VS_REF n={n_read_kmers} identical={n_identical} alt_only={n_alt_only} first_alt={first_alt_global:?} last_alt={last_alt_global:?} rc_hits_ref={n_rc_hits_ref} key_mismatch_same_bytes={n_key_mismatch_same_bytes}"
        );
        assert_eq!(
            n_key_mismatch_same_bytes, 0,
            "byte-identical windows must hash-equal as KmerKey"
        );

        assert!(
            !seed_ids.is_empty(),
            "ALT-only read 85-mers must exist in the RAW graph"
        );
        let component = undirected_component(&graph, &seed_ids);
        let n_ref_in = component
            .iter()
            .filter(|v| graph.ref_nodes.contains(v))
            .count();
        let mut ordered: Vec<usize> = component.iter().copied().collect();
        ordered.sort_unstable();
        eprintln!(
            "CA_COMPONENT vertices={} ref_nodes={n_ref_in} seed_alt_vertices={}",
            component.len(),
            seed_ids.len()
        );

        let mut n_in_ref_set = 0usize;
        for &v in &ordered {
            let km = graph.kmer_at(v);
            let in_ref_set = ref_set.contains(km);
            if in_ref_set {
                n_in_ref_set += 1;
            }
            let ins = graph.incoming_nodes(v);
            let outs = graph.outgoing_nodes(v);
            let in_w: Vec<u32> = ins
                .iter()
                .map(|&p| graph.edge_support(p, v).unwrap_or(0))
                .collect();
            let out_w: Vec<u32> = outs
                .iter()
                .map(|&t| graph.edge_support(v, t).unwrap_or(0))
                .collect();
            eprintln!(
                "  V{v} kmer={} in_ref_set={in_ref_set} is_ref_node={} in={:?} w={in_w:?} out={:?} w={out_w:?} seq={}",
                compact(km),
                graph.ref_nodes.contains(&v),
                ins,
                outs,
                String::from_utf8_lossy(km)
            );
        }
        eprintln!(
            "COMPONENT_KMERS_IN_REF_SET={n_in_ref_set}/{n}",
            n = component.len()
        );
        assert_eq!(n_ref_in, 0, "C/A component must have 0 REF nodes (6R.22)");
        assert_eq!(
            n_in_ref_set, 0,
            "no ALT-component 85-mer equals any padded-REF 85-mer"
        );

        let expected_before = off_ca.checked_sub(K);
        let expected_after = off_gc + 1;
        for (label, start) in [
            ("REF_BEFORE_CA", expected_before),
            (
                "REF_AT_CA",
                off_ca
                    .checked_sub(K - 1)
                    .filter(|&s| s + K <= ref_bases.len()),
            ),
            (
                "REF_AFTER_GC",
                Some(expected_after).filter(|&s| s + K <= ref_bases.len()),
            ),
        ] {
            match start {
                Some(st) if st + K <= ref_bases.len() => {
                    let km = &ref_bases[st..st + K];
                    let vid = lookup_vertex(&graph, km);
                    eprintln!(
                        "BOUNDARY {label} ref_off={st} in_graph={} is_ref_node={} kmer={}",
                        vid.is_some(),
                        vid.map(|id| graph.ref_nodes.contains(&id)).unwrap_or(false),
                        compact(km)
                    );
                }
                _ => eprintln!("BOUNDARY {label} ref window OOB"),
            }
        }

        if let Some((_, first_i)) = first_alt_global {
            if first_i > 0 {
                eprintln!(
                    "BOUNDARY first_alt_kmer_start={first_i} previous_kmer_would_start={}",
                    first_i - 1
                );
            } else {
                eprintln!(
                    "BOUNDARY first_alt_kmer_start=0 (HQ segment begins already ALT-only; no left REF kmer on the read)"
                );
            }
        }

        let mut ca_only = Vec::new();
        if off_ca < ref_bases.len() {
            let first = off_ca.saturating_sub(K - 1);
            let last = off_ca.min(ref_bases.len() - K);
            for st in first..=last {
                let mut km = ref_bases[st..st + K].to_vec();
                km[off_ca - st] = b'A';
                ca_only.push(km);
            }
        }

        let c_ca = classify_control(
            "synth_CA",
            &reference,
            &[synth_read(ref_bases, &[(off_ca, b'A')])],
            &ca_only,
        );
        let c_ca_tc = classify_control(
            "synth_CA_TC",
            &reference,
            &[synth_read(ref_bases, &[(off_ca, b'A'), (off_tc, b'C')])],
            &seed_kmers_overlapping(
                &synth_read(ref_bases, &[(off_ca, b'A'), (off_tc, b'C')]).bases,
                &[off_ca, off_tc],
            ),
        );
        let trip_read = synth_read(ref_bases, &[(off_ca, b'A'), (off_tc, b'C'), (off_gc, b'C')]);
        let c_trip = classify_control(
            "synth_CA_TC_GC",
            &reference,
            &[trip_read.clone()],
            &seed_kmers_overlapping(&trip_read.bases, &[off_ca, off_tc, off_gc]),
        );
        let real_seeds: Vec<Vec<u8>> = {
            let mut s = Vec::new();
            for ar in &assembly_reads {
                if ar.bases.len() < K {
                    continue;
                }
                for i in 0..=ar.bases.len() - K {
                    let km = ar.bases[i..i + K].to_vec();
                    if !ref_set.contains(&km) {
                        s.push(km);
                    }
                }
            }
            s.sort();
            s.dedup();
            s
        };
        let c_real = classify_control("real_reads", &reference, &assembly_reads, &real_seeds);

        for c in [&c_ca, &c_ca_tc, &c_trip, &c_real] {
            eprintln!(
                "CONTROL[{}] read_kmers={} identical_to_REF={} alt_only={} first_alt={:?} last_alt={:?} island_n={} island_ref_nodes={} REF_connected={}",
                c.label,
                c.n_read_kmers,
                c.n_identical_to_ref,
                c.n_alt_only,
                c.first_alt,
                c.last_alt,
                c.island_vertices,
                c.island_ref_nodes,
                c.connected
            );
        }

        assert!(
            c_ca.n_identical_to_ref > 0,
            "full-REF synthetic C→A must share flanking 85-mers with REF"
        );
        assert!(
            c_ca.connected,
            "lone C→A on full padded REF must form a REF-connected bubble at k=85"
        );
        assert!(
            c_trip.connected,
            "C→A+T→C+G→C on full padded REF must remain REF-connected at k=85"
        );
        assert_eq!(n_identical, 0);
        assert_eq!(n_rc_hits_ref, 0);
        assert!(
            !c_real.connected,
            "real 2-read haplotype remains a REF-free island (6R.22 topology)"
        );
    }
}

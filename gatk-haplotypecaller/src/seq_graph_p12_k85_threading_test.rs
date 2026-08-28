//! 6R.11 TEST-ONLY: why k=85 threading leaves the P12 TATG island disconnected.
//! Does not change production assembler behavior, the P12 waiver, or W-H1.

#[cfg(test)]
mod traces {
    use crate::assembly::{AssemblyGraph, AssemblyGraphParams, AssemblyRead};
    use crate::read_event_discovery::P12_CLUSTER_TTC_START;
    use crate::read_projection::query_index_at_reference_position;
    use crate::read_threading_assembler::ReadThreadingAssemblerArgs;
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
    use rust_htslib::bam;
    use rust_htslib::bam::record::CigarString;
    use std::collections::{HashMap, HashSet};
    use std::path::Path;

    const K: usize = 85;
    const ALT_CORE: &[u8] = b"TATGTG";
    const REF_CORE: &[u8] = b"TTCATG";
    const ALT_WIN: &[u8] = b"CTTTTATGTGATGTAT";
    const REF_WIN: &[u8] = b"CTTTTTCATGATGTAT";
    const MIN_Q: u8 = 10;

    const REAL_P12_ACTIVE_START: u64 = P12_CLUSTER_TTC_START - 96;
    const REAL_P12_ACTIVE_END: u64 = P12_CLUSTER_TTC_START + 76;
    const REAL_P12_ATG_START: u64 = P12_CLUSTER_TTC_START + 3;
    const REAL_P12_WIN_LO: u64 = P12_CLUSTER_TTC_START - 4;
    const REAL_P12_WIN_HI: u64 = P12_CLUSTER_TTC_START + 11;

    fn k85_params(start_only_existing: bool) -> AssemblyGraphParams {
        AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_from_usize(K).expect("k=85"),
            min_base_quality: MIN_Q,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: 32,
            max_haplotype_bases: 4096,
            start_threading_only_at_existing_vertex: start_only_existing,
        }
    }

    fn unique_dna(len: usize, salt: u64) -> String {
        const BASES: &[u8] = b"ACGT";
        let mut s = String::with_capacity(len);
        let mut x = salt;
        for _ in 0..len {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            s.push(BASES[((x >> 62) as usize) % 4] as char);
        }
        s
    }

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    fn usable_segments(bases: &[u8], quals: &[u8], k: usize, min_q: u8) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut last_good: Option<usize> = None;
        for end in 0..=bases.len() {
            let unusable = end == bases.len()
                || quals[end] < min_q
                || !matches!(bases[end], b'A' | b'C' | b'G' | b'T' | b'N');
            if unusable {
                if let Some(start) = last_good {
                    if end - start >= k {
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

    fn kmers_in(seq: &[u8], k: usize) -> HashSet<Vec<u8>> {
        if seq.len() < k {
            return HashSet::new();
        }
        (0..=seq.len() - k)
            .map(|i| seq[i..i + k].to_vec())
            .collect()
    }

    fn non_unique_in(seq: &[u8], k: usize) -> HashSet<Vec<u8>> {
        let mut seen = HashSet::new();
        let mut dup = HashSet::new();
        if seq.len() < k {
            return dup;
        }
        for i in 0..=seq.len() - k {
            let w = seq[i..i + k].to_vec();
            if !seen.insert(w.clone()) {
                dup.insert(w);
            }
        }
        dup
    }

    /// GATK `findStart`: `i < stop - kmerSize` (last k-mer cannot be the start).
    fn find_start(
        seq: &[u8],
        k: usize,
        only_existing: bool,
        unique_in_graph: &HashSet<Vec<u8>>,
        non_unique: &HashSet<Vec<u8>>,
    ) -> Option<usize> {
        let last = seq.len().saturating_sub(k);
        for i in 0..last {
            let kmer = seq[i..i + k].to_vec();
            let ok = if only_existing {
                unique_in_graph.contains(&kmer)
            } else {
                !non_unique.contains(&kmer)
            };
            if ok {
                return Some(i);
            }
        }
        None
    }

    fn rt_cc(g: &AssemblyGraph) -> usize {
        let n = g.node_count();
        if n == 0 {
            return 0;
        }
        let mut adj = vec![Vec::new(); n];
        for e in g.edges_sorted() {
            adj[e.from].push(e.to);
            adj[e.to].push(e.from);
        }
        let mut seen = vec![false; n];
        let mut c = 0usize;
        for i in 0..n {
            if seen[i] {
                continue;
            }
            c += 1;
            let mut stack = vec![i];
            seen[i] = true;
            while let Some(v) = stack.pop() {
                for &w in &adj[v] {
                    if !seen[w] {
                        seen[w] = true;
                        stack.push(w);
                    }
                }
            }
        }
        c
    }

    fn vertex_for(g: &AssemblyGraph, kmer: &[u8]) -> Option<usize> {
        g.vertex_id_for_kmer(kmer)
            .or_else(|| g.nodes().iter().position(|n| n.kmer.as_ref() == kmer))
    }

    fn component_of(g: &AssemblyGraph, start: usize) -> HashSet<usize> {
        let mut seen = HashSet::new();
        let mut stack = vec![start];
        seen.insert(start);
        while let Some(v) = stack.pop() {
            for w in g.outgoing_nodes(v).into_iter().chain(g.incoming_nodes(v)) {
                if seen.insert(w) {
                    stack.push(w);
                }
            }
        }
        seen
    }

    fn format_cigar(rec: &bam::Record) -> String {
        rec.cigar()
            .iter()
            .map(|c| format!("{c}"))
            .collect::<String>()
    }

    fn identical_run_left(
        rec: &bam::Record,
        pad_start: u64,
        first_var_1based: u64,
        ref_bytes: &[u8],
    ) -> usize {
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let seq = rec.seq().as_bytes();
        let mut n = 0usize;
        let mut rp = first_var_1based;
        loop {
            if rp <= pad_start {
                break;
            }
            rp -= 1;
            let ref_idx = (rp - pad_start) as usize;
            if ref_idx >= ref_bytes.len() {
                break;
            }
            let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, (rp - 1) as i64)
            else {
                break;
            };
            if qi >= seq.len() || seq[qi] != ref_bytes[ref_idx] {
                break;
            }
            n += 1;
            if n > 2000 {
                break;
            }
        }
        n
    }

    fn identical_run_right(
        rec: &bam::Record,
        pad_start: u64,
        last_var_1based: u64,
        ref_bytes: &[u8],
    ) -> usize {
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let seq = rec.seq().as_bytes();
        let mut n = 0usize;
        let mut rp = last_var_1based;
        loop {
            rp += 1;
            let ref_idx = (rp - pad_start) as usize;
            if ref_idx >= ref_bytes.len() {
                break;
            }
            let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, (rp - 1) as i64)
            else {
                break;
            };
            if qi >= seq.len() || seq[qi] != ref_bytes[ref_idx] {
                break;
            }
            n += 1;
            if n > 2000 {
                break;
            }
        }
        n
    }

    struct LoadedP12 {
        reference: AssemblyRead,
        reads: Vec<AssemblyRead>,
        ref_bytes: Vec<u8>,
        finalized: Vec<bam::Record>,
        pad_start: u64,
    }

    fn load_real_p12() -> Option<LoadedP12> {
        use crate::assembly_region_finalize::{
            assembly_reference_read, finalize_region_reads_for_assembly,
            gatk_min_tail_quality_for_assembly, padded_reference_loc, records_to_assembly_reads,
        };
        use crate::read_model::ReadFilterParams;
        use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
        use crate::walker_traversal::{
            flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
        };
        use gatk_core::reference::{
            parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
        };

        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        let dict = SequenceDictionary::from_fasta_path(&ref_path).ok()?;
        let interval = format!("2:{REAL_P12_ACTIVE_START}-{REAL_P12_ACTIVE_END}");
        let specs = parse_intervals_cli_string(&dict, &interval).ok()?;
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_path,
            &bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .ok()?;
        let regions = flatten_assembly_regions(&walk);
        let region = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= P12_CLUSTER_TTC_START
                && r.end.get() >= REAL_P12_ATG_START
        })?;
        let mut ref_cache = ReferenceWindowCache::new(ref_path.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).ok()?;
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let reads = records_to_assembly_reads(&finalized);
        let (pad_start, _) = padded_reference_loc(region, &dict);
        let ref_bytes = reference.bases.clone();
        Some(LoadedP12 {
            reference,
            reads,
            ref_bytes,
            finalized,
            pad_start,
        })
    }

    #[test]
    fn six_r11_production_pins_unchanged() {
        let abc = include_str!("assembly_based_caller.rs");
        assert!(abc.contains("assembler.use_seq_graph = false;"));
        assert!(abc.contains("assembler.dangling_java_exact = true;"));
        assert_eq!(
            ReadThreadingAssemblerArgs::default().recover_dangling_branches,
            true
        );
        assert!(
            !AssemblyGraphParams::default().start_threading_only_at_existing_vertex,
            "production recover=true => start_threading_only_at_existing_vertex=false"
        );
    }

    #[test]
    fn six_r11_synthetic_shared_85mer_connects_at_k85() {
        let left = unique_dna(120, 7);
        let right = unique_dna(120, 11);
        let ref_seq = format!("{left}TTCA{right}");
        let alt_seq = format!("{left}TATG{right}");
        let reference = read(&ref_seq, 30);
        assert!(
            !crate::read_threading_graph::reference_has_non_unique_kmers(&reference, K),
            "synthetic REF must have unique 85-mers"
        );
        let shared: HashSet<_> = kmers_in(ref_seq.as_bytes(), K)
            .intersection(&kmers_in(alt_seq.as_bytes(), K))
            .cloned()
            .collect();
        assert!(
            !shared.is_empty(),
            "control: unique flanks must produce shared 85-mers"
        );
        let alt_core_start = left.len();
        let overlapping_alt: Vec<Vec<u8>> = (0..=alt_seq.len() - K)
            .filter(|&i| i < alt_core_start + 4 && i + K > alt_core_start)
            .map(|i| alt_seq.as_bytes()[i..i + K].to_vec())
            .collect();
        assert!(!overlapping_alt.is_empty());
        for only_existing in [false, true] {
            let (g, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
                &reference,
                &[read(&alt_seq, 30)],
                &k85_params(only_existing),
            )
            .expect("thread");
            let src = g.reference_source_vertex().expect("src");
            let src_cc = component_of(&g, src);
            let mut attached = 0usize;
            for km in &overlapping_alt {
                let Some(id) = vertex_for(&g, km) else {
                    panic!(
                        "overlapping ALT 85-mer missing from graph only_existing={only_existing}"
                    );
                };
                if src_cc.contains(&id) {
                    attached += 1;
                }
            }
            assert_eq!(rt_cc(&g), 1, "shared flanks must yield cc=1");
            assert!(
                attached == overlapping_alt.len(),
                "shared 85-mers must attach the TTCA/TATG bubble to REF only_existing={only_existing} attached={attached}/{}",
                overlapping_alt.len()
            );
        }
    }

    #[test]
    fn six_r11_synthetic_no_shared_85mer_stays_disconnected() {
        let left = unique_dna(120, 13);
        let right = unique_dna(120, 17);
        let ref_seq = format!("{left}TTCA{right}");
        let alt_left = unique_dna(50, 19);
        let alt_right = unique_dna(50, 23);
        let alt_seq = format!("{alt_left}TATG{alt_right}");
        assert!(alt_seq.len() >= K);
        let shared: HashSet<_> = kmers_in(ref_seq.as_bytes(), K)
            .intersection(&kmers_in(alt_seq.as_bytes(), K))
            .cloned()
            .collect();
        assert!(
            shared.is_empty(),
            "control: unrelated ALT DNA must share no 85-mer"
        );
        let reference = read(&ref_seq, 30);
        let (g_false, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &reference,
            &[read(&alt_seq, 30)],
            &k85_params(false),
        )
        .expect("thread false");
        let (g_true, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &reference,
            &[read(&alt_seq, 30)],
            &k85_params(true),
        )
        .expect("thread true");
        let src_f = g_false.reference_source_vertex().expect("src");
        let src_cc = component_of(&g_false, src_f);
        let tatg_false = g_false.nodes().iter().find(|n| {
            n.kmer
                .as_ref()
                .windows(ALT_CORE.len())
                .any(|w| w == ALT_CORE)
        });
        assert!(tatg_false.is_some(), "false mode should create ALT island");
        assert!(
            !src_cc.contains(&tatg_false.unwrap().id),
            "no shared 85-mer: TATG must stay off the REF component"
        );
        assert!(rt_cc(&g_false) >= 2);
        let tatg_true = g_true.nodes().iter().any(|n| {
            n.kmer
                .as_ref()
                .windows(ALT_CORE.len())
                .any(|w| w == ALT_CORE)
        });
        assert!(
            !tatg_true,
            "start_only_existing=true must drop a read with no graph k-mer"
        );
        assert_eq!(rt_cc(&g_true), 1);
    }

    #[test]
    fn six_r11_real_p12_k85_shared_kmer_and_threading_ab() {
        let Some(p12) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let k = K;
        let ref_kmers = kmers_in(&p12.ref_bytes, k);
        let mut non_unique = non_unique_in(&p12.ref_bytes, k);
        for r in &p12.reads {
            for (a, b) in usable_segments(&r.bases, &r.base_quals, k, MIN_Q) {
                non_unique.extend(non_unique_in(&r.bases[a..b], k));
            }
        }
        let unique_ref: HashSet<Vec<u8>> = ref_kmers.difference(&non_unique).cloned().collect();

        eprintln!("=== 6R.11 REAL P12 k=85 read-threading diagnosis ===");
        eprintln!(
            "pad_start={} ref_len={} reads={} ref_kmers={} unique_ref_kmers={} non_unique={}",
            p12.pad_start,
            p12.ref_bytes.len(),
            p12.finalized.len(),
            ref_kmers.len(),
            unique_ref.len(),
            non_unique.len()
        );
        let off = (P12_CLUSTER_TTC_START - p12.pad_start) as usize;
        let lo = off.saturating_sub(4);
        let hi = (off + 12).min(p12.ref_bytes.len());
        eprintln!(
            "REF_CLUSTER_WINDOW {} (offset {off})",
            String::from_utf8_lossy(&p12.ref_bytes[lo..hi])
        );
        eprintln!(
            "left_85mer_requires_read_start<={}",
            P12_CLUSTER_TTC_START.saturating_sub(K as u64)
        );

        let mut any_shared_from_tatg_read = false;
        let mut tatg_read_count = 0usize;
        let mut ttca_read_count = 0usize;
        let mut island_head_from: Option<String> = None;

        for rec in &p12.finalized {
            let seq = rec.seq().as_bytes();
            let quals = rec.qual().to_vec();
            let cigar = CigarString(rec.cigar().iter().copied().collect());
            let qname = String::from_utf8_lossy(rec.qname()).into_owned();
            let aln_start_1 = rec.pos() + 1;
            let win = window_at_cluster(rec, &cigar, &seq);
            let has_tatg = seq.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
                || win
                    .as_ref()
                    .is_some_and(|w| w.windows(4).any(|x| x == b"TATG"));
            let has_ttca = seq.windows(REF_CORE.len()).any(|w| w == REF_CORE)
                || win
                    .as_ref()
                    .is_some_and(|w| w.windows(4).any(|x| x == b"TTCA"));
            if has_ttca {
                ttca_read_count += 1;
            }
            if has_ttca && !has_tatg {
                let segs = usable_segments(&seq, &quals, k, MIN_Q);
                let shared = segs
                    .iter()
                    .map(|&(a, b)| kmers_in(&seq[a..b], k).intersection(&ref_kmers).count())
                    .sum::<usize>();
                eprintln!(
                    "TTCA_READ qname={qname} start={aln_start_1} cigar={} len={} window={} shared_kmers={shared}",
                    format_cigar(rec),
                    seq.len(),
                    win.as_ref().map(|w| String::from_utf8_lossy(w).into_owned()).unwrap_or_else(|| "-".into()),
                );
            }
            if !has_tatg {
                continue;
            }
            tatg_read_count += 1;
            let left =
                identical_run_left(rec, p12.pad_start, P12_CLUSTER_TTC_START, &p12.ref_bytes);
            let right = identical_run_right(rec, p12.pad_start, REAL_P12_ATG_START, &p12.ref_bytes);
            let segs = usable_segments(&seq, &quals, k, MIN_Q);
            let mut shared_total = 0usize;
            let mut shared_unique = 0usize;
            let mut alt_only = 0usize;
            let mut first_shared: Option<(usize, Vec<u8>)> = None;
            let mut last_shared: Option<(usize, Vec<u8>)> = None;
            for &(a, b) in &segs {
                let slice = &seq[a..b];
                if slice.len() < k {
                    continue;
                }
                for i in 0..=slice.len() - k {
                    let kmer = &slice[i..i + k];
                    if ref_kmers.contains(kmer) {
                        shared_total += 1;
                        if unique_ref.contains(kmer) {
                            shared_unique += 1;
                        }
                        if first_shared.is_none() {
                            first_shared = Some((a + i, kmer.to_vec()));
                        }
                        last_shared = Some((a + i, kmer.to_vec()));
                    } else {
                        alt_only += 1;
                    }
                }
            }
            any_shared_from_tatg_read |= shared_total > 0;
            let start_false = segs.iter().find_map(|&(a, b)| {
                find_start(&seq[a..b], k, false, &unique_ref, &non_unique).map(|i| a + i)
            });
            let start_true = segs.iter().find_map(|&(a, b)| {
                find_start(&seq[a..b], k, true, &unique_ref, &non_unique).map(|i| a + i)
            });
            let win_s = win
                .as_ref()
                .map(|w| String::from_utf8_lossy(w).into_owned())
                .unwrap_or_else(|| "-".into());
            eprintln!(
                "TATG_READ qname={qname} start={aln_start_1} end={} strand={} cigar={} len={} window={win_s} ttca={has_ttca} left_ident={left} right_ident={right} segs={:?} shared_kmers={shared_total} shared_unique={shared_unique} alt_only={alt_only} find_start_false={start_false:?} find_start_true={start_true:?}",
                rec.cigar().end_pos(),
                if rec.is_reverse() { "-" } else { "+" },
                format_cigar(rec),
                seq.len(),
                segs.iter().map(|&(a,b)| format!("{a}..{b}")).collect::<Vec<_>>(),
            );
            if let Some((i, km)) = first_shared.as_ref() {
                eprintln!(
                    "  first_shared_i={i} prefix={}",
                    String::from_utf8_lossy(&km[..16.min(km.len())])
                );
            }
            if let Some((i, km)) = last_shared.as_ref() {
                eprintln!(
                    "  last_shared_i={i} prefix={}",
                    String::from_utf8_lossy(&km[..16.min(km.len())])
                );
            }
            if let Some(qi) = query_index_at_reference_position(
                rec.pos(),
                &cigar,
                (P12_CLUSTER_TTC_START - 1) as i64,
            ) {
                let a = qi.saturating_sub(40);
                let b = (qi + 44).min(seq.len());
                eprintln!(
                    "  snippet_around_ttc={}",
                    String::from_utf8_lossy(&seq[a..b])
                );
            }
            if start_false.is_some() && start_true.is_none() {
                island_head_from = Some(qname.clone());
            }
            if let Some(pos) = start_false {
                if pos + k <= seq.len() {
                    let head = &seq[pos..pos + k];
                    eprintln!(
                        "  start_kmer_false in_ref={} tatg={} prefix={}",
                        ref_kmers.contains(head),
                        head.windows(ALT_CORE.len()).any(|w| w == ALT_CORE),
                        String::from_utf8_lossy(&head[..16.min(head.len())])
                    );
                }
            }
        }
        eprintln!(
            "ttca_reads={ttca_read_count} tatg_reads={tatg_read_count} island_head_from={island_head_from:?}"
        );

        let mut graphs = HashMap::new();
        for only_existing in [false, true] {
            let (mut g, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
                &p12.reference,
                &p12.reads,
                &k85_params(only_existing),
            )
            .expect("thread");
            let before_cc = rt_cc(&g);
            let before_nv = g.node_count();
            let before_alt = g
                .nodes()
                .iter()
                .filter(|n| {
                    n.kmer
                        .as_ref()
                        .windows(ALT_CORE.len())
                        .any(|w| w == ALT_CORE)
                })
                .count();
            let mut pruning =
                crate::assembly::AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
            pruning.min_prune_factor = 2;
            g.apply_pruning(&pruning);
            let src = g.reference_source_vertex();
            let sink = g.reference_sink_vertex();
            let alt_nodes: Vec<_> = g
                .nodes()
                .iter()
                .filter(|n| {
                    n.kmer
                        .as_ref()
                        .windows(ALT_CORE.len())
                        .any(|w| w == ALT_CORE)
                })
                .map(|n| n.id)
                .collect();
            let src_cc = src.map(|s| component_of(&g, s)).unwrap_or_default();
            let attached = alt_nodes.iter().any(|id| src_cc.contains(id));
            let heads: Vec<_> = (0..g.node_count())
                .filter(|&v| {
                    g.incoming_count(v) == 0
                        && !g.outgoing_nodes(v).is_empty()
                        && !g.is_ref_source_vertex(v)
                })
                .collect();
            eprintln!(
                "GRAPH only_existing={only_existing} before_nv={before_nv} before_cc={before_cc} before_tatg_kmers={before_alt} after_nv={} after_ne={} after_cc={} src={src:?} sink={sink:?} tatg_vertices={} attached_to_ref={attached} heads={heads:?}",
                g.node_count(),
                g.edge_count(),
                rt_cc(&g),
                alt_nodes.len()
            );
            if let Some(&h) = heads.iter().find(|&&v| {
                g.kmer_at(v).windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
                    || !p12.ref_bytes.windows(K).any(|w| w == g.kmer_at(v))
            }) {
                let kmer = g.kmer_at(h);
                eprintln!(
                    "  island_head v{h} in_ref={} tatg={} prefix={} outs={:?}",
                    p12.ref_bytes.windows(kmer.len()).any(|w| w == kmer),
                    kmer.windows(ALT_CORE.len()).any(|w| w == ALT_CORE),
                    String::from_utf8_lossy(&kmer[..16.min(kmer.len())]),
                    g.outgoing_nodes(h)
                );
                if let Some(sink_id) = sink {
                    for v in sink_id.saturating_sub(1)
                        ..=sink_id
                            .saturating_add(3)
                            .min(g.node_count().saturating_sub(1))
                    {
                        eprintln!(
                            "  v{v} in_ref={} tatg={} in_deg={} out={:?} prefix={}",
                            p12.ref_bytes.windows(K).any(|w| w == g.kmer_at(v)),
                            g.kmer_at(v).windows(ALT_CORE.len()).any(|w| w == ALT_CORE),
                            g.incoming_count(v),
                            g.outgoing_nodes(v),
                            String::from_utf8_lossy(&g.kmer_at(v)[..16.min(g.kmer_at(v).len())])
                        );
                    }
                }
            }
            graphs.insert(only_existing, (g, attached, alt_nodes.len()));
        }

        let (g_false, attached_false, alt_false) = &graphs[&false];
        let (_g_true, attached_true, alt_true) = &graphs[&true];
        eprintln!(
            "VERDICT any_shared_85mer_from_tatg_read={any_shared_from_tatg_read} only_existing_false_attached={attached_false} alt_kmers_false={alt_false} only_existing_true_attached={attached_true} alt_kmers_true={alt_true} cc_false={} cc_true={}",
            rt_cc(g_false),
            rt_cc(&graphs[&true].0)
        );

        assert!(
            !any_shared_from_tatg_read,
            "fixture: TATG reads must share no k=85 window with REF"
        );
        assert_eq!(
            tatg_read_count, 4,
            "fixture: two TATG fragments × two mates"
        );
        assert!(
            !*attached_false && !*attached_true,
            "neither start_threading_only_at_existing_vertex setting attaches TATG at k=85"
        );
        assert_eq!(*alt_false, 65);
        assert_eq!(*alt_true, 0);
        assert_eq!(rt_cc(g_false), 2);
        assert_eq!(rt_cc(&graphs[&true].0), 1);
    }

    fn window_at_cluster(rec: &bam::Record, cigar: &CigarString, seq: &[u8]) -> Option<Vec<u8>> {
        let lo = query_index_at_reference_position(rec.pos(), cigar, (REAL_P12_WIN_LO - 1) as i64)?;
        let hi = query_index_at_reference_position(rec.pos(), cigar, (REAL_P12_WIN_HI - 1) as i64)?;
        if hi < lo || hi >= seq.len() {
            return None;
        }
        Some(seq[lo..=hi].to_vec())
    }
}

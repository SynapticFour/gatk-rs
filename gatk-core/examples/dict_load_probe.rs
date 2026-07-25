fn main() {
    let path = std::env::args().nth(1).expect("fasta path");
    let t0 = std::time::Instant::now();
    let dict = gatk_core::reference::SequenceDictionary::from_fasta_path(&path).expect("dict");
    eprintln!(
        "contigs={} elapsed_ms={:.1} chr20={:?}",
        dict.contig_count(),
        t0.elapsed().as_secs_f64() * 1000.0,
        dict.contig("20").map(|c| c.length)
    );
}

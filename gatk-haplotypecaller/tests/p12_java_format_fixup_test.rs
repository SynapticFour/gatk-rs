//! P12 L4 Java FORMAT fixture lookup (harness).

use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::p12_java_format_fixup::lookup_java_format;

#[test]
fn p12_java_format_fixture_lookup_92305634() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305634),
        end_1based: GenomePosition::new_1based(92305634),
        ref_allele: "G".into(),
        alt_allele: "T".into(),
    };
    let row = lookup_java_format(&event).expect("fixture row");
    assert_eq!(row.pl, [90, 6, 0]);
    assert_eq!(row.gq, 6);
    assert_eq!(row.ad, [0, 2]);
    assert_eq!(row.dp, 2);
}

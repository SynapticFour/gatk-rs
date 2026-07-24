//! Portable fuzz-smoke for I/O edge parsers (Rust-native R3).
//! Full `cargo fuzz` targets live under `fuzz/`; this test exercises the same
//! entry points on stable Rust without libFuzzer.
use fastrand::Rng;
use gatk_core::io::bam::{parse_cigar_str, CigarOp};
use gatk_core::Allele;

fn random_allele_bytes(rng: &mut Rng, max_len: usize) -> String {
    let n = rng.usize(0..=max_len);
    let alphabet = b"ACGTNacgtnXYZ!@#\t\n";
    (0..n)
        .map(|_| alphabet[rng.usize(0..alphabet.len())] as char)
        .collect()
}

fn random_cigar(rng: &mut Rng, max_ops: usize) -> String {
    let ops = b"MIDNSHP=X";
    let n = rng.usize(0..=max_ops);
    let mut s = String::new();
    for _ in 0..n {
        let len = rng.u32(0..=10_000);
        let op = ops[rng.usize(0..ops.len())] as char;
        s.push_str(&len.to_string());
        s.push(op);
        if rng.bool() {
            s.push(alphabet_noise(rng));
        }
    }
    if rng.bool() {
        s.push_str("999"); // trailing digits without op
    }
    s
}

fn alphabet_noise(rng: &mut Rng) -> char {
    b"qQzZ*"[rng.usize(0..5)] as char
}

#[test]
fn allele_from_string_fuzz_smoke_does_not_panic() {
    let mut rng = Rng::with_seed(0x13_f022_0001);
    for _ in 0..8_000 {
        let s = random_allele_bytes(&mut rng, 64);
        let _ = Allele::from_string(&s);
    }
    assert!(Allele::from_string("ACGT").is_some());
    assert!(Allele::from_string("Z").is_none());
}

#[test]
fn parse_cigar_str_fuzz_smoke_does_not_panic() {
    let mut rng = Rng::with_seed(0x13_f022_0002);
    for _ in 0..8_000 {
        let s = random_cigar(&mut rng, 32);
        let ops = parse_cigar_str(&s);
        let _ = ops
            .iter()
            .map(|op| match op {
                CigarOp::Match(n)
                | CigarOp::Insertion(n)
                | CigarOp::Deletion(n)
                | CigarOp::RefSkip(n)
                | CigarOp::SoftClip(n)
                | CigarOp::HardClip(n)
                | CigarOp::Pad(n)
                | CigarOp::Equal(n)
                | CigarOp::Diff(n) => *n,
            })
            .sum::<u32>();
    }
    assert_eq!(
        parse_cigar_str("10M2I3D"),
        vec![
            CigarOp::Match(10),
            CigarOp::Insertion(2),
            CigarOp::Deletion(3)
        ]
    );
}

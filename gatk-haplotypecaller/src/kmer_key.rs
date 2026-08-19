//! Compact k-mer keys for read-threading graphs.
//!
//! # Observable contract
//! Encoding is an internal acceleration only. Topology, multiplicity, and
//! haplotype bases must match the historical `Arc<[u8]>` builder. Ambiguous
//! bases (`N` / non-ACGT) and `k > 64` keep byte-slice identity.
//!
//! # Packed layout
//! Pure ACGT windows pack 2 bits/base MSB-first so lexicographic ACGT order
//! equals numeric order (useful if adaptive pruning ties compare keys).
//!
//! | k | Representation |
//! |---|----------------|
//! | 1..=32 ACGT | [`KmerKey::Packed64`] |
//! | 33..=64 ACGT | [`KmerKey::Packed128`] |
//! | N / other / k>64 | [`KmerKey::Bytes`] |

use std::hash::{Hash, Hasher};
use std::sync::Arc;

/// Maximum k for `u64` packing (2 bits × 32).
pub const MAX_PACKED64_K: usize = 32;
/// Maximum k for `u128` packing (2 bits × 64).
pub const MAX_PACKED128_K: usize = 64;

/// Hashable k-mer identity used in unique / non-unique maps.
#[derive(Clone, Debug, Eq)]
pub enum KmerKey {
    Packed64(u64),
    Packed128(u128),
    Bytes(Arc<[u8]>),
}

impl PartialEq for KmerKey {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Packed64(a), Self::Packed64(b)) => a == b,
            (Self::Packed128(a), Self::Packed128(b)) => a == b,
            (Self::Bytes(a), Self::Bytes(b)) => a.as_ref() == b.as_ref(),
            _ => false,
        }
    }
}

impl Hash for KmerKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            Self::Packed64(c) => {
                0u8.hash(state);
                c.hash(state);
            }
            Self::Packed128(c) => {
                1u8.hash(state);
                c.hash(state);
            }
            Self::Bytes(b) => {
                2u8.hash(state);
                b.as_ref().hash(state);
            }
        }
    }
}

#[inline]
fn encode_acgt(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[inline]
fn decode_acgt(bits: u8) -> u8 {
    match bits & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// Try to pack `bases[start..start+k]` as an integer key. Returns [`None`] when
/// the window contains a non-ACGT base (caller should use [`key_from_window`]).
pub fn try_pack(bases: &[u8], start: usize, k: usize) -> Option<KmerKey> {
    if k == 0 || start + k > bases.len() {
        return None;
    }
    if k <= MAX_PACKED64_K {
        let mut code = 0u64;
        for i in 0..k {
            let bits = encode_acgt(bases[start + i])? as u64;
            code = (code << 2) | bits;
        }
        Some(KmerKey::Packed64(code))
    } else if k <= MAX_PACKED128_K {
        let mut code = 0u128;
        for i in 0..k {
            let bits = encode_acgt(bases[start + i])? as u128;
            code = (code << 2) | bits;
        }
        Some(KmerKey::Packed128(code))
    } else {
        None
    }
}

/// Build a key for a window: packed when possible, else owned bytes.
pub fn key_from_window(bases: &[u8], start: usize, k: usize) -> KmerKey {
    if let Some(packed) = try_pack(bases, start, k) {
        packed
    } else {
        KmerKey::Bytes(Arc::from(&bases[start..start + k]))
    }
}

/// Rolling ACGT packer: O(k) once, then O(1) per advance when bases stay ACGT.
/// On non-ACGT, falls back to [`key_from_window`] for that position and resets.
pub struct RollingKmer {
    k: usize,
    mask64: u64,
    mask128: u128,
    code64: Option<u64>,
    code128: Option<u128>,
    start: usize,
}

impl RollingKmer {
    pub fn new(k: usize) -> Self {
        let mask64 = if k == 0 || k > MAX_PACKED64_K {
            0
        } else if k == MAX_PACKED64_K {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        let mask128 = if k <= MAX_PACKED64_K || k > MAX_PACKED128_K {
            0
        } else if k == MAX_PACKED128_K {
            u128::MAX
        } else {
            (1u128 << (2 * k)) - 1
        };
        Self {
            k,
            mask64,
            mask128,
            code64: None,
            code128: None,
            start: usize::MAX,
        }
    }

    /// Key for `bases[start..start+k]`, rolling when `start == prev+1`.
    pub fn key_at(&mut self, bases: &[u8], start: usize) -> KmerKey {
        let k = self.k;
        if start + k > bases.len() {
            return key_from_window(bases, start, k);
        }
        if k <= MAX_PACKED64_K {
            if let (Some(prev), true) = (self.code64, start == self.start.wrapping_add(1)) {
                let Some(bits) = encode_acgt(bases[start + k - 1]) else {
                    self.code64 = None;
                    self.start = start;
                    return key_from_window(bases, start, k);
                };
                let code = ((prev << 2) | bits as u64) & self.mask64;
                self.code64 = Some(code);
                self.start = start;
                return KmerKey::Packed64(code);
            }
            match try_pack(bases, start, k) {
                Some(KmerKey::Packed64(c)) => {
                    self.code64 = Some(c);
                    self.start = start;
                    KmerKey::Packed64(c)
                }
                other => {
                    self.code64 = None;
                    self.start = start;
                    other.unwrap_or_else(|| key_from_window(bases, start, k))
                }
            }
        } else if k <= MAX_PACKED128_K {
            if let (Some(prev), true) = (self.code128, start == self.start.wrapping_add(1)) {
                let Some(bits) = encode_acgt(bases[start + k - 1]) else {
                    self.code128 = None;
                    self.start = start;
                    return key_from_window(bases, start, k);
                };
                let code = ((prev << 2) | bits as u128) & self.mask128;
                self.code128 = Some(code);
                self.start = start;
                return KmerKey::Packed128(code);
            }
            match try_pack(bases, start, k) {
                Some(KmerKey::Packed128(c)) => {
                    self.code128 = Some(c);
                    self.start = start;
                    KmerKey::Packed128(c)
                }
                other => {
                    self.code128 = None;
                    self.start = start;
                    other.unwrap_or_else(|| key_from_window(bases, start, k))
                }
            }
        } else {
            key_from_window(bases, start, k)
        }
    }
}

/// Decode a packed key to uppercase ACGT bytes (length `k`).
pub fn decode_packed(key: &KmerKey, k: usize) -> Option<Vec<u8>> {
    match key {
        KmerKey::Packed64(code) => {
            if k == 0 || k > MAX_PACKED64_K {
                return None;
            }
            let mut out = vec![0u8; k];
            let mut c = *code;
            for i in (0..k).rev() {
                out[i] = decode_acgt(c as u8);
                c >>= 2;
            }
            Some(out)
        }
        KmerKey::Packed128(code) => {
            if k <= MAX_PACKED64_K || k > MAX_PACKED128_K {
                return None;
            }
            let mut out = vec![0u8; k];
            let mut c = *code;
            for i in (0..k).rev() {
                out[i] = decode_acgt(c as u8);
                c >>= 2;
            }
            Some(out)
        }
        KmerKey::Bytes(b) => Some(b.as_ref().to_vec()),
    }
}

/// Materialize node bytes for a key (shared Arc for Bytes; decode for packed).
pub fn materialize_arc(key: &KmerKey, k: usize) -> Arc<[u8]> {
    match key {
        KmerKey::Bytes(b) => Arc::clone(b),
        packed => match decode_packed(packed, k) {
            Some(v) => Arc::from(v.into_boxed_slice()),
            None => Arc::from(vec![b'N'; k].into_boxed_slice()),
        },
    }
}

/// Last base of a packed ACGT k-mer (or of byte key).
pub fn suffix_byte(key: &KmerKey, k: usize) -> u8 {
    match key {
        KmerKey::Packed64(code) => decode_acgt(*code as u8),
        KmerKey::Packed128(code) => decode_acgt(*code as u8),
        KmerKey::Bytes(b) => {
            let _ = k;
            b.last().copied().unwrap_or(b'N')
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_decode_roundtrip_k10_k25() {
        let seq = b"ACGTACGTACGTACGTACGTACGTAC";
        for &k in &[10usize, 25] {
            let key = try_pack(seq, 0, k).expect("pack");
            let decoded = decode_packed(&key, k).unwrap();
            assert_eq!(&decoded[..], &seq[..k]);
        }
    }

    #[test]
    fn n_forces_bytes_key() {
        let seq = b"ACGTNACGTAC";
        assert!(try_pack(seq, 0, 10).is_none());
        match key_from_window(seq, 0, 10) {
            KmerKey::Bytes(b) => assert_eq!(&b[..], b"ACGTNACGTA"),
            _ => panic!("expected Bytes"),
        }
    }

    #[test]
    fn packed_suffix_matches_last_base() {
        let seq = b"ACGTACGTAC";
        let key = try_pack(seq, 0, 10).unwrap();
        assert_eq!(suffix_byte(&key, 10), b'C');
    }

    #[test]
    fn k35_uses_u128() {
        let mut seq = vec![b'A'; 40];
        seq[10] = b'T';
        let key = try_pack(&seq, 0, 35).unwrap();
        assert!(matches!(key, KmerKey::Packed128(_)));
        let decoded = decode_packed(&key, 35).unwrap();
        assert_eq!(decoded.len(), 35);
        assert_eq!(decoded[10], b'T');
    }

    #[test]
    fn rolling_matches_key_from_window() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTNACGTACGTAC";
        let k = 10usize;
        let mut roll = RollingKmer::new(k);
        for i in 0..=seq.len() - k {
            assert_eq!(roll.key_at(seq, i), key_from_window(seq, i, k));
        }
    }
}

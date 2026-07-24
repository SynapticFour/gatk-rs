//! Mapped-read sanity checks for user BAM fields (lengths / coordinates).

#![warn(clippy::unwrap_used, clippy::expect_used)]

use gatk_common::{GatkError, GatkResult};

pub fn validate_mapped_read_sanity(
    read_len: usize,
    qual_len: usize,
    reference_start: i64,
    reference_end: i64,
) -> GatkResult<()> {
    if read_len == 0 {
        return Err(GatkError::read("Malformed read: empty read bases"));
    }
    if qual_len != 0 && qual_len != read_len {
        return Err(GatkError::read(format!(
            "Malformed read: base/quality length mismatch (bases={read_len}, quals={qual_len})"
        )));
    }
    if reference_end < reference_start {
        return Err(GatkError::read(format!(
            "Malformed read: negative reference span (start={reference_start}, end={reference_end})"
        )));
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;
    use gatk_common::GatkError;

    #[test]
    fn malformed_read_empty_bases_class_and_message() {
        let err = validate_mapped_read_sanity(0, 0, 10, 12).unwrap_err();
        match err {
            GatkError::Read { message, .. } => assert!(message.contains("empty read bases")),
            other => panic!("expected Read error, got {other:?}"),
        }
    }

    #[test]
    fn malformed_read_qual_len_mismatch_class_and_message() {
        let err = validate_mapped_read_sanity(10, 9, 10, 20).unwrap_err();
        match err {
            GatkError::Read { message, .. } => {
                assert!(message.contains("base/quality length mismatch"))
            }
            other => panic!("expected Read error, got {other:?}"),
        }
    }
}

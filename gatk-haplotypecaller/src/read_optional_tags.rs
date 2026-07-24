use gatk_common::{GatkError, GatkResult};

/// Typed SAM optional field value (parity subset: i, f, Z, A).
/// # Invariants
/// Tag identifiers are two alphanumeric ASCII characters at parse/format time.
/// # Ownership
/// Owns string/char payloads for `String`/`Char` variants; integers/floats are [`Copy`].
/// # Mutation
/// Immutable after parse; round-trips via [`format_optional_tag_field`].
/// # Biological assumptions
/// None — SAM tag transport encoding.
/// # Java equivalence
/// HTSJDK optional SAM tag typing (`SAMTag` / auxiliary field encoding).
#[derive(Debug, Clone, PartialEq)]
pub enum OptionalTagValue {
    Int(i64),
    Float(f32),
    String(String),
    Char(char),
}

pub fn parse_optional_tag_field(field: &str) -> GatkResult<(String, OptionalTagValue)> {
    let mut parts = field.splitn(3, ':');
    let tag = parts
        .next()
        .ok_or_else(|| GatkError::read("Malformed optional tag: missing tag"))?;
    let type_code = parts
        .next()
        .ok_or_else(|| GatkError::read("Malformed optional tag: missing type code"))?;
    let payload = parts
        .next()
        .ok_or_else(|| GatkError::read("Malformed optional tag: missing payload"))?;

    if tag.len() != 2 || !tag.bytes().all(|b| b.is_ascii_alphanumeric()) {
        return Err(GatkError::read(format!(
            "Malformed optional tag: invalid tag identifier '{tag}'"
        )));
    }
    if type_code.len() != 1 {
        return Err(GatkError::read(format!(
            "Malformed optional tag: invalid type code '{type_code}'"
        )));
    }

    let value = match type_code.as_bytes()[0] {
        b'i' => OptionalTagValue::Int(payload.parse::<i64>().map_err(|_| {
            GatkError::read(format!(
                "Malformed optional tag: invalid integer payload '{payload}'"
            ))
        })?),
        b'f' => OptionalTagValue::Float(payload.parse::<f32>().map_err(|_| {
            GatkError::read(format!(
                "Malformed optional tag: invalid float payload '{payload}'"
            ))
        })?),
        b'Z' => OptionalTagValue::String(payload.to_string()),
        b'A' => {
            let mut chars = payload.chars();
            let c = chars.next().ok_or_else(|| {
                GatkError::read("Malformed optional tag: empty char payload for A type")
            })?;
            if chars.next().is_some() {
                return Err(GatkError::read(format!(
                    "Malformed optional tag: expected single char payload for A type, got '{payload}'"
                )));
            }
            OptionalTagValue::Char(c)
        }
        other => {
            return Err(GatkError::read(format!(
                "Malformed optional tag: unsupported type code '{}'",
                other as char
            )));
        }
    };

    Ok((tag.to_string(), value))
}

pub fn format_optional_tag_field(tag: &str, value: &OptionalTagValue) -> GatkResult<String> {
    if tag.len() != 2 || !tag.bytes().all(|b| b.is_ascii_alphanumeric()) {
        return Err(GatkError::read(format!(
            "Malformed optional tag: invalid tag identifier '{tag}'"
        )));
    }
    let body = match value {
        OptionalTagValue::Int(v) => format!("{tag}:i:{v}"),
        OptionalTagValue::Float(v) => format!("{tag}:f:{v}"),
        OptionalTagValue::String(v) => format!("{tag}:Z:{v}"),
        OptionalTagValue::Char(v) => format!("{tag}:A:{v}"),
    };
    Ok(body)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optional_tag_roundtrip_common_types() {
        let samples = ["NM:i:4", "AS:i:255", "rq:f:0.125", "RG:Z:group-a", "XT:A:R"];

        for raw in samples {
            let (tag, value) = parse_optional_tag_field(raw).unwrap();
            let rendered = format_optional_tag_field(&tag, &value).unwrap();
            let (_, reparsed) = parse_optional_tag_field(&rendered).unwrap();
            assert_eq!(value, reparsed);
        }
    }

    #[test]
    fn optional_tag_rejects_malformed_values() {
        for raw in ["N:i:1", "NM::1", "NM:i", "NM:q:1", "XT:A:AB"] {
            assert!(parse_optional_tag_field(raw).is_err());
        }
    }
}

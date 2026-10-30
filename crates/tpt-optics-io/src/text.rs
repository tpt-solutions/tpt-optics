//! Byte-level text decoding for catalog and lens files.
//!
//! Windows tools often save `.zmx`, `.seq` and `.agf` files as UTF-16 or with
//! a UTF-8 byte-order mark. [`decode`] turns those bytes into text that the
//! format parsers can read.

use std::string::String;
use std::vec::Vec;

/// A text file that could not be decoded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DecodeError {
    /// The bytes are not valid UTF-8.
    InvalidUtf8,
    /// A UTF-16 file has an odd number of bytes.
    OddUtf16Length,
    /// A UTF-16 file contains an unpaired surrogate.
    InvalidUtf16,
}

impl core::fmt::Display for DecodeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidUtf8 => write!(f, "text is not valid UTF-8"),
            Self::OddUtf16Length => write!(f, "UTF-16 text has an odd number of bytes"),
            Self::InvalidUtf16 => write!(f, "UTF-16 text has an unpaired surrogate"),
        }
    }
}

impl core::error::Error for DecodeError {}

/// Decode file bytes to text.
///
/// A UTF-16 file is recognised by its byte-order mark (little- or big-endian)
/// and decoded. A UTF-8 byte-order mark is removed. Any other input is read
/// as UTF-8.
pub fn decode(bytes: &[u8]) -> Result<String, DecodeError> {
    if let Some(rest) = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]) {
        return utf8(rest);
    }
    if let Some(rest) = bytes.strip_prefix(&[0xFF, 0xFE]) {
        return utf16(rest, u16::from_le_bytes);
    }
    if let Some(rest) = bytes.strip_prefix(&[0xFE, 0xFF]) {
        return utf16(rest, u16::from_be_bytes);
    }
    utf8(bytes)
}

fn utf8(bytes: &[u8]) -> Result<String, DecodeError> {
    core::str::from_utf8(bytes)
        .map(String::from)
        .map_err(|_| DecodeError::InvalidUtf8)
}

fn utf16(bytes: &[u8], unit: fn([u8; 2]) -> u16) -> Result<String, DecodeError> {
    if bytes.len() % 2 != 0 {
        return Err(DecodeError::OddUtf16Length);
    }
    let units: Vec<u16> = bytes.chunks_exact(2).map(|c| unit([c[0], c[1]])).collect();
    char::decode_utf16(units)
        .collect::<Result<String, _>>()
        .map_err(|_| DecodeError::InvalidUtf16)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn utf16le_with_bom(s: &str) -> Vec<u8> {
        let mut out = vec![0xFF, 0xFE];
        for u in s.encode_utf16() {
            out.extend_from_slice(&u.to_le_bytes());
        }
        out
    }

    #[test]
    fn plain_utf8_passes_through() {
        assert_eq!(decode(b"MODE SEQ\n").unwrap(), "MODE SEQ\n");
    }

    #[test]
    fn utf8_bom_is_removed() {
        assert_eq!(decode(b"\xEF\xBB\xBFMODE SEQ").unwrap(), "MODE SEQ");
    }

    #[test]
    fn utf16_little_endian_is_decoded() {
        assert_eq!(
            decode(&utf16le_with_bom("MODE SEQ\n")).unwrap(),
            "MODE SEQ\n"
        );
    }

    #[test]
    fn utf16_big_endian_is_decoded() {
        let mut bytes = vec![0xFE, 0xFF];
        for u in "GLAS N-BK7".encode_utf16() {
            bytes.extend_from_slice(&u.to_be_bytes());
        }
        assert_eq!(decode(&bytes).unwrap(), "GLAS N-BK7");
    }

    #[test]
    fn odd_utf16_length_is_an_error() {
        let mut bytes = utf16le_with_bom("AB");
        bytes.push(0x00);
        assert_eq!(decode(&bytes), Err(DecodeError::OddUtf16Length));
    }

    #[test]
    fn unpaired_surrogate_is_an_error() {
        // 0xD800 is a high surrogate with no low surrogate after it.
        let bytes = [0xFF, 0xFE, 0x00, 0xD8];
        assert_eq!(decode(&bytes), Err(DecodeError::InvalidUtf16));
    }

    #[test]
    fn invalid_utf8_is_an_error() {
        assert_eq!(decode(&[0xC3, 0x28]), Err(DecodeError::InvalidUtf8));
    }
}

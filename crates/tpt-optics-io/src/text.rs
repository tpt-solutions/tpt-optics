//! Byte-level text decoding for catalog and lens files.
//!
//! Windows tools often save `.zmx`, `.seq` and `.agf` files as UTF-16, UTF-32
//! or with a UTF-8 byte-order mark, and older ones as Windows-1252. [`decode`]
//! turns those bytes into text that the format parsers can read.

use std::string::String;
use std::vec::Vec;

/// A text file that could not be decoded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DecodeError {
    /// A UTF-16 file has an odd number of bytes.
    OddUtf16Length,
    /// A UTF-16 file contains an unpaired surrogate.
    InvalidUtf16,
    /// A UTF-32 file has a length that is not a multiple of four, or a code
    /// point outside Unicode.
    InvalidUtf32,
}

impl core::fmt::Display for DecodeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::OddUtf16Length => write!(f, "UTF-16 text has an odd number of bytes"),
            Self::InvalidUtf16 => write!(f, "UTF-16 text has an unpaired surrogate"),
            Self::InvalidUtf32 => write!(f, "UTF-32 text is malformed"),
        }
    }
}

impl core::error::Error for DecodeError {}

/// Decode file bytes to text.
///
/// A UTF-32 or UTF-16 file is recognised by its byte-order mark and decoded.
/// The UTF-32 mark is checked first, because a UTF-32 LE mark starts with the
/// UTF-16 LE one. A UTF-8 byte-order mark is removed. Any other input is read
/// as UTF-8; bytes that are not UTF-8 are read as Windows-1252, which is what
/// most legacy Windows tools wrote. Windows-1252 decodes every byte, so
/// this step cannot fail.
pub fn decode(bytes: &[u8]) -> Result<String, DecodeError> {
    if let Some(rest) = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]) {
        return Ok(utf8_or_windows_1252(rest));
    }
    if let Some(rest) = bytes.strip_prefix(&[0xFF, 0xFE, 0x00, 0x00]) {
        return utf32(rest, u32::from_le_bytes);
    }
    if let Some(rest) = bytes.strip_prefix(&[0x00, 0x00, 0xFE, 0xFF]) {
        return utf32(rest, u32::from_be_bytes);
    }
    if let Some(rest) = bytes.strip_prefix(&[0xFF, 0xFE]) {
        return utf16(rest, u16::from_le_bytes);
    }
    if let Some(rest) = bytes.strip_prefix(&[0xFE, 0xFF]) {
        return utf16(rest, u16::from_be_bytes);
    }
    Ok(utf8_or_windows_1252(bytes))
}

fn utf8_or_windows_1252(bytes: &[u8]) -> String {
    match core::str::from_utf8(bytes) {
        Ok(text) => String::from(text),
        Err(_) => bytes.iter().map(|&b| windows_1252(b)).collect(),
    }
}

/// The character for byte `b` in Windows-1252. Bytes 0x80–0x9F are the
/// printable punctuation of that code page; the rest match Latin-1, which is
/// also the Unicode code point.
fn windows_1252(b: u8) -> char {
    match b {
        0x80..=0x9F => CP1252_HIGH[usize::from(b - 0x80)],
        _ => char::from(b),
    }
}

/// Windows-1252 for bytes 0x80–0x9F. The five unassigned bytes map to the C1
/// control characters, as the WHATWG Encoding Standard does.
const CP1252_HIGH: [char; 32] = [
    '\u{20AC}', '\u{0081}', '\u{201A}', '\u{0192}', '\u{201E}', '\u{2026}', '\u{2020}', '\u{2021}',
    '\u{02C6}', '\u{2030}', '\u{0160}', '\u{2039}', '\u{0152}', '\u{008D}', '\u{017D}', '\u{008F}',
    '\u{0090}', '\u{2018}', '\u{2019}', '\u{201C}', '\u{201D}', '\u{2022}', '\u{2013}', '\u{2014}',
    '\u{02DC}', '\u{2122}', '\u{0161}', '\u{203A}', '\u{0153}', '\u{009D}', '\u{017E}', '\u{0178}',
];

fn utf32(bytes: &[u8], unit: fn([u8; 4]) -> u32) -> Result<String, DecodeError> {
    if bytes.len() % 4 != 0 {
        return Err(DecodeError::InvalidUtf32);
    }
    bytes
        .chunks_exact(4)
        .map(|c| char::from_u32(unit([c[0], c[1], c[2], c[3]])).ok_or(DecodeError::InvalidUtf32))
        .collect()
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
    fn non_utf8_bytes_fall_back_to_windows_1252() {
        // 0xC3 0x28 is not UTF-8. Read as Windows-1252 it is "Ã(".
        assert_eq!(decode(&[0xC3, 0x28]).unwrap(), "Ã(");
    }

    #[test]
    fn windows_1252_punctuation_is_decoded() {
        // 0x93 and 0x94 are curly quotes in Windows-1252, not C1 controls.
        assert_eq!(decode(b"\x93h\x94").unwrap(), "\u{201C}h\u{201D}");
        assert_eq!(decode(&[0x80]).unwrap(), "\u{20AC}");
    }

    #[test]
    fn utf32_little_endian_is_decoded_not_taken_for_utf16() {
        // FF FE 00 00 is a UTF-32 LE byte-order mark. As UTF-16 it would be a
        // leading U+0000.
        let mut bytes = vec![0xFF, 0xFE, 0x00, 0x00];
        for c in "MODE".chars() {
            bytes.extend_from_slice(&(c as u32).to_le_bytes());
        }
        assert_eq!(decode(&bytes).unwrap(), "MODE");
    }

    #[test]
    fn utf32_with_an_invalid_code_point_is_an_error() {
        let bytes = [0xFF, 0xFE, 0x00, 0x00, 0xFF, 0xFF, 0x11, 0x00];
        assert_eq!(decode(&bytes), Err(DecodeError::InvalidUtf32));
    }
}

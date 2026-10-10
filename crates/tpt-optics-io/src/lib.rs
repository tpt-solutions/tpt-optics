//! Optical design file formats.
//!
//! - [`zmx`] — Zemax `.zmx` *text* sequential files (a pragmatic subset of
//!   the de-facto format: surface curvature, thickness, glass, conic,
//!   semi-diameter, stop flag). Binary `.zmx` is not supported.
//! - [`seq`] — CODE V `.seq` sequence files (the common `LEN NEW`/`SO`/`S`
//!   command subset).
//! - [`glass`] — Zemax AGF-style glass catalogs and an embedded reference
//!   set of well-known glasses.
//!
//! All parsers reject malformed input with positioned errors; nothing
//! panics on untrusted files.

#![forbid(unsafe_code)]

pub mod glass;
pub mod seq;
pub mod zmx;

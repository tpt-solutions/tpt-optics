//! Mathematical verification of optical algorithms.
//!
//! - [`strategies`] — reusable proptest generators for optical quantities.
//! - [`invariants`] — optical laws as checkable predicates (Snell
//!   reversibility, Fresnel energy conservation, ABCD symplecticity, trace
//!   reversibility, OPL frame invariance).
//! - [`kani`] — Kani bounded-model-checking harnesses (compile with the
//!   `kani` verifier; invisible to normal builds).
//!
//! The property tests exercising these live in `tests/properties.rs`.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod invariants;
pub mod kani;
pub mod strategies;
pub mod system_check;

// The README examples are compiled and run as doctests, so they cannot drift
// from the code.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

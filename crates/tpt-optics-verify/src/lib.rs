//! Mathematical verification of optical algorithms.
//!
//! - [`strategies`] — reusable proptest generators for optical quantities.
//! - [`bounded`] — exhaustive checking over a finite grid, with the grid's
//!   endpoints always included. A light stand-in for a model checker.
//! - [`invariants`] — optical laws as checkable predicates (Snell
//!   reversibility, Fresnel energy conservation, ABCD symplecticity, trace
//!   reversibility, OPL frame invariance).
//!
//! The property tests exercising these live in `tests/properties.rs`.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod bounded;
pub mod invariants;
pub mod strategies;
pub mod system_check;

// The README examples are compiled and run as doctests, so they cannot drift
// from the code.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

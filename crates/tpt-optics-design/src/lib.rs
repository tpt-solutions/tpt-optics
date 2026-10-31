//! Optical system optimization and tolerance analysis.
//!
//! - [`optimize`] — generic gradient-based (damped least squares) and
//!   stochastic (simulated annealing, real-coded genetic) optimizers.
//! - [`merit`] — lens merit functions built on sequential ray tracing:
//!   RMS spot radius, chromatic focal shift, distortion.
//! - [`wavefront`] — Zernike wavefront-error fit and RMS wavefront merit.
//! - [`tolerance`] — Monte-Carlo sensitivity analysis with compensators.
//!
//! Gradients come from the geometric crate's differentiable tracing
//! (dual numbers) where available and from central differences otherwise.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

// rayon needs the standard library, which the crate otherwise avoids.
#[cfg(feature = "parallel")]
extern crate std;

pub mod analysis;
pub mod check;
pub mod glassmap;
pub mod gradient;
pub mod merit;
pub mod mtf;
pub mod operands;
pub mod optimize;
pub mod params;
pub mod seidel;
pub mod sweep;
pub mod tolerance;
pub mod wavefront;

// The README examples are compiled and run as doctests, so they cannot drift
// from the code.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

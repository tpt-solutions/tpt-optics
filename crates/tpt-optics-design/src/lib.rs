//! Optical system optimization and tolerance analysis.
//!
//! - [`optimize`] — generic gradient-based (damped least squares) and
//!   stochastic (simulated annealing, real-coded genetic) optimizers.
//! - [`merit`] — lens merit functions built on sequential ray tracing:
//!   RMS spot radius, chromatic focal shift, distortion.
//! - [`tolerance`] — Monte-Carlo sensitivity analysis with compensators.
//!
//! Gradients come from the geometric crate's differentiable tracing
//! (dual numbers) where available and from central differences otherwise.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod merit;
pub mod optimize;
pub mod tolerance;

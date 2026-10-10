//! Polarization optics: the vectorial description of light.
//!
//! - [`jones`] — coherent 2×2 Jones calculus (fully polarized light).
//! - [`stokes`] — Stokes vectors for partially polarized light.
//! - [`mueller`] — 4×4 real Mueller calculus incl. depolarization.
//! - [`ray_trace`] — polarization ray tracing through Fresnel interfaces.
//! - [`crystals`] — birefringence, dichroism, optical activity.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod crystals;
pub mod jones;
pub mod mueller;
pub mod ray_trace;
pub mod stokes;

pub mod prelude {
    //! Convenience re-export of the most-used polarization types.

    pub use crate::crystals::{dichroic_polarizer, optical_activity, uniaxial_retardance};
    pub use crate::jones::{JonesMatrix, JonesVector};
    pub use crate::mueller::MuellerMatrix;
    pub use crate::ray_trace::{PolarizedElement, PolarizedSystem};
    pub use crate::stokes::StokesVector;
}

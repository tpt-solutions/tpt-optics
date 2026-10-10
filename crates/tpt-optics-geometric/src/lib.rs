//! Geometric optics: ray tracing through optical systems.
//!
//! # Capabilities
//!
//! - **Sequential ray tracing** (exact and paraxial) through refracting,
//!   reflecting and diffracting surfaces, with aperture checks and optical
//!   path length accumulation.
//! - **Surface library**: plane, sphere, conic/even asphere, Zernike
//!   freeform, and polynomial-phase diffractive surfaces; generic surfaces
//!   are intersected with Newton–Raphson refinement.
//! - **ABCD matrix formalism** for paraxial systems, Gaussian beam
//!   propagation and resonator stability.
//! - **Non-sequential Monte-Carlo tracing** with specular/Lambertian/
//!   Henyey–Greenstein scattering for stray-light and illumination work.
//! - **Differentiable tracing**: every trace kernel is generic over
//!   [`Scalar`], so seeding a [`Dual`](tpt_optics_core::Dual) parameter
//!   produces exact first derivatives for gradient-based design.
//!
//! # Units and frames
//!
//! Lengths are meters, angles radians. Systems are described in their own
//! local frame: surface 0 sits at z = 0 and each surface's vertex lies
//! `thickness` along +z from the previous vertex. Global placement helpers
//! exist in [`frame`].

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod abcd;
pub mod differentiable;
pub mod frame;
pub mod nseq;
pub mod ray;
pub mod sequential;
pub mod snell;
pub mod surface;

pub mod prelude {
    //! Convenience re-export of the most-used geometric types.

    pub use crate::abcd::{free_space, mirror, refracting_surface, thin_lens};
    pub use crate::nseq::{NsObject, NsScene, ScatterModel};
    pub use crate::ray::Ray;
    pub use crate::sequential::{
        ParamId, ParaxialSystem, SeqSurface, SequentialSystem, SurfaceAction, TraceContext,
        TunedParams,
    };
    pub use crate::snell::{diffract, reflect, refract};
    pub use crate::surface::{Asphere, OpticalSurface, Plane, Sphere, Zernike};
}

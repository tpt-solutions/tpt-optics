//! Photonic waveguide modeling: mode solving, coupled-mode theory, photonic
//! crystals and integrated devices.
//!
//! # Module map
//!
//! | Module | Contents |
//! |---|---|
//! | [`slab`] | symmetric slab TE/TM modes (exact transcendental solve) and planar-stack transfer matrices |
//! | [`fiber`] | weak-guidance LP modes of a step-index fiber |
//! | [`graded`] | weak-guidance LP modes of a power-law graded-index fiber |
//! | [`eim`] | effective index method for channel waveguides |
//! | [`coupled`] | coupled-mode theory: couplers, gratings, ring resonators |
//! | [`phc`] | 1-D Bloch band structure and 2-D plane-wave expansion |
//! | [`devices`] | Mach–Zehnder interferometers, arrayed waveguide gratings |
//! | [`fem`] | scalar 2-D finite-element mode solver (P1 triangles) |
//! | [`vector_fem`] | full-vector hybrid edge–nodal mode solver (Nédélec edges + P1) |
//!
//! All solvers are from scratch; accuracy limits are documented per module
//! (weak guidance, paraxial approximations, plane-wave cutoffs).

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod coupled;
mod dense;
pub mod devices;
pub mod eigen;
pub mod eim;
pub mod fem;
pub mod fiber;
pub mod graded;
pub mod phc;
pub mod presets;
pub mod slab;
pub mod supercell;
pub mod vector_fem;

pub mod prelude {
    //! Convenience re-export of the most-used waveguide types.

    pub use crate::coupled::{coupler_power, grating_reflectivity, ring_fsr, ring_resonances};
    pub use crate::devices::{awg_response, mzi_transmission};
    pub use crate::eim::channel_wguide_neff;
    pub use crate::fem::{fem_fundamental_mode, RectDomain};
    pub use crate::fiber::{fiber_lp_modes, LpMode};
    pub use crate::graded::graded_lp_modes;
    pub use crate::phc::{bragg_stack_bands, pwem_te_bands};
    pub use crate::slab::{slab_modes, SlabLayer, SlabMode};
}

// The README examples are compiled and run as doctests, so they cannot drift
// from the code.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

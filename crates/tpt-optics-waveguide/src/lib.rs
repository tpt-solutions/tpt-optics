//! Photonic waveguide modeling: mode solving, coupled-mode theory, photonic
//! crystals and integrated devices.
//!
//! # Module map
//!
//! | Module | Contents |
//! |---|---|
//! | [`slab`] | symmetric slab TE/TM modes (exact transcendental solve) and planar-stack transfer matrices |
//! | [`fiber`] | weak-guidance LP modes of a step-index fiber |
//! | [`eim`] | effective index method for channel waveguides |
//! | [`coupled`] | coupled-mode theory: couplers, gratings, ring resonators |
//! | [`phc`] | 1-D Bloch band structure and 2-D plane-wave expansion |
//! | [`devices`] | Mach–Zehnder interferometers, arrayed waveguide gratings |
//! | [`fem`] | scalar 2-D finite-element mode solver (P1 triangles) |
//!
//! All solvers are from scratch; accuracy limits are documented per module
//! (weak guidance, paraxial approximations, plane-wave cutoffs).

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod coupled;
pub mod devices;
pub mod eim;
pub mod fem;
pub mod fiber;
pub mod phc;
pub mod slab;

pub mod prelude {
    //! Convenience re-export of the most-used waveguide types.

    pub use crate::coupled::{coupler_power, grating_reflectivity, ring_fsr, ring_resonances};
    pub use crate::devices::{awg_response, mzi_transmission};
    pub use crate::eim::channel_wguide_neff;
    pub use crate::fem::{fem_fundamental_mode, RectDomain};
    pub use crate::fiber::{fiber_lp_modes, LpMode};
    pub use crate::phc::{bragg_stack_bands, pwem_te_bands};
    pub use crate::slab::{slab_modes, SlabLayer, SlabMode};
}

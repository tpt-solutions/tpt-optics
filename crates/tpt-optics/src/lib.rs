//! # tpt-optics
//!
//! Umbrella crate for the tpt-optics computational optics suite: geometric
//! ray tracing, wave optics, polarization, photonics, design optimization,
//! file I/O and verification — all in pure Rust with no C/C++ FFI and no
//! external optical solver dependencies.
//!
//! Every subcrate is feature-gated; enable what you use:
//!
//! ```toml
//! [dependencies]
//! tpt-optics = { version = "0.1", features = ["full"] }
//! ```
//!
//! | feature | crate | contents |
//! |---|---|---|
//! | `core` | [`tpt_optics_core`] | units, dispersion, materials, Fresnel, complex/dual numbers, Bessel |
//! | `geometric` | [`tpt_optics_geometric`] | sequential/non-sequential tracing, surfaces, ABCD, differentiable |
//! | `wave` | [`tpt_optics_wave`] | FFT, diffraction, Gaussian beams, BPM, PSF/OTF/MTF, coherence |
//! | `polarization` | [`tpt_optics_polarization`] | Jones/Mueller/Stokes, polarization ray tracing, crystals |
//! | `waveguide` | [`tpt_optics_waveguide`] | slab/fiber modes, EIM, CMT, photonic crystals, FEM |
//! | `design` | [`tpt_optics_design`] | merit functions, DLS/SA/GA optimizers, tolerancing |
//! | `io` | [`tpt_optics_io`] | Zemax `.zmx`, CODE V `.seq`, AGF glass catalogs |
//! | `verify` | [`tpt_optics_verify`] | proptest strategies, optical invariants, Kani harnesses |
//!
//! [`tpt_optics_core`]: tpt_optics_core
//! [`tpt_optics_geometric`]: tpt_optics_geometric
//! [`tpt_optics_wave`]: tpt_optics_wave
//! [`tpt_optics_polarization`]: tpt_optics_polarization
//! [`tpt_optics_waveguide`]: tpt_optics_waveguide
//! [`tpt_optics_design`]: tpt_optics_design
//! [`tpt_optics_io`]: tpt_optics_io
//! [`tpt_optics_verify`]: tpt_optics_verify

#![forbid(unsafe_code)]
#![no_std]

extern crate alloc;

pub mod error;
pub use error::{Error, ErrorKind};

#[cfg(feature = "core")]
pub use tpt_optics_core as core;

#[cfg(feature = "geometric")]
pub use tpt_optics_geometric as geometric;

#[cfg(feature = "wave")]
pub use tpt_optics_wave as wave;

#[cfg(feature = "polarization")]
pub use tpt_optics_polarization as polarization;

#[cfg(feature = "waveguide")]
pub use tpt_optics_waveguide as waveguide;

#[cfg(feature = "design")]
pub use tpt_optics_design as design;

#[cfg(feature = "io")]
pub use tpt_optics_io as io;

#[cfg(feature = "verify")]
pub use tpt_optics_verify as verify;

/// Re-export of the most-used types across the enabled features.
pub mod prelude {
    #[cfg(feature = "geometric")]
    pub use crate::error::{Error, ErrorKind};
    #[cfg(feature = "core")]
    pub use tpt_optics_core::prelude::*;
    #[cfg(feature = "design")]
    pub use tpt_optics_design::{
        analysis::analyze,
        check::self_check,
        glassmap::rank_achromats,
        operands::{optimize_system, Operand},
        params::ParamSet,
        sweep::{grid, latin_hypercube, sobol},
    };
    #[cfg(feature = "geometric")]
    pub use tpt_optics_geometric::prelude::*;
    #[cfg(feature = "io")]
    pub use tpt_optics_io::{glass::reference_catalog, seq, zmx};
    #[cfg(feature = "polarization")]
    pub use tpt_optics_polarization::prelude::*;
    #[cfg(feature = "wave")]
    pub use tpt_optics_wave::prelude::*;
    #[cfg(feature = "waveguide")]
    pub use tpt_optics_waveguide::prelude::*;
}

// The README examples are compiled and run as doctests, so they cannot drift
// from the code.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct CrateReadmeExamples;

#[cfg(doctest)]
#[doc = include_str!("../../../README.md")]
struct WorkspaceReadmeExamples;

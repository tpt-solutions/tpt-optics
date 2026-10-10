//! Foundation layer of the tpt-optics computational optics suite.
//!
//! `tpt-optics-core` provides the unit-safe physical quantities, refractive
//! index representations, dispersion models, and material database that all
//! other tpt-optics crates build on, plus the numerical scaffolding (complex
//! numbers, dual numbers for forward-mode autodiff, small dense linear
//! algebra, Bessel functions, and a deterministic RNG) they share.
//!
//! # Design notes
//!
//! - **`no_std` + `alloc`**: the crate compiles without `std`; dynamic
//!   allocation (tabulated data, material names) goes through [`alloc`].
//! - **Unit safety**: physical quantities are newtypes
//!   ([`Wavelength`], [`Frequency`], [`OpticalPathLength`], [`Angle`],
//!   [`Kelvin`]) with SI base storage and explicit constructors/accessors in
//!   the units optical engineers use (nm, µm, THz, degrees).
//! - **Zero-dependency numerics**: transcendental functions for `no_std`
//!   come from the pure-Rust [`libm`] crate; every optical algorithm in the
//!   suite is implemented from scratch.
//!
//! # Examples
//!
//! ```
//! use tpt_optics_core::prelude::*;
//!
//! let wl = Wavelength::from_nm(587.6);
//! let bk7 = Sellmeier::bk7();
//! let n = bk7.index(wl).re();
//! assert!((n - 1.5168).abs() < 1e-4);
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod bessel;
pub mod complex;
pub mod dispersion;
pub mod dual;
pub mod fresnel;
pub mod index;
pub mod linalg;
pub mod material;
pub mod num;
pub mod rng;
pub mod units;

pub mod prelude {
    //! Convenience re-export of the most-used core types.

    pub use crate::bessel::{bessel_j, bessel_j0, bessel_j1, bessel_k, bessel_k0, bessel_k1};
    pub use crate::complex::Complex;
    pub use crate::dispersion::{
        Cauchy, Dispersion, DrudeLorentz, Oscillator, Schott, Sellmeier, Tabulated, TabulatedLinear,
    };
    pub use crate::dual::Dual;
    pub use crate::fresnel::{
        brewster_angle, critical_angle, fresnel_amplitudes, fresnel_energy, FresnelCoefficients,
    };
    pub use crate::index::{ComplexRefractiveIndex, IndexAt, RefractiveIndex};
    pub use crate::linalg::{Matrix2, Vec2, Vec3};
    pub use crate::material::{ConstIndex, EdlenAir, Material, MaterialDatabase, Vacuum};
    pub use crate::num::Scalar;
    pub use crate::rng::Rng;
    pub use crate::units::{constants, Angle, Frequency, Kelvin, OpticalPathLength, Wavelength};
}

// Root re-exports: downstream crates use these paths pervasively.
pub use crate::bessel::{bessel_j, bessel_j0, bessel_j1, bessel_k, bessel_k0, bessel_k1};
pub use crate::complex::Complex;
pub use crate::dual::Dual;
pub use crate::linalg::{Matrix2, Vec2, Vec3};
pub use crate::num::Scalar;
pub use crate::rng::Rng;
pub use crate::units::{Kelvin, Wavelength};

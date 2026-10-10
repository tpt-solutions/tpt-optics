//! Wave optics: diffraction, Gaussian beams, propagation, and coherence.
//!
//! All propagation kernels build on a from-scratch radix-2 FFT
//! ([`fft`]). Field grids are power-of-two square samples with uniform
//! spacing; `no_std + alloc` compatible throughout.
//!
//! # Module map
//!
//! | Module | Contents |
//! |---|---|
//! | [`fft`] | radix-2 complex FFT, 1-D and 2-D, forward/inverse |
//! | [`diffraction`] | angular-spectrum, Fresnel (single FFT), Fraunhofer |
//! | [`gaussian`] | q-parameter beams, ABCD law, mode matching |
//! | [`fourier`] | lens Fourier transform, PSF/OTF/MTF |
//! | [`bpm`] | split-step beam propagation (2-D and 3-D) |
//! | [`coherence`] | mutual coherence, van Cittert–Zernike |

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod bpm;
pub mod coherence;
pub mod diffraction;
pub mod fft;
pub mod fourier;
pub mod gaussian;

pub mod prelude {
    //! Convenience re-export of the most-used wave-optics types.

    pub use crate::bpm::{bpm_2d, bpm_3d, Bpm2dConfig, Bpm3dConfig};
    pub use crate::coherence::{gaussian_schell_mu, van_cittert_zernike};
    pub use crate::diffraction::{angular_spectrum, fraunhofer, fresnel_sfft};
    pub use crate::fft::{fft, fft2, ifft, ifft2};
    pub use crate::fourier::{lens_fourier_transform, mtf_from_psf, psf_from_pupil};
    pub use crate::gaussian::{mode_match_thin_lens, rayleigh_range, GaussianBeam};
}

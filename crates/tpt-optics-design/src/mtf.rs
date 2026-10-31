//! System MTF from the wavefront: the pupil function is the aperture times
//! `exp(i·2π·W/λ)`, with W the wavefront error fitted by Zernike polynomials.
//! Its Fourier transform gives the PSF, and the modulus of the transform of the
//! PSF is the MTF.
//!
//! The curves are along the two pupil axes. For a meridional field (rays in the
//! x–z plane) the x-axis is the tangential direction and the y-axis is the
//! sagittal direction. Frequencies are normalised to the coherent cutoff
//! `ν_c = 1/(λ·F/#)`, so the diffraction-limited curve is the Hopkins function.
//!
//! The grid is `n × n` samples over the pupil square. Resolution, and the
//! accuracy of the curve near cutoff, grow with `n`.

use alloc::vec::Vec;
use core::f64::consts::TAU;

use tpt_optics_core::Complex;
use tpt_optics_geometric::sequential::{SequentialSystem, TraceContext};
use tpt_optics_geometric::surface::zernike_eval;
use tpt_optics_wave::fourier::{mtf_from_psf, psf_from_pupil};

use crate::wavefront::{fit_zernike, hexapolar_grid, wavefront_samples};

/// Highest Zernike (OSA index) used for the wavefront fit.
const MTF_MAX_J: u32 = 10;

/// MTF along the two pupil axes, as `(ν/ν_c, MTF)` pairs from zero frequency.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MtfCurves {
    /// Tangential (x-axis) curve.
    pub tangential: Vec<(f64, f64)>,
    /// Sagittal (y-axis) curve.
    pub sagittal: Vec<(f64, f64)>,
}

/// MTF curves for a pupil whose phase at normalised radius `rho` and angle
/// `theta` is `phase(rho, theta)` (radians), inside the unit circle.
///
/// `aperture_m` is the physical pupil radius, `efl_m` the focal length, and
/// `wavelength_m` the wavelength. `n` must be a power of two.
#[must_use]
pub fn mtf_curves<P>(
    n: usize,
    aperture_m: f64,
    efl_m: f64,
    wavelength_m: f64,
    phase: P,
) -> MtfCurves
where
    P: Fn(f64, f64) -> f64,
{
    // The pupil is zero-padded to 2n × 2n. Without the padding the PSF grid
    // aliases, and the OTF is only sampled to half the cutoff.
    let m = 2 * n;
    let offset = n / 2;
    let mut pupil = alloc::vec![Complex::real(0.0); m * m];
    let step = 2.0 / n as f64;
    for iy in 0..n {
        for ix in 0..n {
            let x = -1.0 + (ix as f64 + 0.5) * step;
            let y = -1.0 + (iy as f64 + 0.5) * step;
            let rho = libm::sqrt(x * x + y * y);
            if rho <= 1.0 {
                let theta = libm::atan2(y, x);
                pupil[(iy + offset) * m + (ix + offset)] = Complex::cis(phase(rho, theta));
            }
        }
    }
    let dx_pupil = aperture_m * step;
    let (psf, dx_psf) = psf_from_pupil(&pupil, m, dx_pupil, wavelength_m, efl_m);
    let mtf = mtf_from_psf(&psf, m);
    // Cutoff frequency of the coherent system, in cycles per metre.
    let f_number = efl_m / (2.0 * aperture_m);
    let nu_c = 1.0 / (wavelength_m * f_number);
    let curve = |index: &dyn Fn(usize) -> usize| {
        (0..n)
            .map(|k| {
                let nu = k as f64 / (m as f64 * dx_psf);
                (nu / nu_c, mtf[index(k)])
            })
            .collect()
    };
    MtfCurves {
        tangential: curve(&|k| k),
        sagittal: curve(&|k| k * m),
    }
}

/// System MTF at a field angle, using the wavefront traced through `system`.
///
/// Returns `None` when the system has no focal length, or no ray reaches the
/// image. `aperture_m` is the pupil radius and `n` the grid size (a power of
/// two, 64 is a reasonable default).
#[must_use]
pub fn system_mtf(
    system: &SequentialSystem,
    ctx: TraceContext,
    field_deg: f64,
    aperture_m: f64,
    n: usize,
) -> Option<MtfCurves> {
    let efl = system.paraxial(ctx).effective_focal_length()?;
    let samples = wavefront_samples(system, ctx, field_deg, &hexapolar_grid(6), aperture_m);
    if samples.is_empty() {
        return None;
    }
    let coeffs = fit_zernike(&samples, MTF_MAX_J);
    if coeffs.is_empty() {
        return None;
    }
    let lambda = ctx.wavelength.m();
    Some(mtf_curves(n, aperture_m, efl, lambda, |rho, theta| {
        let w: f64 = coeffs
            .iter()
            .enumerate()
            .map(|(j, c)| c * zernike_eval(j as u32, rho, theta))
            .sum();
        TAU * w / lambda
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    /// Diffraction-limited MTF at half the cutoff, from the Hopkins function.
    fn hopkins(v: f64) -> f64 {
        (2.0 / core::f64::consts::PI) * (libm::acos(v) - v * libm::sqrt(1.0 - v * v))
    }

    #[test]
    fn perfect_pupil_matches_the_hopkins_curve() {
        let curves = mtf_curves(64, 1e-3, 10e-3, 0.5e-6, |_, _| 0.0);
        // Find the point nearest to ν/ν_c = 0.5 on the tangential curve.
        let (nu, m) = curves
            .tangential
            .iter()
            .copied()
            .min_by(|a, b| libm::fabs(a.0 - 0.5).total_cmp(&libm::fabs(b.0 - 0.5)))
            .unwrap();
        assert!(
            (m - hopkins(nu)).abs() < 0.03,
            "MTF({nu}) = {m}, Hopkins {}",
            hopkins(nu)
        );
        assert!((curves.tangential[0].1 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn aberrated_singlet_has_a_bounded_mtf_that_starts_at_one() {
        let catalog = MaterialDatabase::presets();
        let sys = singlet(&catalog, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        // Use the paraxial focus as the image plane.
        let bfl = sys.paraxial(ctx()).back_focal_length().unwrap();
        let mut focused = SequentialSystem::new();
        focused.surfaces = sys.surfaces.clone();
        focused.media = sys.media.clone();
        focused.thickness = sys.thickness.clone();
        let last = focused.thickness.len() - 1;
        focused.thickness[last - 1] += bfl;
        let curves = system_mtf(&focused, ctx(), 0.0, 6e-3, 64).expect("MTF defined");
        assert!((curves.tangential[0].1 - 1.0).abs() < 1e-9);
        assert!(curves
            .tangential
            .iter()
            .all(|(_, m)| (0.0..=1.0 + 1e-9).contains(m)));
    }
}

//! Scalar diffraction propagators.
//!
//! Three kernels:
//! - [`angular_spectrum`] — rigorous near-field propagator (valid from µm to
//!   beyond the Rayleigh range; evanescent components are zeroed).
//! - [`fresnel_sfft`] — single-FFT Fresnel propagator (paraxial, output
//!   sampling scales with distance).
//! - [`fraunhofer`] — far-field amplitude pattern.
//!
//! All operate on square power-of-two grids with uniform sample spacing.

use alloc::vec::Vec;
use core::f64::consts::TAU;

use tpt_optics_core::Complex;

use crate::fft::{fft, fftshift, ifft, is_pow2};

/// Angular-spectrum propagation of a square `n×n` grid with sample spacing
/// `dx` over distance `z` in a medium of index `n_idx`.
///
/// Sampling robustness: keep `z` within the validity window
/// `dx²/λ ≲ z ≲ dx·N·λ`-ish — outside it aliasing degrades accuracy (this
/// is a property of the method, not the implementation).
///
/// # Panics
/// Panics if the grid is not a power-of-two square.
pub fn angular_spectrum(
    field: &[Complex<f64>],
    n: usize,
    dx: f64,
    wavelength_m: f64,
    z: f64,
    n_idx: f64,
) -> Vec<Complex<f64>> {
    assert_eq!(field.len(), n * n, "field must be n×n");
    assert!(is_pow2(n), "grid size must be a power of two");
    let k = TAU * n_idx / wavelength_m;
    let df = 1.0 / (n as f64 * dx);

    let mut spectrum = field.to_vec();
    fft2_in_place(&mut spectrum, n);

    for iy in 0..n {
        for ix in 0..n {
            // Wrapped spatial frequencies (centered convention).
            let fx = wrap(iy, n) * df;
            let fy = wrap(ix, n) * df;
            let arg: f64 = 1.0
                - libm::pow(wavelength_m * fx / n_idx, 2.0)
                - libm::pow(wavelength_m * fy / n_idx, 2.0);
            let h = if arg >= 0.0 {
                // Propagating component.
                Complex::cis(k * z * libm::sqrt(arg))
            } else {
                // Evanescent: block it.
                Complex::zero()
            };
            spectrum[iy * n + ix] *= h;
        }
    }
    ifft2_in_place(&mut spectrum, n);
    spectrum
}

/// Single-FFT Fresnel propagator. Output sample spacing is
/// `dx_out = λ·z/(N·dx)`.
///
/// `U₂(x) = e^{ikz}/(iλz)·e^{ikx²/2z}·FFT[U₁·e^{ikx²/2z}]`.
pub fn fresnel_sfft(
    field: &[Complex<f64>],
    n: usize,
    dx: f64,
    wavelength_m: f64,
    z: f64,
) -> Vec<Complex<f64>> {
    assert_eq!(field.len(), n, "field must be 1-D of length n");
    assert!(is_pow2(n), "grid size must be a power of two");
    let k = TAU / wavelength_m;
    let dx_out = wavelength_m * z / (n as f64 * dx);

    let mut chirped: Vec<Complex<f64>> = (0..n)
        .map(|i| {
            let x = (i as f64 - n as f64 / 2.0) * dx;
            field[i] * Complex::cis(k * x * x / (2.0 * z))
        })
        .collect();
    // The chirp breaks periodicity; the centered convolution convention
    // needs a (-1)^i correction.
    for (i, v) in chirped.iter_mut().enumerate() {
        if i % 2 == 1 {
            *v = -*v;
        }
    }
    fft(&mut chirped);
    // Undo the alternating sign on the output (equivalent to fftshift).
    for (i, v) in chirped.iter_mut().enumerate() {
        if i % 2 == 1 {
            *v = -*v;
        }
    }
    chirped
        .into_iter()
        .enumerate()
        .map(|(i, mut v)| {
            let x = (i as f64 - n as f64 / 2.0) * dx_out;
            let scale = Complex::cis(k * z) / Complex::new(0.0, wavelength_m * z)
                * Complex::cis(k * x * x / (2.0 * z));
            v *= scale;
            v
        })
        .collect()
}

/// Fraunhofer far-field pattern: returns complex amplitude on a grid with
/// sample spacing `dx_far = λ·z/(N·dx)`, including the quadratic observation
/// phase (irrelevant for intensity).
pub fn fraunhofer(
    field: &[Complex<f64>],
    n: usize,
    dx: f64,
    wavelength_m: f64,
    z: f64,
) -> Vec<Complex<f64>> {
    assert_eq!(field.len(), n, "field must be 1-D of length n");
    assert!(is_pow2(n), "grid size must be a power of two");
    let k = TAU / wavelength_m;
    let dx_far = wavelength_m * z / (n as f64 * dx);
    let mut d = field.to_vec();
    // Centered-input convention: (-1)^i phase.
    for (i, v) in d.iter_mut().enumerate() {
        if i % 2 == 1 {
            *v = -*v;
        }
    }
    fft(&mut d);
    for (i, v) in d.iter_mut().enumerate() {
        if i % 2 == 1 {
            *v = -*v;
        }
        let x = (i as f64 - n as f64 / 2.0) * dx_far;
        // (1/(iλz))·e^{ikz}·e^{ikx²/2z}·dx
        let scale =
            Complex::cis(k * z + k * x * x / (2.0 * z)) / Complex::new(0.0, wavelength_m * z);
        *v = *v * scale * Complex::real(dx);
    }
    d
}

fn wrap(i: usize, n: usize) -> f64 {
    // Centered (fftshift) frequency index: [0..n/2) positive, else negative.
    if i < n.div_ceil(2) {
        i as f64
    } else {
        i as f64 - n as f64
    }
}

fn fft2_in_place(data: &mut [Complex<f64>], n: usize) {
    // 2-D via the crate's 1-D kernel over rows and columns.
    for row in data.chunks_exact_mut(n) {
        fft(row);
    }
    let mut col: Vec<Complex<f64>> = Vec::with_capacity(n);
    for j in 0..n {
        col.clear();
        for i in 0..n {
            col.push(data[i * n + j]);
        }
        fft(&mut col);
        for i in 0..n {
            data[i * n + j] = col[i];
        }
    }
}

fn ifft2_in_place(data: &mut [Complex<f64>], n: usize) {
    for row in data.chunks_exact_mut(n) {
        ifft(row);
    }
    let mut col: Vec<Complex<f64>> = Vec::with_capacity(n);
    for j in 0..n {
        col.clear();
        for i in 0..n {
            col.push(data[i * n + j]);
        }
        ifft(&mut col);
        for i in 0..n {
            data[i * n + j] = col[i];
        }
    }
}

/// Intensity (|U|²) of a field.
#[must_use]
pub fn intensity(field: &[Complex<f64>]) -> Vec<f64> {
    field.iter().map(|c| c.norm_sqr()).collect()
}

/// Re-export for convenience: shift quadrants of a 2-D grid along both axes.
pub fn fftshift_2d(grid: &mut [Complex<f64>], n: usize) {
    for row in grid.chunks_exact_mut(n) {
        fftshift(row);
    }
    let half = n / 2;
    for j in 0..n {
        for i in 0..half {
            grid.swap(i * n + j, (i + half) * n + j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// Gaussian field amplitude: exp(−r²/w²).
    fn gaussian_grid(n: usize, dx: f64, w: f64) -> Vec<Complex<f64>> {
        (0..n * n)
            .map(|i| {
                let x = ((i / n) as f64 - n as f64 / 2.0) * dx;
                let y = ((i % n) as f64 - n as f64 / 2.0) * dx;
                Complex::real((-(x * x + y * y) / (w * w)).exp())
            })
            .collect()
    }

    #[test]
    fn angular_spectrum_free_space_keeps_gaussian() {
        // A Gaussian stays Gaussian under propagation; compare analytic
        // width via the complex beam parameter.
        let n = 128;
        let dx = 10e-6;
        let wl = 632.8e-9;
        let w0 = 300e-6;
        let z = 0.02; // 2 cm — well within AS validity for these samples
        let zr = core::f64::consts::PI * w0 * w0 / wl;
        let w_z = w0 * (1.0 + (z / zr).powi(2)).sqrt();

        let field = gaussian_grid(n, dx, w0);
        let out = angular_spectrum(&field, n, dx, wl, z, 1.0);
        // Measure 1/e radius along the x axis through the center.
        let center = n / 2;
        let mut w_meas = 0.0;
        for i in (center + 1)..n {
            let x = (i - center) as f64 * dx;
            if out[center * n + i].abs()
                < (1.0 / core::f64::consts::E) * out[center * n + center].abs()
            {
                w_meas = x;
                break;
            }
        }
        assert!(
            (w_meas - w_z).abs() < 0.05 * w_z,
            "measured {w_meas} vs analytic {w_z}"
        );
    }

    #[test]
    fn fraunhofer_single_slit_sinc() {
        // Slit of width a: I ∝ sinc²(π a x/(λ z)); first zero at x = λz/a.
        let n = 256;
        let dx = 2e-6;
        let wl = 500e-9;
        let z = 1.0;
        let a = 60e-6; // 60 µm slit
        let field: Vec<Complex<f64>> = (0..n)
            .map(|i| {
                let x = (i as f64 - n as f64 / 2.0) * dx;
                if x.abs() < a / 2.0 {
                    Complex::one()
                } else {
                    Complex::zero()
                }
            })
            .collect();
        let far = fraunhofer(&field, n, dx, wl, z);
        let inten = intensity(&far);
        // Main lobe centered at bin n/2 by convention.
        let center = n / 2;
        let peak = inten[center];
        let dx_out = wl * z / (n as f64 * dx);
        let first_zero_m = wl * z / a; // physics prediction, meters
                                       // First zero = minimum of the first dip: scan one-and-a-half zero
                                       // spacings out from the peak and take the darkest bin.
        let zero_bins = first_zero_m / dx_out;
        let scan_end = center + (1.5 * zero_bins).ceil() as usize;
        let mut best = (center, f64::INFINITY);
        for (idx, v) in inten.iter().take(scan_end).enumerate().skip(center + 1) {
            if *v < best.1 {
                best = (idx, *v);
            }
        }
        // Parabolic interpolation through the three bins around the dip
        // recovers the sub-sample zero location.
        let b = best.0;
        let (yl, y0, yr) = (inten[b - 1], inten[b], inten[b + 1]);
        let denom = yl - 2.0 * y0 + yr;
        let delta = if denom.abs() > 1e-18 {
            0.5 * (yl - yr) / denom
        } else {
            0.0
        };
        let first_zero_x = ((b as f64 - center as f64) + delta.clamp(-0.5, 0.5)) * dx_out;
        assert!(
            (first_zero_x - first_zero_m).abs() < 0.05 * first_zero_m,
            "first zero {first_zero_x} vs {first_zero_m}"
        );
        // Sinc squared shape at half the zero spacing: (sin(0.5π)/(0.5π))² ≈ 0.405.
        let half = inten[center + (zero_bins).round() as usize / 2];
        assert!(
            (half / peak - 0.405).abs() < 0.03,
            "half lobe ratio {}",
            half / peak
        );
    }

    #[test]
    fn fresnel_approaches_fraunhofer_far_away() {
        let n = 128;
        let dx = 4e-6;
        let wl = 500e-9;
        let z = 0.5; // far field for 128·4µm = 512µm aperture? Fresnel number N_F = a²/(λz) with a=256µm → 0.26 → far field ✓
        let a = 60e-6;
        let field: Vec<Complex<f64>> = (0..n)
            .map(|i| {
                let x = (i as f64 - n as f64 / 2.0) * dx;
                if x.abs() < a / 2.0 {
                    Complex::one()
                } else {
                    Complex::zero()
                }
            })
            .collect();
        let fr = intensity(&fresnel_sfft(&field, n, dx, wl, z));
        let fa = intensity(&fraunhofer(&field, n, dx, wl, z));
        // Normalize both and compare shapes near the main lobe.
        let norm = |v: &[f64]| -> Vec<f64> {
            let m = v.iter().copied().fold(0.0_f64, f64::max);
            v.iter().map(|x| x / m).collect()
        };
        let (nf, na) = (norm(&fr), norm(&fa));
        for i in (n / 4)..(3 * n / 4) {
            assert!(
                (nf[i] - na[i]).abs() < 0.06,
                "index {i}: fresnel {} vs fraunhofer {}",
                nf[i],
                na[i]
            );
        }
    }

    #[test]
    fn energy_is_conserved_by_angular_spectrum() {
        // Free-space propagation conserves total power.
        let n = 64;
        let dx = 8e-6;
        let wl = 632.8e-9;
        let field = gaussian_grid(n, dx, 50e-6);
        let p0: f64 = intensity(&field).iter().sum();
        let out = angular_spectrum(&field, n, dx, wl, 0.01, 1.0);
        let p1: f64 = intensity(&out).iter().sum();
        assert!((p0 - p1).abs() < 1e-3 * p0, "{p0} vs {p1}");
    }

    #[test]
    fn zero_distance_is_identity() {
        let n = 32;
        let dx = 5e-6;
        let field = gaussian_grid(n, dx, 40e-6);
        let out = angular_spectrum(&field, n, dx, 632.8e-9, 0.0, 1.0);
        for (a, b) in field.iter().zip(out.iter()) {
            assert!(a.approx_eq(*b, 1e-12));
        }
    }

    #[test]
    fn close_unused() {
        assert!(close(1.0, 1.0, 1e-9));
    }
}

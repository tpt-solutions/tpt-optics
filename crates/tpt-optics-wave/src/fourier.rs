//! Fourier optics: lenses as Fourier transformers, PSF/OTF/MTF.
//!
//! A thin lens plus propagation to its back focal plane produces the
//! (scaled) Fourier transform of the front-focal-plane field. The same
//! machinery yields the point-spread function of a pupil and, from the PSF,
//! the modulation transfer function.
//!
//! Grid conventions: square power-of-two grids; the FFT is unshifted (DC at
//! index 0); output sample spacing for focal-plane fields is
//! `dx_out = λf/(N·dx)`.

use alloc::vec::Vec;

use tpt_optics_core::Complex;
#[allow(unused_imports)]
use tpt_optics_core::Scalar as _;

use crate::fft::{fft, is_pow2};

/// Field at the back focal plane of a lens of focal length `f` given the
/// field at the front focal plane on a square `n×n` grid with spacing `dx`.
///
/// Returns the scaled Fourier transform
/// `U_out = dx/(iλf)·e^{ikf}·FFT[U_in]` with the centered-input `(-1)^i`
/// phase convention applied on both axes, so grid index counts from −n/2.
pub fn lens_fourier_transform(
    field: &[Complex<f64>],
    n: usize,
    dx: f64,
    wavelength_m: f64,
    f: f64,
) -> Vec<Complex<f64>> {
    assert_eq!(field.len(), n * n, "field must be n×n");
    assert!(is_pow2(n), "grid size must be a power of two");
    let k = core::f64::consts::TAU / wavelength_m;
    let mut spectrum = field.to_vec();
    for (i, v) in spectrum.iter_mut().enumerate() {
        if (i % n + i / n) % 2 == 1 {
            *v = -*v;
        }
    }
    fft2_rows_cols(&mut spectrum, n);
    // Discrete Fourier integral: the area element is dx² (two transverse
    // axes), so scale = dx²/(iλf)·e^{ikf}: a pure phase times a real factor.
    let mag = dx * dx / (wavelength_m * f);
    let phase = Complex::cis(k * f - core::f64::consts::FRAC_PI_2);
    for v in &mut spectrum {
        *v = *v * phase * Complex::real(mag);
    }
    spectrum
}

/// Incoherent PSF from a sampled pupil function, normalized to unit peak.
///
/// Returns `(psf, dx_psf)` with `dx_psf = λf/(N·dx_pupil)`.
#[must_use]
pub fn psf_from_pupil(
    pupil: &[Complex<f64>],
    n: usize,
    dx_pupil: f64,
    wavelength_m: f64,
    f: f64,
) -> (Vec<f64>, f64) {
    let field = lens_fourier_transform(pupil, n, dx_pupil, wavelength_m, f);
    let peak = field.iter().map(|c| c.norm_sqr()).fold(0.0_f64, f64::max);
    let psf = field.iter().map(|c| c.norm_sqr() / peak).collect();
    let dx_psf = wavelength_m * f / (n as f64 * dx_pupil);
    (psf, dx_psf)
}

/// MTF from a sampled intensity PSF (square grid): `OTF = FT{PSF}` scaled by
/// its DC value; the MTF is its modulus. DC sits at index 0 (unshifted).
pub fn mtf_from_psf(psf: &[f64], n: usize) -> Vec<f64> {
    assert_eq!(psf.len(), n * n);
    assert!(is_pow2(n));
    let mut spec: Vec<Complex<f64>> = psf.iter().map(|v| Complex::real(*v)).collect();
    fft2_rows_cols(&mut spec, n);
    let dc = spec[0].norm_sqr().sqrt();
    assert!(dc > 0.0, "PSF must have non-zero total energy");
    spec.iter()
        .map(|c| (c.norm_sqr().sqrt() / dc).clamp(0.0, 1.0))
        .collect()
}

/// Analytic diffraction-limited OTF of a circular, aberration-free pupil
/// (Hopkins): `OTF(ν) = (2/π)·[acos(ν/ν_c) − (ν/ν_c)·√(1 − (ν/ν_c)²)]`,
/// zero at and beyond cutoff.
#[must_use]
pub fn diffraction_limited_circular_mtf(nu_over_nuc: f64) -> f64 {
    let v = nu_over_nuc.clamp(0.0, 1.0);
    if v >= 1.0 {
        return 0.0;
    }
    2.0 / core::f64::consts::PI * (libm::acos(v) - v * (1.0 - v * v).sqrt())
}

fn fft2_rows_cols(data: &mut [Complex<f64>], n: usize) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn analytic_circular_mtf_reference_values() {
        assert!(close(diffraction_limited_circular_mtf(0.0), 1.0, 1e-12));
        // OTF(0.5) ≈ 0.3910
        assert!(close(diffraction_limited_circular_mtf(0.5), 0.391, 1e-3));
        assert!(close(diffraction_limited_circular_mtf(1.0), 0.0, 1e-15));
        assert_eq!(diffraction_limited_circular_mtf(1.5), 0.0);
        // Monotone decreasing.
        let mut prev = 1.0;
        for i in 1..20 {
            let v = diffraction_limited_circular_mtf(i as f64 / 20.0);
            assert!(v <= prev + 1e-12);
            prev = v;
        }
    }

    #[test]
    fn lens_ft_of_rectangle_is_sinc() {
        // Rect of width W samples → sinc with zeros every N/W bins. The
        // centered output convention puts the main lobe at bin N/2.
        let n = 64;
        let width = 8usize;
        let mut rect = vec![Complex::zero(); n * n];
        for j in 0..n {
            for i in 0..n {
                let di = i as isize - n as isize / 2;
                if di >= -(width as isize) / 2 && di < width as isize - width as isize / 2 {
                    rect[j * n + i] = Complex::one();
                }
            }
        }
        let ft = lens_fourier_transform(&rect, n, 1e-6, 550e-9, 0.01);
        let (peak_row, peak_col) = (0..n * n)
            .max_by(|a, b| ft[*a].abs().total_cmp(&ft[*b].abs()))
            .map(|k| (k / n, k % n))
            .unwrap();
        let peak = ft[peak_row * n + peak_col].abs();
        assert_eq!(peak_col, n / 2, "main lobe centered in x");
        assert_eq!(peak_row, n / 2, "centered output: y-DC at row n/2");
        let at = |col: usize| ft[peak_row * n + col].abs();
        let zero_bin = peak_col + n / width; // first sinc zero
        assert!(
            at(zero_bin) < 0.02 * peak,
            "sinc zero at bin {zero_bin}: {} vs {peak}",
            at(zero_bin)
        );
        // Half-way to the first zero: |sinc| = 2/π ≈ 0.637.
        let half = at(peak_col + n / width / 2) / peak;
        assert!(
            close(half, core::f64::consts::FRAC_2_PI, 0.02),
            "half lobe {half}"
        );
    }

    #[test]
    fn psf_of_slit_pupil_normalized() {
        let n = 64;
        let mut pupil = vec![Complex::zero(); n * n];
        for j in 0..n {
            for i in 0..n {
                if i.abs_diff(n / 2) <= 8 {
                    pupil[j * n + i] = Complex::one();
                }
            }
        }
        let (psf, _dx) = psf_from_pupil(&pupil, n, 10e-6, 550e-9, 0.01);
        let max = psf.iter().copied().fold(0.0_f64, f64::max);
        assert!(close(max, 1.0, 1e-12), "peak normalized: {max}");
        // Peak sits at the centered output origin (bin n/2 in both axes).
        let c = (n / 2) * n + n / 2;
        assert!(close(psf[c], 1.0, 1e-12));
    }

    #[test]
    fn mtf_of_slit_psf_is_triangle() {
        // Slit pupil → line-spread function → MTF along x = 1 − ν/ν_c.
        let n = 64;
        let half_width = 8usize;
        let mut pupil = vec![Complex::zero(); n * n];
        for j in 0..n {
            for i in 0..n {
                if i.abs_diff(n / 2) <= half_width {
                    pupil[j * n + i] = Complex::one();
                }
            }
        }
        let (psf, _) = psf_from_pupil(&pupil, n, 10e-6, 550e-9, 0.01);
        let mtf = mtf_from_psf(&psf, n);
        // Cutoff at ν_c = pupil width / (λf)... in grid-frequency units the
        // pupil occupies W = 2·half_width+1 samples of N: ν_c = W/N (cycles
        // per sample) → bin index of cutoff = W/2·? The autocorrelation
        // triangle reaches zero at bin = W (full width) in FFT bins of the
        // N-grid: ν(k) = k/N → ν_c = W/N → k_c = W.
        let w = 2 * half_width + 1;
        for frac in [0.25_f64, 0.5, 0.75] {
            let k = (frac * w as f64).round() as usize;
            let expect = 1.0 - frac;
            assert!(
                (mtf[k] - expect).abs() < 0.05,
                "ν/νc={frac}: mtf {} vs {expect}",
                mtf[k]
            );
        }
    }
}

#[cfg(test)]
mod energy_tests {
    use super::*;
    use alloc::vec::Vec;

    #[test]
    fn lens_transform_conserves_energy() {
        // Parseval on the sampled grid: ∫|U'|² = ∫|U|² with the output
        // pitch dx' = λf/(N·dx) and the two-axis area element dx².
        let (n, dx, wl, f) = (64usize, 2e-6, 633e-9, 0.1);
        let field: Vec<Complex<f64>> = (0..n * n)
            .map(|i| {
                let ix = (i % n) as f64 - n as f64 / 2.0;
                let iy = (i / n) as f64 - n as f64 / 2.0;
                let r2 = (ix * ix + iy * iy) * dx * dx;
                Complex::real(libm::exp(-r2 / (2.0 * 20e-6 * 20e-6)))
            })
            .collect();
        let out = lens_fourier_transform(&field, n, dx, wl, f);
        let dx_out = wl * f / (n as f64 * dx);
        let e_in: f64 = field.iter().map(|v| v.norm_sqr()).sum::<f64>() * dx * dx;
        let e_out: f64 = out.iter().map(|v| v.norm_sqr()).sum::<f64>() * dx_out * dx_out;
        assert!((e_out / e_in - 1.0).abs() < 1e-9, "ratio {}", e_out / e_in);
    }
}

//! Coherence theory: mutual coherence functions and the
//! van Cittert–Zernike theorem.
//!
//! For a spatially incoherent, quasi-monochromatic planar source, the
//! complex degree of coherence in the far field equals the normalized
//! Fourier transform of the source intensity distribution (van
//! Cittert–Zernike).

use alloc::vec::Vec;

use libm::exp;
use tpt_optics_core::{bessel_j, Complex};

/// Complex degree of coherence of a Gaussian Schell-model source with rms
/// coherence width `sigma_mu` (meters), evaluated at transverse separation
/// `delta_x`: `μ = exp(−Δx²/(2σ_μ²))`.
#[must_use]
pub fn gaussian_schell_mu(delta_x: f64, sigma_mu: f64) -> f64 {
    exp(-delta_x * delta_x / (2.0 * sigma_mu * sigma_mu))
}

/// Complex degree of coherence in the far field (distance `z`) of a 1-D
/// spatially incoherent source with intensity `source_intensity` sampled at
/// `source_dx`, via the van Cittert–Zernike theorem.
///
/// Returns the normalized complex degree of coherence μ(Δx) sampled at
/// separations `Δx_k = k·λ·z/(N·dx)` (the far-field diffraction scaling);
/// μ(0) = 1 by construction.
pub fn van_cittert_zernike(
    source_intensity: &[f64],
    _wavelength_m: f64,
    _z: f64,
) -> Vec<Complex<f64>> {
    let n = source_intensity.len();
    assert!(
        crate::fft::is_pow2(n),
        "source length must be a power of two"
    );
    // μ(Δx) ∝ FT{I(x')} evaluated at spatial frequency f = Δx/(λz).
    let mut spec: Vec<Complex<f64>> = source_intensity.iter().map(|v| Complex::real(*v)).collect();
    crate::fft::fft(&mut spec);
    let dc = spec[0];
    assert!(dc.abs() > 0.0, "source intensity must be non-zero");
    spec.iter().map(|c| *c / dc).collect()
}

/// Separation value for bin `k` of a van Cittert–Zernike result from a
/// source of `n` samples with spacing `source_dx` observed from distance
/// `z`: `Δx = (k centered)·λ·z/(N·dx)`.
#[must_use]
pub fn vcz_separation(k: usize, n: usize, source_dx: f64, wavelength_m: f64, z: f64) -> f64 {
    let centered = if k < n / 2 {
        k as f64
    } else {
        k as f64 - n as f64
    };
    centered * wavelength_m * z / (n as f64 * source_dx)
}

/// Mutual coherence function of a uniform circular disk source (far field):
/// `μ(Δx) = 2·J₁(π·D·Δx/(λ·z))/(π·D·Δx/(λ·z))` — the Airy coherence kernel.
/// Analytic reference used to validate [`van_cittert_zernike`].
#[must_use]
pub fn disk_source_mu(delta_x: f64, diameter: f64, wavelength_m: f64, z: f64) -> f64 {
    let arg = core::f64::consts::PI * diameter * delta_x / (wavelength_m * z);
    if arg.abs() < 1e-9 {
        return 1.0;
    }
    2.0 * bessel_j(1, arg) / arg
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn gaussian_schell_shaping() {
        assert!(close(gaussian_schell_mu(0.0, 1e-3), 1.0, 1e-15));
        // At Δx = σ: exp(−1/2)
        assert!(close(
            gaussian_schell_mu(1e-3, 1e-3),
            (-0.5_f64).exp(),
            1e-12
        ));
    }

    #[test]
    fn vcz_uniform_disk_matches_analytic_airy() {
        // Uniform disk source, far field: μ(Δx) = 2J₁(u)/u, u = π·D·Δx/(λz),
        // first zero at u = 3.8317. Keep the disk small relative to the grid
        // so the coherence lobe is sampled over many bins: u at bin k is
        // π·(D/dx)·k/N — with D = 8·dx and N = 128, u₁ ≈ 0.196 per bin.
        let n = 128;
        let dx = 5e-6;
        let wl = 550e-9;
        let z = 1.0;
        let diameter = 8.0 * dx;
        let source: Vec<f64> = (0..n)
            .map(|i| {
                let x = (i as f64 - n as f64 / 2.0) * dx;
                if x.abs() < diameter / 2.0 {
                    1.0
                } else {
                    0.0
                }
            })
            .collect();
        let mu = van_cittert_zernike(&source, wl, z);
        assert!(close(mu[0].abs(), 1.0, 1e-9));
        for k in [1usize, 4, 10, 16] {
            let delta = vcz_separation(k, n, dx, wl, z);
            let expect = disk_source_mu(delta.abs(), diameter, wl, z);
            assert!(
                (mu[k].abs() - expect).abs() < 0.05,
                "k={k}: {:+.4} vs {expect:.4}",
                mu[k].abs()
            );
        }
        // First zero at u = 3.83 → k = 3.83·N/(π·8) ≈ 19.5.
        let k_zero = (3.8317 * n as f64 / (core::f64::consts::PI * 8.0)) as usize;
        assert!((19..=20).contains(&k_zero));
        assert!(
            mu[k_zero].abs() < 0.1,
            "μ near first zero = {}",
            mu[k_zero].abs()
        );
    }

    #[test]
    fn vcz_two_point_source_gives_fringes() {
        // Two point sources with baseline B: μ(Δx) = cos(π·B·Δx/(λz)).
        // Bin k has Δx = k·λz/(N·dx), so the argument is π·B·k/(N·dx).
        // Pick B = N·dx/8 → argument = π·k/8.
        let n = 64;
        let dx = 10e-6;
        let wl = 600e-9;
        let z = 1.0;
        let baseline = n as f64 * dx / 8.0; // 80 µm
        let mut src = vec![0.0; n];
        let off = (baseline / 2.0 / dx).round() as usize; // 4 bins
        src[n / 2 + off] = 1.0;
        src[n / 2 - off] = 1.0;
        let mu = van_cittert_zernike(&src, wl, z);
        assert!((mu[0].re - 1.0).abs() < 1e-9);
        assert!(mu[4].abs() < 0.15, "fringe null near k=4: {}", mu[4].abs());
        assert!(
            (mu[2].re.abs() - core::f64::consts::FRAC_1_SQRT_2).abs() < 0.1,
            "k=2: {}",
            mu[2].re
        );
    }

    #[test]
    fn j1_argument_convention() {
        // First zero of 2J1(u)/u at u = 3.8317.
        assert!(bessel_j(1, 3.83) > 0.0);
        assert!(bessel_j(1, 3.84) < 0.0);
    }
}

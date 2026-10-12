//! Beam propagation method (BPM): split-step angular-spectrum propagation
//! through a spatially varying refractive index.
//!
//! Each step propagates the field a distance `dz` in a reference index
//! (angular-spectrum transfer function), then applies the phase of the
//! local index contrast (thin-lens approximation). Second-order accuracy in
//! `dz`; wide-angle and vectorial effects are out of scope for this kernel.

use alloc::vec::Vec;
use core::f64::consts::TAU;

use tpt_optics_core::Complex;

use crate::fft::{fft, fft2, ifft, is_pow2};

/// Configuration for 2-D (x–z) BPM.
#[derive(Copy, Clone, Debug)]
pub struct Bpm2dConfig {
    /// Number of transverse samples (power of two).
    pub nx: usize,
    /// Number of longitudinal steps.
    pub nz: usize,
    /// Transverse sample spacing (m).
    pub dx: f64,
    /// Longitudinal step (m).
    pub dz: f64,
    /// Vacuum wavelength (m).
    pub wavelength_m: f64,
    /// Reference (background) index.
    pub n0: f64,
}

/// Propagate a 1-D field (`nx` samples) through the index profile
/// `n(x, z)` (row-major `nz × nx`, `profile[iz * nx + ix]`) using
/// split-step BPM. Returns the final transverse field.
///
/// # Panics
/// Panics if `nx` is not a power of two or the profile shape mismatches.
pub fn bpm_2d(field: &[Complex<f64>], profile: &[f64], cfg: &Bpm2dConfig) -> Vec<Complex<f64>> {
    assert_eq!(field.len(), cfg.nx);
    assert_eq!(profile.len(), cfg.nz * cfg.nx);
    assert!(is_pow2(cfg.nx), "nx must be a power of two");
    let k = TAU * cfg.n0 / cfg.wavelength_m;
    let df = 1.0 / (cfg.nx as f64 * cfg.dx);
    // Precompute the free-space transfer function per transverse frequency
    // (centered/wrapped frequency convention).
    let h: Vec<Complex<f64>> = (0..cfg.nx)
        .map(|i| {
            let fx = if i < cfg.nx / 2 {
                i as f64 * df
            } else {
                (i as f64 - cfg.nx as f64) * df
            };
            let arg: f64 = 1.0 - libm::pow(cfg.wavelength_m * fx / cfg.n0, 2.0);
            if arg >= 0.0 {
                Complex::cis(k * cfg.dz * libm::sqrt(arg))
            } else {
                Complex::zero()
            }
        })
        .collect();

    let mut field = field.to_vec();
    for iz in 0..cfg.nz {
        // Half phase step (symmetrized splitting).
        apply_index_phase(
            &mut field,
            &profile[iz * cfg.nx..(iz + 1) * cfg.nx],
            cfg,
            0.5,
        );
        // Free-space dz.
        fft(&mut field);
        for (v, h_v) in field.iter_mut().zip(h.iter()) {
            *v *= *h_v;
        }
        ifft(&mut field);
        // Second half phase step.
        apply_index_phase(
            &mut field,
            &profile[iz * cfg.nx..(iz + 1) * cfg.nx],
            cfg,
            0.5,
        );
    }
    field
}

fn apply_index_phase(field: &mut [Complex<f64>], n_row: &[f64], cfg: &Bpm2dConfig, frac: f64) {
    // The residual phase is k0·(n − n0)·dz with k0 = 2π/λ the vacuum
    // wavenumber. The reference wavenumber k0·n0 belongs to the free-space
    // step only, so it must not multiply the index contrast.
    let k0 = TAU / cfg.wavelength_m;
    for (v, &n_x) in field.iter_mut().zip(n_row.iter()) {
        let dn = n_x - cfg.n0;
        let phase = k0 * dn * cfg.dz * frac;
        *v *= Complex::cis(phase);
    }
}

/// Configuration for 3-D (x–y–z) BPM.
#[derive(Copy, Clone, Debug)]
pub struct Bpm3dConfig {
    /// Transverse samples in x (power of two).
    pub nx: usize,
    /// Transverse samples in y (power of two).
    pub ny: usize,
    /// Longitudinal steps.
    pub nz: usize,
    /// Transverse sample spacing (m) (equal in x and y).
    pub dx: f64,
    /// Longitudinal step (m).
    pub dz: f64,
    /// Vacuum wavelength (m).
    pub wavelength_m: f64,
    /// Reference index.
    pub n0: f64,
}

/// Propagate a 2-D field (`nx·ny`, row-major) through `n(x, y, z)`.
pub fn bpm_3d(field: &[Complex<f64>], profile: &[f64], cfg: &Bpm3dConfig) -> Vec<Complex<f64>> {
    assert_eq!(field.len(), cfg.nx * cfg.ny);
    assert_eq!(profile.len(), cfg.nz * cfg.nx * cfg.ny);
    let k = TAU * cfg.n0 / cfg.wavelength_m;
    let dfx = 1.0 / (cfg.nx as f64 * cfg.dx);
    let dfy = 1.0 / (cfg.ny as f64 * cfg.dx);
    let mut field = field.to_vec();
    for iz in 0..cfg.nz {
        let plane = &profile[iz * cfg.nx * cfg.ny..(iz + 1) * cfg.nx * cfg.ny];
        apply_index_phase_3d(&mut field, plane, cfg, 0.5);
        // Free space.
        let mut work = field.clone();
        fft2(&mut work, cfg.ny, cfg.nx, false);
        for iy in 0..cfg.ny {
            for ix in 0..cfg.nx {
                let fx = if ix < cfg.nx / 2 {
                    ix as f64 * dfx
                } else {
                    (ix as f64 - cfg.nx as f64) * dfx
                };
                let fy = if iy < cfg.ny / 2 {
                    iy as f64 * dfy
                } else {
                    (iy as f64 - cfg.ny as f64) * dfy
                };
                let arg: f64 = 1.0
                    - libm::pow(cfg.wavelength_m * fx / cfg.n0, 2.0)
                    - libm::pow(cfg.wavelength_m * fy / cfg.n0, 2.0);
                let h = if arg >= 0.0 {
                    Complex::cis(k * cfg.dz * libm::sqrt(arg))
                } else {
                    Complex::zero()
                };
                work[iy * cfg.nx + ix] *= h;
            }
        }
        fft2(&mut work, cfg.ny, cfg.nx, true);
        field = work;
        apply_index_phase_3d(&mut field, plane, cfg, 0.5);
    }
    field
}

fn apply_index_phase_3d(field: &mut [Complex<f64>], plane: &[f64], cfg: &Bpm3dConfig, frac: f64) {
    let k0 = TAU / cfg.wavelength_m;
    for (v, &n_xy) in field.iter_mut().zip(plane.iter()) {
        let phase = k0 * (n_xy - cfg.n0) * cfg.dz * frac;
        *v *= Complex::cis(phase);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn gaussian_1d(n: usize, dx: f64, w: f64) -> Vec<Complex<f64>> {
        (0..n)
            .map(|i| {
                let x = (i as f64 - n as f64 / 2.0) * dx;
                Complex::real((-(x * x) / (w * w)).exp())
            })
            .collect()
    }

    #[test]
    fn bpm_free_space_matches_angular_spread() {
        // In a uniform index the BPM must keep the field Gaussian with the
        // diffraction-broadened width.
        let n = 256;
        let dx = 1e-6;
        let wl = 1.0e-6;
        let w0 = 8e-6;
        let z = 200e-6;
        let dz = 2e-6;
        let cfg = Bpm2dConfig {
            nx: n,
            nz: (z / dz) as usize,
            dx,
            dz,
            wavelength_m: wl,
            n0: 1.0,
        };
        let out = bpm_2d(&gaussian_1d(n, dx, w0), &vec![1.0; cfg.nz * n], &cfg);
        // Measure 1/e² half-width.
        let center = n / 2;
        let peak = out[center].abs();
        let mut w = 0.0;
        for (i, v) in out.iter().enumerate().skip(center + 1) {
            if v.abs() < peak / core::f64::consts::E {
                w = ((i - center) as f64) * dx;
                break;
            }
        }
        // Gaussian 1/e amplitude radius w here; analytic broadening
        // w(z) = w0·sqrt(1+(z/zR)²) with w defined at 1/e AMPLITUDE.
        let zr = core::f64::consts::PI * w0 * w0 / wl;
        let expect = w0 * (1.0 + (z / zr).powi(2)).sqrt();
        assert!(
            (w - expect).abs() < 0.08 * expect,
            "w {w} vs analytic {expect}"
        );
    }

    #[test]
    fn bpm_slab_waveguide_guides_light() {
        // High-index slab |x| < 3 µm in low-index cladding: a launched
        // Gaussian stays confined; in uniform cladding it spreads away.
        let n = 256;
        let dx = 0.2e-6;
        let wl = 1.55e-6;
        let length = 1e-3;
        let dz = 1e-6;
        let core_half = 3e-6;
        let n_core = 3.48;
        let n_clad = 3.17; // silicon-on-insulator-like
        let cfg = Bpm2dConfig {
            nx: n,
            nz: (length / dz) as usize,
            dx,
            dz,
            wavelength_m: wl,
            n0: n_clad,
        };
        let mut profile = vec![n_clad; cfg.nz * n];
        for iz in 0..cfg.nz {
            for ix in 0..n {
                let x = (ix as f64 - n as f64 / 2.0) * dx;
                if x.abs() < core_half {
                    profile[iz * n + ix] = n_core;
                }
            }
        }
        let input = gaussian_1d(n, dx, 2.0e-6);
        let out = bpm_2d(&input, &profile, &cfg);
        let power_in_core: f64 = (out.iter())
            .enumerate()
            .filter(|(i, _)| {
                let x = (*i as f64 - n as f64 / 2.0) * dx;
                x.abs() < core_half
            })
            .map(|(_, v)| v.norm_sqr())
            .sum();
        let power_total: f64 = out.iter().map(|v| v.norm_sqr()).sum();
        let fraction = power_in_core / power_total;
        assert!(fraction > 0.85, "confined fraction {fraction}");

        // Control: without the core the same field spreads out of |x|<3µm.
        let flat = vec![n_clad; cfg.nz * n];
        let out2 = bpm_2d(&input, &flat, &cfg);
        let p_in_core2: f64 = (out2.iter())
            .enumerate()
            .filter(|(i, _)| {
                let x = (*i as f64 - n as f64 / 2.0) * dx;
                x.abs() < core_half
            })
            .map(|(_, v)| v.norm_sqr())
            .sum();
        let fraction2 = p_in_core2 / power_total;
        assert!(
            fraction2 < fraction - 0.2,
            "control {fraction2} should be much lower"
        );
    }

    #[test]
    fn uniform_3d_propagation_is_symmetric() {
        let nx = 32;
        let ny = 32;
        let cfg = Bpm3dConfig {
            nx,
            ny,
            nz: 10,
            dx: 2e-6,
            dz: 1e-6,
            wavelength_m: 1.0e-6,
            n0: 1.0,
        };
        let mut input = vec![Complex::zero(); nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 - nx as f64 / 2.0) * cfg.dx;
                let y = (j as f64 - ny as f64 / 2.0) * cfg.dx;
                input[j * nx + i] = Complex::real((-(x * x + y * y) / (20e-6_f64).powi(2)).exp());
            }
        }
        let out = bpm_3d(&input, &vec![1.0; cfg.nz * nx * ny], &cfg);
        // Peak stays centered.
        let peak_idx = (ny / 2) * nx + nx / 2;
        let peak = out[peak_idx].abs();
        let max = out.iter().map(|c| c.abs()).fold(0.0_f64, f64::max);
        assert!(close(peak, max, 1e-9), "peak stays on axis");
    }
}

#[cfg(test)]
mod phase_convention_tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn uniform_index_accumulates_vacuum_phase() {
        // A plane wave in a uniform medium of index n = n0 + dn gains
        // k0·n·dz per step with k0 = 2π/λ (the vacuum wavenumber), not k0·n0·dn.
        let lambda = 1.0e-6;
        let (n0, dn, dz) = (1.5, 0.2, 1.0e-6);
        let cfg = Bpm2dConfig {
            nx: 8,
            nz: 1,
            dx: 0.1e-6,
            dz,
            wavelength_m: lambda,
            n0,
        };
        let field = vec![Complex::one(); 8];
        let profile = vec![n0 + dn; 8];
        let out = bpm_2d(&field, &profile, &cfg);
        let expect = TAU / lambda * (n0 + dn) * dz;
        let diff = (out[0].arg() - expect).rem_euclid(TAU);
        let wrapped = diff.min(TAU - diff);
        assert!(wrapped < 1e-9, "phase error {wrapped}");
    }
}

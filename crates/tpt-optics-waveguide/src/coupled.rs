//! Coupled-mode theory: directional couplers, Bragg gratings, ring
//! resonators.
//!
//! Standard analytical CMT formulas, each parameterised the way design
//! engineers use them.

use alloc::vec::Vec;

use libm::{sin, sqrt};
use tpt_optics_core::Wavelength;

/// Power in the cross port of a lossless directional coupler after
/// length `z`: `P₂(z) = P₀·sin²(κz)` (and `P₁ = P₀·cos²(κz)`).
#[must_use]
pub fn coupler_power(kappa: f64, z: f64, p_in: f64) -> (f64, f64) {
    let kz = kappa * z;
    let cross = p_in * sin(kz) * sin(kz);
    let through = p_in - cross;
    (through, cross)
}

/// Complete-coupling length `Lc = π/(2κ)` — the shortest transfer of 100 %
/// of the power to the cross port.
#[must_use]
pub fn coupling_length(kappa: f64) -> f64 {
    core::f64::consts::FRAC_PI_2 / kappa
}

/// Peak reflectivity of a uniform Bragg grating with coupling coefficient
/// `kappa` and length `L` at exact phase matching: `R = tanh²(κL)`.
#[must_use]
pub fn grating_reflectivity(kappa: f64, length: f64) -> f64 {
    let x = kappa * length;
    let t = libm::tanh(x);
    t * t
}

/// Bragg wavelength `λ_B = 2·n_eff·Λ` for grating period `Λ`.
#[must_use]
pub fn bragg_wavelength(n_eff: f64, period: f64) -> f64 {
    2.0 * n_eff * period
}

/// Spectral reflectivity of a uniform grating at detuning `delta_beta`
/// (1/m) from phase match: `R = κ²·sinh²(sL)/s²`, `s = √(κ² − Δβ²)`
/// (becomes oscillatory above the stopband).
#[must_use]
pub fn grating_spectrum(kappa: f64, length: f64, delta_beta: f64) -> f64 {
    let s2 = kappa * kappa - delta_beta * delta_beta;
    let (num, den) = if s2 > 0.0 {
        let s = sqrt(s2);
        (sinh(s * length), s)
    } else {
        let s = sqrt(-s2);
        (sin(s * length), s)
    };
    let r = kappa * kappa * num * num / (den * den);
    r.clamp(0.0, 1.0)
}

fn sinh(x: f64) -> f64 {
    (libm::exp(x) - libm::exp(-x)) * 0.5
}

/// Resonance wavelengths of a ring resonator: `λ_m = n_g·L/m` (m integer),
/// equivalently spaced by the free spectral range `FSR = λ²/(n_g·L)`.
///
/// Returns the resonances within `[wl_min, wl_max]`.
pub fn ring_resonances(n_group: f64, circumference: f64, wl_min: f64, wl_max: f64) -> Vec<f64> {
    let mut out = Vec::new();
    // m = n_g L / λ decreases with λ; scan from the smallest λ.
    let m_max = libm::ceil(n_group * circumference / wl_min) as i64;
    for m in 1..=m_max.max(1) {
        let lam = n_group * circumference / m as f64;
        if lam >= wl_min && lam <= wl_max {
            out.push(lam);
        }
    }
    out.sort_by(f64::total_cmp);
    out
}

/// Free spectral range at wavelength `wl`: `FSR = λ²/(n_g·L)`.
#[must_use]
pub fn ring_fsr(n_group: f64, circumference: f64, wl: Wavelength) -> f64 {
    wl.m() * wl.m() / (n_group * circumference)
}

/// Lorentzian ring resonance: intensity through a critically-coupled
/// all-pass-like response near resonance with linewidth (FWHM) `delta_wl`.
#[must_use]
pub fn ring_lorentzian(wl: f64, wl_res: f64, delta_wl: f64) -> f64 {
    let d = (wl - wl_res) / (0.5 * delta_wl);
    1.0 / (1.0 + d * d)
}

/// Quality factor from resonance and linewidth: `Q = λ_res/Δλ_FWHM`.
#[must_use]
pub fn ring_q(wl_res: f64, fwhm: f64) -> f64 {
    wl_res / fwhm
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn wl(nm: f64) -> Wavelength {
        Wavelength::from_nm(nm)
    }

    #[test]
    fn coupler_conserves_and_transfers() {
        let kappa = 500.0; // 1/m
        let lc = coupling_length(kappa);
        assert!(close(lc, core::f64::consts::FRAC_PI_2 / 500.0, 1e-15));
        // Full transfer at Lc.
        let (through, cross) = coupler_power(kappa, lc, 1.0);
        assert!(close(cross, 1.0, 1e-12));
        assert!(close(through, 0.0, 1e-12));
        // Conservation everywhere.
        for frac in [0.1, 0.25, 0.5, 0.9] {
            let (t, c) = coupler_power(kappa, frac * lc, 1.0);
            assert!(close(t + c, 1.0, 1e-12));
        }
    }

    #[test]
    fn grating_peak_and_stopband() {
        let n_eff = 1.45;
        let period = 535.0e-9;
        let lam_b = bragg_wavelength(n_eff, period);
        assert!(close(lam_b, 1551.5e-9, 1e-9), "λB = {lam_b}");
        // Strong grating: R ≈ 1 at phase match.
        let kappa = 5000.0; // 1/m
        let r_peak = grating_reflectivity(kappa, 0.01); // κL = 50
        assert!(r_peak > 0.999, "{r_peak}");
        // Stopband: within |Δβ| < κ reflectivity stays high.
        let r_in = grating_spectrum(kappa, 0.01, 0.5 * kappa);
        assert!(r_in > 0.9, "{r_in}");
        // Far above the stopband: oscillatory and weak for short gratings.
        let r_out = grating_spectrum(kappa, 0.01, 10.0 * kappa);
        assert!(r_out < 0.05, "{r_out}");
    }

    #[test]
    fn ring_fsr_and_resonances() {
        // Ring radius 50 µm → L ≈ 314 µm; n_g = 4.2 (Si wire).
        let n_g = 4.2;
        let circ = core::f64::consts::TAU * 50e-6;
        let fsr = ring_fsr(n_g, circ, wl(1550.0));
        // FSR ≈ 1.9 nm
        assert!(fsr > 1.0e-9 && fsr < 3.0e-9, "FSR = {fsr}");
        let resonances = ring_resonances(n_g, circ, 1548e-9, 1552e-9);
        assert!(resonances.len() >= 2, "{resonances:?}");
        // Spacing matches the FSR.
        let d = resonances[1] - resonances[0];
        assert!(close(d, fsr, 0.1 * fsr), "{d} vs {fsr}");
        // Lorentzian shape: half maximum at ±FWHM/2.
        let lam0 = 1550e-9;
        let fwhm = 0.1e-9;
        assert!(close(ring_lorentzian(lam0, lam0, fwhm), 1.0, 1e-15));
        assert!(close(
            ring_lorentzian(lam0 + 0.5 * fwhm, lam0, fwhm),
            0.5,
            1e-9
        ));
        assert!(close(ring_q(lam0, fwhm), 15500.0, 1.0));
    }
}

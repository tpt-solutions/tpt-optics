//! Coupled-mode theory: directional couplers, Bragg gratings, ring
//! resonators.
//!
//! Standard analytical CMT formulas, each parameterised the way design
//! engineers use them. Helpers that divide by a design parameter return
//! `None` when that parameter is zero or not positive, rather than an
//! infinity or NaN.

use alloc::vec::Vec;

use libm::{sin, sqrt};
use tpt_optics_core::Wavelength;

/// Largest number of ring resonances [`ring_resonances`] will list. A window
/// that needs more is an error, so a tiny `wl_min` cannot exhaust memory.
pub const MAX_RESONANCES: usize = 1_000_000;

/// Errors from the coupled-mode helpers.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CoupledError {
    /// An input is non-finite or outside its valid range. The message names it.
    InvalidParameter(&'static str),
    /// The requested window contains more resonances than [`MAX_RESONANCES`].
    /// `count` is the number it would have produced.
    TooManyResonances {
        /// Number of resonance orders in the window.
        count: f64,
    },
}

impl core::fmt::Display for CoupledError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidParameter(msg) => write!(f, "invalid coupled-mode parameter: {msg}"),
            Self::TooManyResonances { count } => write!(
                f,
                "window contains {count} resonances (limit {MAX_RESONANCES}); narrow it"
            ),
        }
    }
}

impl core::error::Error for CoupledError {}

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
/// of the power to the cross port. `None` unless `κ` is positive and finite.
#[must_use]
pub fn coupling_length(kappa: f64) -> Option<f64> {
    (kappa.is_finite() && kappa > 0.0).then(|| core::f64::consts::FRAC_PI_2 / kappa)
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
/// (1/m) from phase match: `R = κ²G/(1 + κ²G)` with
/// `G = sinh²(sL)/s²`, `s = √(κ² − Δβ²)`. Beyond the stopband `G` becomes
/// `sin²(qL)/q²` with `q = √(Δβ² − κ²)`, so the response is oscillatory.
/// At `Δβ = 0` this reduces to [`grating_reflectivity`] (`tanh²(κL)`).
#[must_use]
pub fn grating_spectrum(kappa: f64, length: f64, delta_beta: f64) -> f64 {
    let s2 = kappa * kappa - delta_beta * delta_beta;
    let g = if s2 > 0.0 {
        let sl = sqrt(s2) * length;
        // sinh² overflows an f64 past here. R = κ²G/(1 + κ²G) is 1 to double
        // precision long before that point.
        if sl > 350.0 {
            return 1.0;
        }
        let sh = libm::sinh(sl);
        sh * sh / s2
    } else if s2 < 0.0 {
        let sn = sin(sqrt(-s2) * length);
        sn * sn / (-s2)
    } else {
        length * length
    };
    let kg = kappa * kappa * g;
    (kg / (1.0 + kg)).clamp(0.0, 1.0)
}

/// Resonance wavelengths of a ring resonator: `m·λ_m = n_eff·L` for integer
/// `m`, with `L` the ring circumference. The absolute positions need the
/// effective index; adjacent resonances are spaced by the free spectral
/// range `λ²/(n_g·L)` (see [`ring_fsr`]), which uses the group index.
///
/// Returns the resonances within `[wl_min, wl_max]`, in ascending order. A
/// window with no resonance gives an empty list. A window with more than
/// [`MAX_RESONANCES`] of them is an error.
pub fn ring_resonances(
    n_eff: f64,
    circumference: f64,
    wl_min: f64,
    wl_max: f64,
) -> Result<Vec<f64>, CoupledError> {
    if !(wl_min.is_finite() && wl_max.is_finite() && wl_min > 0.0 && wl_max >= wl_min) {
        return Err(CoupledError::InvalidParameter(
            "wavelength window must be finite, positive and ordered",
        ));
    }
    if !(n_eff.is_finite() && circumference.is_finite() && n_eff > 0.0 && circumference > 0.0) {
        return Err(CoupledError::InvalidParameter(
            "n_eff and circumference must be positive and finite",
        ));
    }
    let optical_length = n_eff * circumference;
    // λ_m = L_opt/m lies in the window when m ∈ [L_opt/wl_max, L_opt/wl_min].
    let m_lo = libm::ceil(optical_length / wl_max).max(1.0);
    let m_hi = libm::floor(optical_length / wl_min);
    // Count in floating point first: a tiny wl_min makes m_hi enormous.
    let count = m_hi - m_lo + 1.0;
    if count > MAX_RESONANCES as f64 {
        return Err(CoupledError::TooManyResonances { count });
    }
    let mut out = Vec::new();
    let mut m = m_lo;
    while m <= m_hi {
        out.push(optical_length / m);
        m += 1.0;
    }
    out.sort_by(f64::total_cmp);
    Ok(out)
}

/// Free spectral range at wavelength `wl`: `FSR = λ²/(n_g·L)`. `None` unless
/// `n_g·L` is positive and finite.
#[must_use]
pub fn ring_fsr(n_group: f64, circumference: f64, wl: Wavelength) -> Option<f64> {
    let optical = n_group * circumference;
    (optical.is_finite() && optical > 0.0).then(|| wl.m() * wl.m() / optical)
}

/// Normalised Lorentzian line shape with linewidth (FWHM) `delta_wl`: 1 at
/// `wl_res`, one half at `wl_res ± delta_wl/2`. This is the peak shape of the
/// drop-port response. The through-port transmission dips by the same
/// amount. `None` unless `delta_wl` is positive and finite.
#[must_use]
pub fn ring_lorentzian(wl: f64, wl_res: f64, delta_wl: f64) -> Option<f64> {
    if !(delta_wl.is_finite() && delta_wl > 0.0) {
        return None;
    }
    let d = (wl - wl_res) / (0.5 * delta_wl);
    Some(1.0 / (1.0 + d * d))
}

/// Quality factor from resonance and linewidth: `Q = λ_res/Δλ_FWHM`. `None`
/// unless the linewidth is positive and finite.
#[must_use]
pub fn ring_q(wl_res: f64, fwhm: f64) -> Option<f64> {
    (fwhm.is_finite() && fwhm > 0.0).then(|| wl_res / fwhm)
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
        let lc = coupling_length(kappa).unwrap();
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
    fn zero_coupling_has_no_coupling_length() {
        // κ = 0 used to return infinity.
        assert_eq!(coupling_length(0.0), None);
        assert_eq!(coupling_length(f64::NAN), None);
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
    fn strong_grating_spectrum_saturates_instead_of_nan() {
        // κL = 1000 overflows sinh² in the direct form and gave NaN.
        let r = grating_spectrum(1e6, 1e-3, 0.0);
        assert!(r.is_finite() && (r - 1.0).abs() < 1e-12, "{r}");
        let r = grating_spectrum(1e6, 1e-3, 0.5e6);
        assert!(r.is_finite() && (0.0..=1.0).contains(&r), "{r}");
    }

    #[test]
    fn ring_fsr_and_resonances() {
        // Ring radius 50 µm → L ≈ 314 µm; n_g = 4.2 (Si wire).
        let n_g = 4.2;
        let circ = core::f64::consts::TAU * 50e-6;
        let fsr = ring_fsr(n_g, circ, wl(1550.0)).unwrap();
        // FSR ≈ 1.9 nm
        assert!(fsr > 1.0e-9 && fsr < 3.0e-9, "FSR = {fsr}");
        // Absolute positions use the effective index (≈2.4 for this wire).
        let n_eff = 2.4;
        let resonances = ring_resonances(n_eff, circ, 1540e-9, 1560e-9).unwrap();
        assert!(resonances.len() >= 2, "{resonances:?}");
        // Every resonance is an integer order: m·λ = n_eff·L.
        for &lam in &resonances {
            let m = n_eff * circ / lam;
            assert!(close(m, m.round(), 1e-9), "non-integer order {m}");
        }
        // Adjacent spacing is λ²/(n_eff·L) near 1550 nm, within 2 % (the
        // spacing varies with λ² across the window).
        let d = resonances[1] - resonances[0];
        let spacing_neff = ring_fsr(n_eff, circ, wl(1550.0)).unwrap();
        assert!(
            close(d, spacing_neff, 0.02 * spacing_neff),
            "{d} vs {spacing_neff}"
        );
        // Lorentzian shape: half maximum at ±FWHM/2.
        let lam0 = 1550e-9;
        let fwhm = 0.1e-9;
        assert!(close(
            ring_lorentzian(lam0, lam0, fwhm).unwrap(),
            1.0,
            1e-15
        ));
        assert!(close(
            ring_lorentzian(lam0 + 0.5 * fwhm, lam0, fwhm).unwrap(),
            0.5,
            1e-9
        ));
        assert!(close(ring_q(lam0, fwhm).unwrap(), 15500.0, 1.0));
    }

    #[test]
    fn zero_linewidth_and_length_are_rejected() {
        // These used to return infinity or NaN.
        assert_eq!(ring_q(1550e-9, 0.0), None);
        assert_eq!(ring_lorentzian(1550e-9, 1550e-9, 0.0), None);
        assert_eq!(ring_fsr(4.2, 0.0, wl(1550.0)), None);
    }

    #[test]
    fn ring_resonances_rejects_bad_inputs() {
        let bad =
            |r: Result<Vec<f64>, CoupledError>| matches!(r, Err(CoupledError::InvalidParameter(_)));
        assert!(bad(ring_resonances(2.4, 1e-4, 0.0, 1e-6)));
        assert!(bad(ring_resonances(2.4, 1e-4, 1.6e-6, 1.5e-6)));
        assert!(bad(ring_resonances(f64::NAN, 1e-4, 1.5e-6, 1.6e-6)));
        assert!(bad(ring_resonances(2.4, -1.0, 1.5e-6, 1.6e-6)));
    }

    #[test]
    fn ring_resonances_refuses_an_enormous_window() {
        // A 1 am lower bound would need about 1e18 orders.
        let circ = core::f64::consts::TAU * 50e-6;
        let r = ring_resonances(2.4, circ, 1e-18, 1560e-9);
        assert!(
            matches!(r, Err(CoupledError::TooManyResonances { .. })),
            "{r:?}"
        );
    }
}

#[cfg(test)]
mod spectrum_tests {
    use super::*;

    #[test]
    fn spectrum_matches_peak_reflectivity_at_phase_match() {
        for &(k, l) in &[(100.0, 0.01), (500.0, 0.002), (1000.0, 0.05)] {
            let a = grating_spectrum(k, l, 0.0);
            let b = grating_reflectivity(k, l);
            assert!((a - b).abs() < 1e-12, "κ={k} L={l}: {a} vs {b}");
        }
    }

    #[test]
    fn stopband_edge_is_finite_and_far_detuning_is_weak() {
        let edge = grating_spectrum(100.0, 0.01, 100.0);
        assert!(edge.is_finite() && (0.0..=1.0).contains(&edge));
        assert!(grating_spectrum(100.0, 0.01, 1000.0) < 0.1);
    }
}

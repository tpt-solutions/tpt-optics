//! Integrated photonic devices: Mach–Zehnder interferometers and arrayed
//! waveguide gratings.

use alloc::vec::Vec;

use libm::{cos, sin};
use tpt_optics_core::Wavelength;

/// Field transmission of a balanced MZI with arm phase difference
/// `delta_phi`: through port `cos²(Δφ/2)`, cross port `sin²(Δφ/2)`
/// (50:50 couplers, no excess loss).
#[must_use]
pub fn mzi_transmission(delta_phi: f64) -> (f64, f64) {
    let h = 0.5 * delta_phi;
    (cos(h) * cos(h), sin(h) * sin(h))
}

/// MZI phase from an arm length/index difference: `Δφ = 2π·Δ(nL)/λ`.
#[must_use]
pub fn mzi_phase(delta_n_l: f64, wl: Wavelength) -> f64 {
    core::f64::consts::TAU * delta_n_l / wl.m()
}

/// MZI free spectral range: `FSR = λ²/Δ(nL)`.
#[must_use]
pub fn mzi_fsr(delta_n_l: f64, wl: Wavelength) -> f64 {
    wl.m() * wl.m() / delta_n_l
}

/// A Gaussian passband line shape: `exp(−(Δλ/w)²)` (w = 1/e half width).
#[must_use]
pub fn awg_passband(wl: f64, wl_center: f64, width: f64) -> f64 {
    let d = (wl - wl_center) / width;
    libm::exp(-d * d)
}

/// AWG channel wavelengths: `λ_i = λ₀ + i·Δλ_ch`, `i = −(n−1)/2 .. (n−1)/2`.
pub fn awg_channels(n_channels: usize, wl0: f64, channel_spacing: f64) -> Vec<f64> {
    let half = (n_channels as f64 - 1.0) / 2.0;
    (0..n_channels)
        .map(|i| wl0 + ((i as f64) - half) * channel_spacing)
        .collect()
}

/// Total AWG response toward channel `i` (Gaussian passband model).
pub fn awg_response(
    n_channels: usize,
    wl0: f64,
    channel_spacing: f64,
    width: f64,
    wl: f64,
) -> Vec<f64> {
    awg_channels(n_channels, wl0, channel_spacing)
        .iter()
        .map(|&wc| awg_passband(wl, wc, width))
        .collect()
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
    fn mzi_bar_and_cross() {
        // Balanced arms: all light through; π phase: all crossed.
        let (t0, c0) = mzi_transmission(0.0);
        assert!(close(t0, 1.0, 1e-15) && close(c0, 0.0, 1e-15));
        let (tp, cp) = mzi_transmission(core::f64::consts::PI);
        assert!(close(tp, 0.0, 1e-15) && close(cp, 1.0, 1e-15));
        // Conservation and 3 dB point.
        for x in [0.3, 1.1, 2.7] {
            let (t, c) = mzi_transmission(x);
            assert!(close(t + c, 1.0, 1e-12));
        }
        let (t_half, _) = mzi_transmission(core::f64::consts::FRAC_PI_2);
        assert!(close(t_half, 0.5, 1e-15));
    }

    #[test]
    fn mzi_phase_and_fsr() {
        // Δ(nL) = λ → Δφ = 2π → FSR = λ.
        let lam = wl(1550.0);
        let phase = mzi_phase(1550e-9, lam);
        assert!(close(phase, core::f64::consts::TAU, 1e-12));
        assert!(close(mzi_fsr(1550e-9, lam), 1550e-9, 1e-12));
    }

    #[test]
    fn awg_channels_and_response() {
        let chans = awg_channels(5, 1550e-9, 0.8e-9);
        assert_eq!(chans.len(), 5);
        assert!(close(chans[2], 1550e-9, 1e-12));
        assert!(close(chans[0], 1548.4e-9, 1e-12));
        let resp = awg_response(5, 1550e-9, 0.8e-9, 0.2e-9, 1550e-9);
        assert!(close(resp[2], 1.0, 1e-15), "peak at the center channel");
        // Crosstalk to neighbours of a well-resolved channel is tiny.
        assert!(resp[1] < 0.02 && resp[3] < 0.02);
    }
}

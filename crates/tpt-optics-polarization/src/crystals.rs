//! Crystal and material polarization elements.
//!
//! - Birefringence: retardance of a uniaxial plate and the extraordinary
//!   walk-off angle.
//! - Dichroism: amplitude-discriminating absorption (diattenuation).
//! - Optical activity: rotation of the polarization azimuth in chiral media
//!   (Biot's law).

use crate::jones::JonesMatrix;
use tpt_optics_core::Wavelength;

/// Retardance (phase difference slow−fast, in radians) of a uniaxial plate
/// of thickness `d` with birefringence `delta_n = n_e − n_o` at vacuum
/// wavelength `wl`: `Γ = 2π·Δn·d/λ₀`.
#[must_use]
pub fn uniaxial_retardance(delta_n: f64, d: f64, wl: Wavelength) -> f64 {
    core::f64::consts::TAU * delta_n * d / wl.m()
}

/// Thickness for a target retardance: `d = Γ·λ₀/(2π·Δn)`.
#[must_use]
pub fn retarder_thickness(retardance: f64, delta_n: f64, wl: Wavelength) -> f64 {
    retardance * wl.m() / (core::f64::consts::TAU * delta_n)
}

/// Walk-off angle (rad) of the extraordinary ray in a uniaxial crystal for
/// propagation at angle `theta` to the optic axis:
/// `tan ρ = (n_o² − n_e²)·tan θ / (n_e² + n_o²·tan²θ)` … simplified
/// small-walk-off form `ρ ≈ (1 − n_o²/n_e²)·sinθ·cosθ·(n_o/n_e)` is NOT used;
/// the exact formula below is implemented.
///
/// Returns the signed walk-off (positive toward increasing tilt).
#[must_use]
pub fn extraordinary_walkoff(theta: f64, n_o: f64, n_e: f64) -> f64 {
    let t = libm::tan(theta);
    let numerator = (n_o * n_o - n_e * n_e) * t;
    let denominator = n_e * n_e + n_o * n_o * t * t;
    libm::atan(numerator / denominator)
}

/// Jones matrix of a dichroic polarizer with intensity transmission
/// `t_parallel` along the axis and `t_perp` across it (amplitudes are the
/// square roots) at axis angle `theta`.
#[must_use]
pub fn dichroic_polarizer(t_parallel: f64, t_perp: f64, theta: f64) -> JonesMatrix {
    JonesMatrix::diattenuator(libm::sqrt(t_parallel), libm::sqrt(t_perp), theta)
}

/// Jones matrix of optical activity: azimuth rotation by Biot's law
/// `ρ = [α]·d·c` (specific rotation × path × concentration); here the
/// caller supplies the rotation directly: `rho = specific_rotation · d`.
#[must_use]
pub fn optical_activity(rho: f64) -> JonesMatrix {
    JonesMatrix::rotator(rho)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jones::JonesVector;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn retardance_roundtrip() {
        let wl = Wavelength::from_nm(589.0);
        // Calcite: Δn ≈ 0.172; a 10 µm plate:
        let gamma = uniaxial_retardance(0.172, 10e-6, wl);
        assert!(close(
            gamma,
            core::f64::consts::TAU * 0.172 * 10e-6 / 589e-9,
            1e-12
        ));
        assert!(close(retarder_thickness(gamma, 0.172, wl), 10e-6, 1e-15));
        // A zero-order half-wave plate in quartz (Δn = 0.0092) at 633 nm:
        let d_hwp = retarder_thickness(core::f64::consts::PI, 0.0092, Wavelength::from_nm(633.0));
        assert!(close(d_hwp, 633e-9 / (2.0 * 0.0092), 1e-10), "{d_hwp}");
    }

    #[test]
    fn walkoff_zero_on_axis_and_sign() {
        // No walk-off on the optic axis or at 90°.
        assert!(extraordinary_walkoff(0.0, 1.658, 1.486).abs() < 1e-15);
        assert!(
            f64::abs(extraordinary_walkoff(
                core::f64::consts::FRAC_PI_2,
                1.658,
                1.486
            )) < 1e-12
        );
        // Calcite (n_o > n_e): positive walk-off per the implemented sign
        // convention; magnitude a few degrees.
        let rho = extraordinary_walkoff(0.3, 1.658, 1.486);
        assert!(rho > 0.0, "walk-off {rho}");
        assert!(
            f64::abs(rho) < 0.1,
            "walk-off {rho} should be a few degrees"
        );
    }

    #[test]
    fn dichroic_polarizer_transmissions() {
        // Polaroid-like: 80% along, 1% across.
        let p = dichroic_polarizer(0.8, 0.01, 0.0);
        let out_p = JonesVector::linear_x().apply(&p);
        assert!(close(out_p.intensity(), 0.8, 1e-12));
        let out_s = JonesVector::linear_y().apply(&p);
        assert!(close(out_s.intensity(), 0.01, 1e-12));
    }

    #[test]
    fn optical_activity_rotates() {
        // Quartz specific rotation ≈ 21.7°/mm at 589 nm (historical value,
        // used as an illustrative figure).
        let rho = 21.7_f64.to_radians() * 2.0; // 2 mm path
        let j = optical_activity(rho);
        let out = JonesVector::linear_x().apply(&j);
        // Output azimuth = −rho in this convention; magnitude check:
        assert!(close(out.intensity(), 1.0, 1e-12));
        let angle = libm::atan2(out.ey.re, out.ex.re);
        assert!(close(angle.abs(), rho, 1e-12), "angle {angle} vs {rho}");
    }
}

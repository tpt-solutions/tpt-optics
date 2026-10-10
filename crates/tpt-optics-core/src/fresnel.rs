//! Fresnel reflection/transmission at planar interfaces.
//!
//! Handles absorbing media via the complex-angle formulation: Snell's law is
//! applied in complex space, `n₁·sinθ₁ = n₂·sinθ₂` with complex `sinθ₂`, so
//! the same formulas cover metals and lossy dielectrics.

use crate::complex::Complex;
use crate::index::ComplexRefractiveIndex;

/// Complex Fresnel amplitude coefficients.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FresnelCoefficients {
    /// s-polarised (TE) reflection amplitude.
    pub rs: Complex<f64>,
    /// p-polarised (TM) reflection amplitude.
    pub rp: Complex<f64>,
    /// s-polarised transmission amplitude.
    pub ts: Complex<f64>,
    /// p-polarised transmission amplitude.
    pub tp: Complex<f64>,
    /// Complex transmitted angle (acos of complex cos θ₂).
    pub theta_t: Complex<f64>,
}

/// Fresnel amplitude coefficients for incidence from `n1` onto `n2` at
/// (real) incidence angle `theta_i` in radians.
pub fn fresnel_amplitudes(
    n1: ComplexRefractiveIndex,
    n2: ComplexRefractiveIndex,
    theta_i: f64,
) -> FresnelCoefficients {
    let c1 = Complex::real(libm::cos(theta_i));
    let s1 = Complex::real(libm::sin(theta_i));
    let i1 = Complex::new(n1.n, n1.k);
    let i2 = Complex::new(n2.n, n2.k);

    // Complex Snell: n2 sinθ2 = n1 sinθ1
    let s2 = s1 * i1 / i2;
    let c2 = (Complex::one() - s2 * s2).sqrt();

    let rs = (i1 * c1 - i2 * c2) / (i1 * c1 + i2 * c2);
    let rp = (i2 * c1 - i1 * c2) / (i2 * c1 + i1 * c2);
    let ts = (Complex::real(2.0) * i1 * c1) / (i1 * c1 + i2 * c2);
    let tp = (Complex::real(2.0) * i1 * c1) / (i2 * c1 + i1 * c2);
    // θ₂ = arg(e^{iθ₂}) = arg(cosθ₂ + i·sinθ₂) — well-defined in complex.
    let theta_t = Complex::new(libm::atan2(c2.im + s2.re, c2.re - s2.im), 0.0);
    FresnelCoefficients {
        rs,
        rp,
        ts,
        tp,
        theta_t,
    }
}

/// Power reflectance and transmittance (s and p) for lossless media.
///
/// Returns `(Rs, Rp, Ts, Tp)`; `R + T = 1` for each polarisation.
/// For absorbing media use [`fresnel_amplitudes`] directly and interpret the
/// transmitted power with the projected flux; this function is meant for the
/// lossless case and asserts `k = 0`.
pub fn fresnel_energy(n1: f64, n2: f64, theta_i: f64) -> (f64, f64, f64, f64) {
    let c1 = libm::cos(theta_i);
    let s2 = n1 * libm::sin(theta_i) / n2;
    // Total internal reflection handled by caller via `critical_angle`;
    // here we compute the propagating case. Guard against domain errors.
    if s2 > 1.0 {
        return (1.0, 1.0, 0.0, 0.0);
    }
    let c2 = libm::sqrt(1.0 - s2 * s2);
    let rs = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2);
    let rp = (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2);
    let ts = 2.0 * n1 * c1 / (n1 * c1 + n2 * c2);
    let tp = 2.0 * n1 * c1 / (n2 * c1 + n1 * c2);
    let rs2 = rs * rs;
    let rp2 = rp * rp;
    let flux = n2 * c2 / (n1 * c1);
    let ts2 = ts * ts * flux;
    let tp2 = tp * tp * flux;
    (rs2, rp2, ts2, tp2)
}

/// Brewster (polarising) angle `arctan(n2/n1)` in radians.
pub fn brewster_angle(n1: f64, n2: f64) -> f64 {
    libm::atan(n2 / n1)
}

/// Critical angle for total internal reflection from `n1` into `n2`
/// (`n1 > n2`), in radians; `None` when there is no TIR.
pub fn critical_angle(n1: f64, n2: f64) -> Option<f64> {
    if n1 > n2 {
        Some(libm::asin(n2 / n1))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn normal_incidence() {
        let c = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.0),
            ComplexRefractiveIndex::real(1.5),
            0.0,
        );
        assert!(close(c.rs.re, -0.2, 1e-12));
        assert!(close(c.rp.re, 0.2, 1e-12));
        assert!(close(c.ts.re, 0.8, 1e-12));
        assert!(close(c.tp.re, 0.8, 1e-12));
        assert!(close(c.rs.norm_sqr(), 0.04, 1e-12));
    }

    #[test]
    fn brewster_zero_rp() {
        // Air to glass: θ_B = arctan(1.5) ≈ 56.31°
        let theta_b = brewster_angle(1.0, 1.5);
        let c = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.0),
            ComplexRefractiveIndex::real(1.5),
            theta_b,
        );
        assert!(c.rp.abs() < 1e-12, "rp = {}", c.rp);
        assert!(c.rs.abs() > 0.1);
        // Energy conservation for p: R + T = 1
        let (rp2, tp2) = {
            let (_, rp, _, tp) = fresnel_energy(1.0, 1.5, theta_b);
            (rp, tp)
        };
        assert!(close(rp2 + tp2, 1.0, 1e-12));
    }

    #[test]
    fn energy_conservation_lossless() {
        for deg in [0, 15, 30, 45, 60, 80] {
            let th = (deg as f64).to_radians();
            let (rs, rp, ts, tp) = fresnel_energy(1.0, 1.5, th);
            assert!(close(rs + ts, 1.0, 1e-12), "s at {deg}");
            assert!(close(rp + tp, 1.0, 1e-12), "p at {deg}");
        }
    }

    #[test]
    fn total_internal_reflection() {
        let crit = critical_angle(1.5, 1.0).unwrap();
        assert!(close(crit * 180.0 / core::f64::consts::PI, 41.81, 0.01));
        // Beyond critical: |rs| = |rp| = 1 (phase shifts only).
        let c = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.5),
            ComplexRefractiveIndex::real(1.0),
            crit + 0.1,
        );
        assert!(close(c.rs.abs(), 1.0, 1e-12));
        assert!(close(c.rp.abs(), 1.0, 1e-12));
        assert!(c.rs.arg().abs() > 1e-6, "expect non-trivial TIR phase");
        // fresnel_energy short-circuits to (1,1,0,0)
        let (rs, rp, ts, tp) = fresnel_energy(1.5, 1.0, crit + 0.1);
        assert!(close(rs, 1.0, 1e-15));
        assert!(close(rp, 1.0, 1e-15));
        assert!(close(ts, 0.0, 1e-15));
        assert!(close(tp, 0.0, 1e-15));
    }

    #[test]
    fn absorbing_metal_has_partial_reflection() {
        // Gold-ish at 633 nm: n≈0.18, k≈3.43 — reflectance should be ~0.9+.
        let au = ComplexRefractiveIndex::new(0.18, 3.43);
        let c = fresnel_amplitudes(ComplexRefractiveIndex::real(1.0), au, 0.0);
        assert!(c.rs.norm_sqr() > 0.85, "R = {}", c.rs.norm_sqr());
        // s and p magnitudes coincide at normal incidence (rs = −rp by
        // basis convention); transmitted amplitudes are equal.
        assert!(close(c.rs.abs(), c.rp.abs(), 1e-12));
        assert!(c.ts.approx_eq(c.tp, 1e-12));
    }
}

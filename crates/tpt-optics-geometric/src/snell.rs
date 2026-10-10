//! Vector Snell's law, reflection, and diffractive momentum kicks.
//!
//! All functions take a unit surface normal `n` oriented **against** the
//! incident direction (`dot(d, n) < 0`), as produced by
//! [`crate::surface::Hit`].

use tpt_optics_core::{Scalar, Vec2, Vec3};

/// Reflected direction (law of reflection).
pub fn reflect<T: Scalar>(d: Vec3<T>, n: Vec3<T>) -> Vec3<T> {
    // d − 2(d·n)n
    d - n * (T::from_f64(2.0) * d.dot(n))
}

/// Refracted direction via the vector Snell's law.
///
/// Returns `None` on total internal reflection. The result is a unit vector
/// (up to floating-point rounding).
pub fn refract<T: Scalar>(d: Vec3<T>, n: Vec3<T>, n1: f64, n2: f64) -> Option<Vec3<T>> {
    let eta = T::from_f64(n1 / n2);
    let cos_i = -d.dot(n);
    let k = T::ONE - eta * eta * (T::ONE - cos_i * cos_i);
    if k < T::ZERO {
        return None; // total internal reflection
    }
    let t = d * eta + n * (eta * cos_i - k.sqrt());
    t.normalized()
}

/// Deflection by a diffractive (phase) element.
///
/// Conserves tangential optical momentum including the phase kick:
/// `n₂·s₂,ₜ = n₁·s₁,ₜ + (λ₀/2π)·∇Φ`. Returns the transmitted unit
/// direction, or `None` when the diffracted order is evanescent.
///
/// `phase_gradient` is the in-plane phase gradient ∇Φ in rad/m at the hit
/// point (defined in surface coordinates; its normal component is projected
/// out here).
pub fn diffract<T: Scalar>(
    d: Vec3<T>,
    n: Vec3<T>,
    n1: f64,
    n2: f64,
    wavelength_m: f64,
    phase_gradient: Vec2<T>,
) -> Option<Vec3<T>> {
    let kick = T::from_f64(wavelength_m / core::f64::consts::TAU);
    let grad3 = Vec3::new(phase_gradient.x, phase_gradient.y, T::ZERO);
    let cos_i = d.dot(n); // negative under the against-normal convention
                          // Tangential components (v − (v·n)n).
    let d_t = d - n * cos_i;
    let grad_t = grad3 - n * grad3.dot(n);
    // Conserved tangential optical momentum, magnitude n2·sinθ2.
    let p = d_t * T::from_f64(n1) + grad_t * kick;
    let pz2 = T::from_f64(n2 * n2) - p.norm_sqr();
    if pz2 < T::ZERO {
        return None; // evanescent order
    }
    // Continue to the same side the incident ray was heading.
    let sign = if cos_i < T::ZERO { -T::ONE } else { T::ONE };
    let out = p + n * (pz2.sqrt() * sign);
    Some(out / T::from_f64(n2))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn reflection_law() {
        let d = Vec3::new(0.0, 0.0, -1.0);
        let n = Vec3::new(0.0, 0.0, 1.0);
        let r = reflect(d, n);
        assert!(r.approx_eq(Vec3::new(0.0, 0.0, 1.0), 1e-15));
    }

    #[test]
    fn snell_normal_incidence() {
        let d = Vec3::new(0.0, 0.0, -1.0);
        let n = Vec3::new(0.0, 0.0, 1.0);
        let t = refract(d, n, 1.0, 1.5).unwrap();
        // Transmits straight through (against normal): z component −1/n·n...
        // The refracted ray continues along +z crossing the surface:
        assert!(t.approx_eq(Vec3::new(0.0, 0.0, -1.0), 1e-12));
    }

    #[test]
    fn snell_45_degrees() {
        // 45° in air into glass: sin θt = sin45/1.5 → θt ≈ 28.126°
        let s = core::f64::consts::FRAC_1_SQRT_2;
        let d = Vec3::new(s, 0.0, -s); // incident at 45° on the z=0 plane
        let n = Vec3::new(0.0, 0.0, 1.0);
        let t = refract(d, n, 1.0, 1.5).unwrap();
        let theta_t = libm::acos(-t.z);
        assert!(
            close(theta_t.to_degrees(), 28.1255, 1e-3),
            "{}",
            theta_t.to_degrees()
        );
        // Tangential momentum: n1 sinθ1 = n2 sinθ2
        assert!(close(t.x * 1.5, s, 1e-12));
        assert!(close(t.norm(), 1.0, 1e-12));
    }

    #[test]
    fn tir_beyond_critical() {
        // Glass to air at 45° (> 41.8° critical): TIR.
        let s = core::f64::consts::FRAC_1_SQRT_2;
        let d = Vec3::new(s, 0.0, -s);
        let n = Vec3::new(0.0, 0.0, 1.0);
        assert!(refract(d, n, 1.5, 1.0).is_none());
        // Below critical (30°): refracts at sin = 1.5·0.5 = 0.75
        let d2 = Vec3::new(0.5, 0.0, -(3.0_f64).sqrt() / 2.0);
        let t = refract(d2, n, 1.5, 1.0).unwrap();
        assert!(close(t.x, 0.75, 1e-12));
    }

    #[test]
    fn diffraction_grating_deflection() {
        // Linear phase Φ = 2πx/Λ: transverse momentum kick λ/Λ.
        // Normal incidence, Λ = 2λ → sinθ2 = λ/Λ = 0.5 (30°).
        let d = Vec3::new(0.0, 0.0, -1.0);
        let n = Vec3::new(0.0, 0.0, 1.0);
        let lam = 1e-6;
        let grad_x = core::f64::consts::TAU / (2.0 * lam); // rad/m
        let t = diffract(d, n, 1.0, 1.0, lam, Vec2::new(grad_x, 0.0)).unwrap();
        assert!(close(t.x, 0.5, 1e-12), "tx = {}", t.x);
        assert!(close(t.z, -(3.0_f64).sqrt() / 2.0, 1e-12));
        assert!(close(t.norm(), 1.0, 1e-12));
        // Above the propagating-order limit (Λ < λ at normal incidence):
        let grad_x2 = core::f64::consts::TAU / (0.5 * lam);
        assert!(diffract(d, n, 1.0, 1.0, lam, Vec2::new(grad_x2, 0.0)).is_none());
    }
}

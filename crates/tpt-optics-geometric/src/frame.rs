//! Frame transforms: placing surfaces and systems in a global frame.
//!
//! Optical invariants (optical path length, Snell's law) are frame
//! independent — the property tests in `tpt-optics-verify` exploit that by
//! tracing the same system in translated/rotated frames.

use tpt_optics_core::{Scalar, Vec3};

use crate::ray::Ray;

/// Rigid placement: rotation `R` followed by translation `t`
/// (`p_global = R·p_local + t`).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Frame<T> {
    /// Rotation applied to local coordinates.
    pub rotation: Vec3<T>, // Euler angles (rad), Rz·Ry·Rx convention
    /// Translation added after rotation.
    pub translation: Vec3<T>,
}

impl<T: Scalar> Frame<T> {
    /// Identity frame (local == global).
    pub fn identity() -> Self {
        Frame {
            rotation: Vec3::zero(),
            translation: Vec3::zero(),
        }
    }

    /// Pure translation.
    pub fn translation(t: Vec3<T>) -> Self {
        Frame {
            rotation: Vec3::zero(),
            translation: t,
        }
    }

    /// Transform a point into this frame.
    pub fn apply_point(&self, p: Vec3<T>) -> Vec3<T> {
        rotate_euler(self.rotation, p) + self.translation
    }

    /// Transform a direction into this frame (no translation).
    pub fn apply_direction(&self, d: Vec3<T>) -> Vec3<T> {
        rotate_euler(self.rotation, d)
    }

    /// Transform a ray into this frame.
    pub fn apply_ray(&self, ray: &Ray<T>) -> Ray<T> {
        Ray::new(
            self.apply_point(ray.origin),
            self.apply_direction(ray.direction),
        )
    }
}

/// Rotate `v` by Euler angles (ψ, θ, φ) with convention R = Rz(ψ)·Ry(θ)·Rx(φ).
pub fn rotate_euler<T: Scalar>(angles: Vec3<T>, v: Vec3<T>) -> Vec3<T> {
    let (psi, theta, phi) = (angles.x, angles.y, angles.z);
    let (cp, sp) = (psi.cos(), psi.sin());
    let (ct, st) = (theta.cos(), theta.sin());
    let (cf, sf) = (phi.cos(), phi.sin());
    // Rx(φ)
    let v1 = Vec3::new(v.x, v.y * cf - v.z * sf, v.y * sf + v.z * cf);
    // Ry(θ)
    let v2 = Vec3::new(v1.x * ct + v1.z * st, v1.y, -v1.x * st + v1.z * ct);
    // Rz(ψ)
    Vec3::new(v2.x * cp - v2.y * sp, v2.x * sp + v2.y * cp, v2.z)
}

/// Inverse of [`rotate_euler`]: rotate `v` back out of a frame whose
/// rotation is the Euler triple `angles`.
pub fn unrotate_euler<T: Scalar>(angles: Vec3<T>, v: Vec3<T>) -> Vec3<T> {
    // R = Rz(ψ)·Ry(θ)·Rx(φ)  ⇒  R⁻¹ = Rx(−φ)·Ry(−θ)·Rz(−ψ).
    let (psi, theta, phi) = (angles.x, angles.y, angles.z);
    // Rz(−ψ)
    let (cp, sp) = ((-psi).cos(), (-psi).sin());
    let v1 = Vec3::new(v.x * cp - v.y * sp, v.x * sp + v.y * cp, v.z);
    // Ry(−θ)
    let (ct, st) = ((-theta).cos(), (-theta).sin());
    let v2 = Vec3::new(v1.x * ct + v1.z * st, v1.y, -v1.x * st + v1.z * ct);
    // Rx(−φ)
    let (cf, sf) = ((-phi).cos(), (-phi).sin());
    Vec3::new(v2.x, v2.y * cf - v2.z * sf, v2.y * sf + v2.z * cf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translation_only() {
        let f = Frame::translation(Vec3::new(1.0, 2.0, 3.0));
        let p = f.apply_point(Vec3::new(0.5, 0.5, 0.5));
        assert!(p.approx_eq(Vec3::new(1.5, 2.5, 3.5), 1e-15));
    }

    #[test]
    fn rotation_preserves_norm() {
        let f = Frame {
            rotation: Vec3::new(0.3, -1.1, 2.0),
            translation: Vec3::zero(),
        };
        let d = Vec3::new(0.2, -0.5, 0.84);
        let rotated = f.apply_direction(d);
        assert!((rotated.norm() - d.norm()).abs() < 1e-12);
    }

    #[test]
    fn axis_rotations() {
        // Rx(90°): y → z
        let v = rotate_euler(
            Vec3::new(0.0, 0.0, core::f64::consts::FRAC_PI_2),
            Vec3::new(0.0, 1.0, 0.0),
        );
        assert!(v.approx_eq(Vec3::new(0.0, 0.0, 1.0), 1e-12));
    }
}

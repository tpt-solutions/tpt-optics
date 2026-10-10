//! Rays: oriented lines carrying position and propagation direction.

use tpt_optics_core::{Scalar, Vec3};

/// A geometric ray: an origin plus a (not necessarily unit) direction.
///
/// Sequential and non-sequential tracers normalise `direction` on entry;
/// keeping it arbitrary simplifies construction. Optical path length is
/// measured with unit direction.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Ray<T> {
    /// Point where the ray starts.
    pub origin: Vec3<T>,
    /// Direction of propagation (normalised by [`Ray::normalized`]).
    pub direction: Vec3<T>,
}

impl<T: Scalar> Ray<T> {
    /// Construct.
    pub fn new(origin: Vec3<T>, direction: Vec3<T>) -> Self {
        Ray { origin, direction }
    }

    /// A ray starting at `origin` pointing at `target`.
    pub fn from_to(origin: Vec3<T>, target: Vec3<T>) -> Option<Self> {
        let d = target - origin;
        d.normalized().map(|d| Ray {
            origin,
            direction: d,
        })
    }

    /// Copy with a unit direction; `None` for a zero direction.
    pub fn normalized(&self) -> Option<Self> {
        self.direction
            .normalized()
            .map(|d| Ray::new(self.origin, d))
    }

    /// Point at parameter `t`: `origin + t·direction`.
    pub fn at(&self, t: T) -> Vec3<T> {
        self.origin.mul_add(T::ONE, self.direction, t)
    }

    /// Translate the origin without touching the direction.
    pub fn translated(&self, offset: Vec3<T>) -> Self {
        Ray::new(self.origin + offset, self.direction)
    }

    /// Rebase the origin to `point` (typically a hit point on a surface).
    pub fn from_point(&self, point: Vec3<T>) -> Self {
        Ray::new(point, self.direction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_at() {
        let r = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        assert!(r.at(2.5).approx_eq(Vec3::new(0.0, 0.0, 2.5), 1e-15));
    }

    #[test]
    fn normalize() {
        let r = Ray::new(Vec3::zero(), Vec3::new(0.0, 3.0, 4.0));
        let n = r.normalized().unwrap();
        assert!((n.direction.norm() - 1.0).abs() < 1e-15);
        assert!(Ray::new(Vec3::<f64>::zero(), Vec3::zero())
            .normalized()
            .is_none());
    }

    #[test]
    fn from_to() {
        let r = Ray::from_to(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 1.0)).unwrap();
        assert!((r.direction.norm() - 1.0).abs() < 1e-15);
    }
}

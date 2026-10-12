//! Optical surfaces: sag functions, normals and ray intersection.
//!
//! Every surface is described relative to its own vertex frame: the vertex
//! sits at the origin and the optical axis is +z. Sag is positive toward +z.
//!
//! The [`OpticalSurface`] trait is generic over the scalar type, which is
//! what makes the whole ray-trace differentiable: call `sag`/`sag_gradient`
//! with [`Dual`](tpt_optics_core::Dual) coordinates and the derivatives come
//! out exactly.

use alloc::vec::Vec;

use libm::pow;
use tpt_optics_core::{Scalar, Vec2, Vec3};

use crate::ray::Ray;

/// Local frame data attached to a hit: parameter, point and surface normal.
///
/// The normal points toward −z (against incoming rays propagating along +z),
/// i.e. it is oriented against the ray direction as the Snell code expects.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Hit<T> {
    /// Ray parameter at the intersection (`point = ray.at(t)`).
    pub t: T,
    /// Intersection point.
    pub point: Vec3<T>,
    /// Unit normal oriented against the incoming ray.
    pub normal: Vec3<T>,
}

impl<T: Scalar> Hit<T> {
    /// Build a hit from a ray parameter and an upward (+z) normal,
    /// flipping the normal if needed so it opposes the ray.
    pub fn oriented(ray: &Ray<T>, t: T, point: Vec3<T>, upward: Vec3<T>) -> Self {
        let normal = if ray.direction.dot(upward) > T::ZERO {
            -upward
        } else {
            upward
        };
        Hit { t, point, normal }
    }
}

/// A refracting/reflecting surface in its vertex frame.
pub trait OpticalSurface: core::fmt::Debug + Send + Sync {
    /// Sag height z(x, y).
    fn sag<T: Scalar>(&self, x: T, y: T) -> T;

    /// Analytic sag gradient (∂z/∂x, ∂z/∂y).
    fn sag_gradient<T: Scalar>(&self, x: T, y: T) -> Vec2<T>;

    /// Upward unit normal at (x, y): `(−z_x, −z_y, 1)` normalised.
    fn upward_normal<T: Scalar>(&self, x: T, y: T) -> Vec3<T> {
        let g = self.sag_gradient(x, y);
        Vec3::new(-g.x, -g.y, T::ONE)
            .normalized()
            .unwrap_or_else(|| Vec3::new(T::ZERO, T::ZERO, T::ONE))
    }

    /// Nearest intersection of the ray with this surface with `t > t_min`.
    ///
    /// The default implementation seeds from the vertex tangent plane and
    /// refines `F(t) = ray.z(t) − sag(ray.x(t), ray.y(t))` with
    /// Newton–Raphson iterations using the analytic sag gradient. Surfaces
    /// with closed-form quadric intersections (plane, sphere, conic)
    /// override this for exactness and speed.
    fn intersect<T: Scalar>(&self, ray: &Ray<T>, t_min: T) -> Option<Hit<T>> {
        // Seed: tangent plane at vertex (z = 0).
        if ray.direction.z == T::ZERO {
            return None;
        }
        let mut t = -ray.origin.z / ray.direction.z;
        if !t.is_finite() || t <= t_min {
            return None;
        }
        for _ in 0..24 {
            let p = ray.at(t);
            let f = p.z - self.sag(p.x, p.y);
            let g = self.sag_gradient(p.x, p.y);
            // F'(t) = dz_ray/dt − (z_x·dx/dt + z_y·dy/dt)
            let df = ray.direction.z - (g.x * ray.direction.x + g.y * ray.direction.y);
            if df == T::ZERO {
                break;
            }
            let step = f / df;
            t -= step;
            if !t.is_finite() {
                return None;
            }
            if T::to_f64(step).abs() < 1e-15 {
                break;
            }
        }
        if t <= t_min {
            return None;
        }
        let point = ray.at(t);
        // Residual check: reject false convergence.
        let residual = point.z - self.sag(point.x, point.y);
        if T::to_f64(residual.abs()) > 1e-7 {
            return None;
        }
        let upward = self.upward_normal(point.x, point.y);
        Some(Hit::oriented(ray, t, point, upward))
    }
}

/// An infinite plane at z = 0 in its vertex frame.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Plane;

impl OpticalSurface for Plane {
    fn sag<T: Scalar>(&self, _x: T, _y: T) -> T {
        T::ZERO
    }

    fn sag_gradient<T: Scalar>(&self, _x: T, _y: T) -> Vec2<T> {
        Vec2::new(T::ZERO, T::ZERO)
    }

    fn intersect<T: Scalar>(&self, ray: &Ray<T>, t_min: T) -> Option<Hit<T>> {
        if ray.direction.z == T::ZERO {
            return None;
        }
        let t = -ray.origin.z / ray.direction.z;
        if t <= t_min {
            return None;
        }
        Some(Hit::oriented(
            ray,
            t,
            ray.at(t),
            Vec3::new(T::ZERO, T::ZERO, T::ONE),
        ))
    }
}

/// A spherical surface: vertex at the origin, curvature `c = 1/R`.
///
/// `c > 0` is convex toward +z (center of curvature at z = R). This is the
/// sign convention of most lens-design codes.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Sphere {
    /// Curvature 1/R (1/m).
    pub curvature: f64,
}

impl Sphere {
    /// Construct from curvature.
    pub const fn new(curvature: f64) -> Self {
        Sphere { curvature }
    }

    /// Construct from radius of curvature.
    pub const fn from_radius(radius: f64) -> Self {
        Sphere {
            curvature: if radius == 0.0 { 0.0 } else { 1.0 / radius },
        }
    }
}

impl OpticalSurface for Sphere {
    fn sag<T: Scalar>(&self, x: T, y: T) -> T {
        let c = T::from_f64(self.curvature);
        let r2 = x * x + y * y;
        // Exact sag: (c·r²)/(1 + √(1 − c²r²)) == c r²/(1+sqrt(...)); stable
        // form avoids cancellation for large r.
        r2 * c / (T::ONE + (T::ONE - c * c * r2).sqrt())
    }

    fn sag_gradient<T: Scalar>(&self, x: T, y: T) -> Vec2<T> {
        let c = T::from_f64(self.curvature);
        let r2 = x * x + y * y;
        // dz/dr = c·r / sqrt(1 − c²r²); multiply by (x, y)/r.
        let denom = (T::ONE - c * c * r2).sqrt();
        Vec2::new(c * x / denom, c * y / denom)
    }

    fn intersect<T: Scalar>(&self, ray: &Ray<T>, t_min: T) -> Option<Hit<T>> {
        if self.curvature == 0.0 {
            return Plane.intersect(ray, t_min);
        }
        // Quadric |p − C|² = R² with C = (0, 0, R), vertex at the origin.
        let o = ray.origin;
        let d = ray.direction;
        let center = Vec3::new(T::ZERO, T::ZERO, T::from_f64(1.0 / self.curvature));
        let oc = o - center;
        let a = d.norm_sqr();
        let b = T::from_f64(2.0) * oc.dot(d);
        let c0 = oc.norm_sqr() - T::from_f64(1.0 / (self.curvature * self.curvature));
        let disc = b * b - T::from_f64(4.0) * a * c0;
        if disc < T::ZERO || a == T::ZERO {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let inv = T::ONE / (T::from_f64(2.0) * a);
        let t1 = (-b - sqrt_disc) * inv;
        let t2 = (-b + sqrt_disc) * inv;
        // Nearest positive root (t1 <= t2 for a > 0).
        let t = if t1 > t_min {
            t1
        } else if t2 > t_min {
            t2
        } else {
            return None;
        };
        let point = ray.at(t);
        // Radial normal through the center; orient with +z sense.
        let mut up = point - center;
        if up.z < T::ZERO {
            up = -up;
        }
        let up = up.normalized()?;
        Some(Hit::oriented(ray, t, point, up))
    }
}

/// Rotationally symmetric conic/even-asphere surface:
///
/// `z(r) = c·r²/(1 + √(1 − (1+k)·c²·r²)) + Σᵢ aᵢ·r^{2i}`, i = 1..8.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Asphere {
    /// Base curvature 1/R (1/m).
    pub curvature: f64,
    /// Conic constant k (k = −e² for a prolate ellipse of eccentricity e).
    pub conic: f64,
    /// Even polynomial coefficients a₁..a₈ (units m^{−(2i−1)}).
    pub coeffs: [f64; 8],
}

impl Asphere {
    /// A conic surface without polynomial terms.
    pub const fn conic(curvature: f64, conic: f64) -> Self {
        Asphere {
            curvature,
            conic,
            coeffs: [0.0; 8],
        }
    }

    /// A pure sphere.
    pub const fn sphere(curvature: f64) -> Self {
        Asphere::conic(curvature, 0.0)
    }

    /// A plane.
    pub const fn flat() -> Self {
        Asphere::conic(0.0, 0.0)
    }

    /// Paraxial equivalence: for small r, z ≈ c·r²/2, the same as a sphere.
    pub fn effective_curvature(&self) -> f64 {
        self.curvature
    }
}

impl OpticalSurface for Asphere {
    fn sag<T: Scalar>(&self, x: T, y: T) -> T {
        let c = T::from_f64(self.curvature);
        let k = T::from_f64(self.conic);
        let r2 = x * x + y * y;
        let denom = T::ONE + (T::ONE - (T::ONE + k) * c * c * r2).sqrt();
        let mut z = c * r2 / denom;
        let mut rpow = r2;
        for a in &self.coeffs {
            z += T::from_f64(*a) * rpow;
            rpow *= r2;
        }
        z
    }

    fn sag_gradient<T: Scalar>(&self, x: T, y: T) -> Vec2<T> {
        let c = T::from_f64(self.curvature);
        let k = T::from_f64(self.conic);
        let r2 = x * x + y * y;
        // Conic part: ∇z = (c/denom)·(x, y), the singular-free form of
        // (c·r/denom)·(x, y)/r, where denom = √(1 − (1+k)c²r²).
        let denom = (T::ONE - (T::ONE + k) * c * c * r2).sqrt();
        let mut gx = c * x / denom;
        let mut gy = c * y / denom;
        // Polynomial part: ∂/∂x Σ aᵢ r^{2i} = Σ 2i·aᵢ·r^{2i−2}·x.
        let mut r2pow = T::ONE; // r^{2i−2}
        for (i, a) in self.coeffs.iter().enumerate() {
            let coef = T::from_f64(2.0 * (i + 1) as f64 * *a);
            gx += coef * r2pow * x;
            gy += coef * r2pow * y;
            r2pow *= r2;
        }
        Vec2::new(gx, gy)
    }
}

/// Zernike freeform surface: base conic + a sparse sum of standard Zernike
/// polynomial terms (OSA/ANSI single-index j, Fringe-independent).
///
/// `z = z_base(x,y) + Σ_j c_j·Z_j(x/a, y/a)` with normalization
/// `Z_j = N·R_n^m(ρ)·{cos mθ, sin mθ}` and `N = √(2(n+1))` (1 for m = 0).
#[derive(Clone, Debug, PartialEq)]
pub struct Zernike {
    /// Base rotationally symmetric surface.
    pub base: Asphere,
    /// Aperture semi-diameter a (meters) over which ρ = r/a.
    pub aperture: f64,
    /// (OSA index j, coefficient in meters), j ≥ 1 (j = 0 is piston).
    pub terms: Vec<(u32, f64)>,
}

/// Decompose OSA/ANSI index j into radial order n and azimuthal m.
#[must_use]
pub fn zernike_nm(j: u32) -> (u32, i32) {
    let mut n = 0u32;
    let mut j_rem = j;
    while j_rem > n {
        j_rem -= n + 1;
        n += 1;
    }
    let m = 2 * j_rem as i32 - n as i32;
    (n, m)
}

/// Radial polynomial coefficients of R_n^m(ρ) in ascending powers of ρ²
/// (i.e. powers n−2s for s = 0..=(n−m)/2).
#[must_use]
pub fn zernike_radial_coeffs(n: u32, m: u32) -> Vec<f64> {
    // Undefined unless m ≤ n and n − m is even (R_n^m has no such term).
    if m > n || (n - m) % 2 != 0 {
        return Vec::new();
    }
    // |m| — sign handled by cos/sin symmetry.
    let k_max = (n - m) / 2;
    let mut coeffs = Vec::with_capacity(k_max as usize + 1);
    for s in 0..=k_max {
        let sign = if s % 2 == 0 { 1.0 } else { -1.0 };
        let denom = factorial(s) * factorial((n + m) / 2 - s) * factorial((n - m) / 2 - s);
        coeffs.push(sign * factorial(n - s) / denom);
    }
    coeffs
}

/// `k!` as `f64`, computed directly so there is no fixed table to overflow.
fn factorial(k: u32) -> f64 {
    (1..=k).fold(1.0, |acc, i| acc * i as f64)
}

/// Evaluate the standard normalized Zernike j (OSA/ANSI) at (ρ, θ).
///
/// Reference (non-differentiable) version for tests and plotting.
pub fn zernike_eval(j: u32, rho: f64, theta: f64) -> f64 {
    let (n, m) = zernike_nm(j);
    let am = m.unsigned_abs();
    let radial: f64 = zernike_radial_coeffs(n, am)
        .iter()
        .enumerate()
        .map(|(s, c)| c * pow(rho, (n - 2 * s as u32) as f64))
        .sum();
    let norm2 = if m == 0 {
        (n + 1) as f64
    } else {
        2.0 * (n + 1) as f64
    };
    let ang = if m >= 0 {
        (m as f64 * theta).cos()
    } else {
        (am as f64 * theta).sin()
    };
    norm2.sqrt() * radial * ang
}

/// Radial polynomial R_n^m(ρ) as a generic-Scalar power sum.
fn zernike_radial<T: Scalar>(n: u32, m: u32, rho: T) -> T {
    let mut sum = T::ZERO;
    for (s, c) in zernike_radial_coeffs(n, m).iter().enumerate() {
        let power = (n - 2 * s as u32) as i32;
        sum += T::from_f64(*c) * rho.powi(power);
    }
    sum
}

/// Derivative of R_n^m(ρ) wrt ρ, generic over Scalar.
fn zernike_radial_prime<T: Scalar>(n: u32, m: u32, rho: T) -> T {
    let mut sum = T::ZERO;
    for (s, c) in zernike_radial_coeffs(n, m).iter().enumerate() {
        let power = (n - 2 * s as u32) as i32;
        if power != 0 {
            sum += T::from_f64(*c * power as f64) * rho.powi(power - 1);
        }
    }
    sum
}

/// (Z, ∂Z/∂x, ∂Z/∂y) of one Zernike term, generic over Scalar.
///
/// With ρ = r/a and θ = atan2(y, x):
/// `∂Z/∂x = N·[R'(ρ)·cosθ_r·(x/(r·a)) + R·(dAngle/dθ)·(−y/r²)]` and the
/// y-part symmetric.
fn zernike_term<T: Scalar>(n: u32, m: i32, aperture: f64, x: T, y: T) -> (T, T, T) {
    let a = T::from_f64(aperture);
    let r2 = x * x + y * y;
    let r = r2.sqrt();
    let rho = r / a;
    let theta = y.atan2(x);
    let am = m.unsigned_abs();

    let radial = zernike_radial(n, am, rho);
    let radial_p = zernike_radial_prime(n, am, rho);

    let norm2 = if m == 0 {
        (n + 1) as f64
    } else {
        2.0 * (n + 1) as f64
    };
    let norm = T::from_f64(norm2.sqrt());

    // Angular factor and its θ-derivative (cos for m ≥ 0, sin for m < 0).
    let (angle, angle_dtheta) = if m >= 0 {
        let w = T::from_f64(m as f64) * theta;
        (w.cos(), -T::from_f64(m as f64) * w.sin())
    } else {
        let w = T::from_f64(am as f64) * theta;
        (w.sin(), T::from_f64(am as f64) * w.cos())
    };

    let z = norm * radial * angle;
    let zx = norm * (radial_p * angle * x / (r * a) - radial * angle_dtheta * y / r2);
    let zy = norm * (radial_p * angle * y / (r * a) + radial * angle_dtheta * x / r2);
    (z, zx, zy)
}

impl OpticalSurface for Zernike {
    fn sag<T: Scalar>(&self, x: T, y: T) -> T {
        let mut z = self.base.sag(x, y);
        for (j, c) in &self.terms {
            let (n, m) = zernike_nm(*j);
            let (zj, _, _) = zernike_term(n, m, self.aperture, x, y);
            z += T::from_f64(*c) * zj;
        }
        z
    }

    fn sag_gradient<T: Scalar>(&self, x: T, y: T) -> Vec2<T> {
        let mut g = self.base.sag_gradient(x, y);
        if x == T::ZERO && y == T::ZERO {
            // On-axis limit: only first-order (tilt) terms contribute;
            // Z₁ = 2x/a and Z₂ = 2y/a.
            for (j, c) in &self.terms {
                let (n, m) = zernike_nm(*j);
                if n == 1 {
                    let g1 = T::from_f64(2.0 * c / self.aperture);
                    if m > 0 {
                        g.x += g1;
                    } else {
                        g.y += g1;
                    }
                }
            }
            return g;
        }
        for (j, c) in &self.terms {
            let (n, m) = zernike_nm(*j);
            let (_, zx, zy) = zernike_term(n, m, self.aperture, x, y);
            let cj = T::from_f64(*c);
            g.x += cj * zx;
            g.y += cj * zy;
        }
        g
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn sphere_sag_matches_geometry() {
        let s = Sphere::from_radius(10.0);
        // Sag at r = 2: R − √(R² − r²) = 10 − √96
        let expect = 10.0 - (96.0_f64).sqrt();
        assert!(close(s.sag(2.0, 0.0), expect, 1e-12));
        assert!(close(s.sag(0.0, 2.0), expect, 1e-12));
    }

    #[test]
    fn sphere_intersect_axial() {
        let s = Sphere::from_radius(10.0);
        let ray = Ray::new(Vec3::new(0.0, 0.0, -5.0), Vec3::new(0.0, 0.0, 1.0));
        let hit = s.intersect(&ray, 1e-9).unwrap();
        // Front surface at z = 10 − √(100−0) = 0 for r=0 → vertex at z=0,
        // ray from z=−5 hits at z=0, t=5.
        assert!(close(hit.point.z, 0.0, 1e-12));
    }

    #[test]
    fn sphere_intersect_off_axis() {
        let s = Sphere::from_radius(10.0);
        // Parallel ray at height 2 traveling +z.
        let ray = Ray::new(Vec3::new(2.0, 0.0, -5.0), Vec3::new(0.0, 0.0, 1.0));
        let hit = s.intersect(&ray, 1e-9).unwrap();
        let pz = 10.0 - (96.0_f64).sqrt();
        assert!(close(hit.point.z, pz, 1e-10));
        // Normal passes through the center.
        let center = Vec3::new(0.0, 0.0, 10.0);
        let dir = (hit.point - center).normalized().unwrap();
        assert!(dir.approx_eq(-hit.normal, 1e-12) || dir.approx_eq(hit.normal, 1e-12));
    }

    #[test]
    fn conic_sag_parabola_limit() {
        // k = −1 → paraboloid z = c r²/2
        let p = Asphere::conic(0.5, -1.0);
        assert!(close(p.sag(2.0, 0.0), 0.5 * 4.0 / 2.0, 1e-12));
    }

    #[test]
    fn asphere_poly_term() {
        let mut a = Asphere::flat();
        a.coeffs[0] = 1.0; // z = r²
        assert!(close(a.sag(2.0, 0.0), 4.0, 1e-12));
        let g = a.sag_gradient(2.0, 0.0);
        assert!(close(g.x, 4.0, 1e-9), "gx = {}", g.x);
        assert!(close(g.y, 0.0, 1e-12));
    }

    #[test]
    fn asphere_gradient_vs_fd() {
        let mut a = Asphere::conic(0.1, -2.5);
        a.coeffs[2] = 0.001;
        let (x, y) = (0.013, -0.004);
        let g = a.sag_gradient(x, y);
        let h = 1e-9;
        let zx = (a.sag(x + h, y) - a.sag(x - h, y)) / (2.0 * h);
        let zy = (a.sag(x, y + h) - a.sag(x, y - h)) / (2.0 * h);
        assert!((g.x - zx).abs() < 1e-6, "gx {} vs {}", g.x, zx);
        assert!((g.y - zy).abs() < 1e-6, "gy {} vs {}", g.y, zy);
    }

    #[test]
    fn zernike_indices() {
        // OSA indexing: j=4 → n=2,m=0 (defocus); j=5 → n=2,m=±2 45°;
        assert_eq!(zernike_nm(0), (0, 0));
        assert_eq!(zernike_nm(1), (1, -1));
        assert_eq!(zernike_nm(2), (1, 1));
        assert_eq!(zernike_nm(4), (2, 0));
    }

    #[test]
    fn zernike_defocus_is_quadratic() {
        // j=4: Z = √3·(2ρ²−1)
        for rho in [0.0, 0.5, 1.0] {
            let expect = (3.0_f64).sqrt() * (2.0 * rho * rho - 1.0);
            assert!(close(zernike_eval(4, rho, 0.7), expect, 1e-12));
        }
    }

    #[test]
    fn zernike_tilt() {
        // j=2 (n=1, m=+1): Z = 2·ρ·cosθ = 2x/a
        let z = zernike_eval(2, 0.6, 0.0);
        assert!(close(z, 2.0 * 0.6, 1e-12));
        // j=1 (n=1, m=−1): Z = 2·ρ·sinθ = 2y/a
        let z2 = zernike_eval(1, 0.6, core::f64::consts::FRAC_PI_2);
        assert!(close(z2, 2.0 * 0.6, 1e-12));
    }

    #[test]
    fn newton_intersection_general_surface() {
        // Biconic-ish asphere with polynomial terms: use default Newton.
        let mut a = Asphere::conic(0.05, -1.0);
        a.coeffs[1] = 100.0;
        let ray = Ray::new(Vec3::new(0.001, 0.0, -0.01), Vec3::new(0.0, 0.0, 1.0));
        let hit = a.intersect(&ray, 1e-12).unwrap();
        let residual = hit.point.z - a.sag(hit.point.x, hit.point.y);
        assert!(residual.abs() < 1e-12);
    }
}

#[cfg(test)]
mod zernike_range_tests {
    use super::*;

    #[test]
    fn radial_polynomial_is_one_at_the_rim_for_high_order() {
        // R_n^m(1) = 1 for every valid (n, m); n = 25 is beyond the old fixed
        // factorial table (which panicked above 20).
        for &(n, m) in &[(25u32, 3u32), (20, 0), (24, 6)] {
            let sum: f64 = zernike_radial_coeffs(n, m).iter().sum();
            assert!((sum - 1.0).abs() < 1e-6, "R_{n}^{m}(1) = {sum}");
        }
    }

    #[test]
    fn invalid_radial_order_is_empty_not_a_panic() {
        assert!(zernike_radial_coeffs(2, 4).is_empty());
        assert!(zernike_radial_coeffs(4, 1).is_empty());
    }
}

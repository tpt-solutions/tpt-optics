//! Small dense linear algebra generic over [`Scalar`].
//!
//! Only the shapes the optics kernels need are provided: 2/3-vectors, a 2×2
//! matrix (ABCD paraxial formalism), and helpers. Heavier algebra lives in
//! the specialised crates.

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::complex::Complex;
use crate::num::Scalar;

/// 2-vector.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vec2<T> {
    /// First component.
    pub x: T,
    /// Second component.
    pub y: T,
}

impl<T: Scalar> Vec2<T> {
    /// Construct.
    #[inline]
    pub const fn new(x: T, y: T) -> Self {
        Vec2 { x, y }
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, o: Self) -> T {
        self.x * o.x + self.y * o.y
    }

    /// Squared length.
    #[inline]
    pub fn norm_sqr(self) -> T {
        self.dot(self)
    }

    /// Length.
    #[inline]
    pub fn norm(self) -> T {
        self.norm_sqr().sqrt()
    }
}

impl<T: Scalar> Add for Vec2<T> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Vec2::new(self.x + o.x, self.y + o.y)
    }
}

impl<T: Scalar> Sub for Vec2<T> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Vec2::new(self.x - o.x, self.y - o.y)
    }
}

impl<T: Scalar> Mul<T> for Vec2<T> {
    type Output = Self;
    #[inline]
    fn mul(self, s: T) -> Self {
        Vec2::new(self.x * s, self.y * s)
    }
}

impl<T: Scalar> Div<T> for Vec2<T> {
    type Output = Self;
    #[inline]
    fn div(self, s: T) -> Self {
        let inv = T::ONE / s;
        Vec2::new(self.x * inv, self.y * inv)
    }
}

impl<T: Scalar> AddAssign for Vec2<T> {
    #[inline]
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl<T: Scalar> SubAssign for Vec2<T> {
    #[inline]
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o;
    }
}

/// 3-vector with the usual vector-algebra operations.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vec3<T> {
    /// X component.
    pub x: T,
    /// Y component.
    pub y: T,
    /// Z component.
    pub z: T,
}

impl Vec3<f64> {
    /// Approximate componentwise equality.
    pub fn approx_eq(self, o: Self, tol: f64) -> bool {
        (self.x - o.x).abs() < tol && (self.y - o.y).abs() < tol && (self.z - o.z).abs() < tol
    }
}

impl<T: Scalar> Vec3<T> {
    /// Construct.
    #[inline]
    pub const fn new(x: T, y: T, z: T) -> Self {
        Vec3 { x, y, z }
    }

    /// The zero vector.
    #[inline]
    pub fn zero() -> Self {
        Vec3::new(T::ZERO, T::ZERO, T::ZERO)
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, o: Self) -> T {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product `self × o`.
    #[inline]
    pub fn cross(self, o: Self) -> Self {
        Vec3::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Squared length.
    #[inline]
    pub fn norm_sqr(self) -> T {
        self.dot(self)
    }

    /// Length.
    #[inline]
    pub fn norm(self) -> T {
        self.norm_sqr().sqrt()
    }

    /// Unit-length copy; `None` for the zero vector.
    #[inline]
    pub fn normalized(self) -> Option<Self> {
        let n = self.norm();
        if n > T::ZERO {
            Some(self / n)
        } else {
            None
        }
    }

    /// Componentwise scaling by a scalar.
    #[inline]
    pub fn scale(self, s: T) -> Self {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }

    /// Linear combination `a·self + b·o`.
    #[inline]
    pub fn mul_add(self, a: T, o: Self, b: T) -> Self {
        Vec3::new(
            self.x * a + o.x * b,
            self.y * a + o.y * b,
            self.z * a + o.z * b,
        )
    }

    /// Squared Euclidean distance to `o`.
    #[inline]
    pub fn dist_sqr(self, o: Self) -> T {
        (self - o).norm_sqr()
    }

    /// Apply a scalar function to each component.
    #[inline]
    pub fn map<U: Scalar>(self, f: impl Fn(T) -> U) -> Vec3<U> {
        Vec3::new(f(self.x), f(self.y), f(self.z))
    }

    /// Componentwise negation (also available via [`Neg`]).
    #[inline]
    pub fn flip(self) -> Self {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

impl<T: Scalar> Add for Vec3<T> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl<T: Scalar> Sub for Vec3<T> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl<T: Scalar> Mul<T> for Vec3<T> {
    type Output = Self;
    #[inline]
    fn mul(self, s: T) -> Self {
        self.scale(s)
    }
}

impl<T: Scalar> Div<T> for Vec3<T> {
    type Output = Self;
    #[inline]
    fn div(self, s: T) -> Self {
        self.scale(T::ONE / s)
    }
}

impl<T: Scalar> Neg for Vec3<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        self.flip()
    }
}

impl<T: Scalar> AddAssign for Vec3<T> {
    #[inline]
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl<T: Scalar> SubAssign for Vec3<T> {
    #[inline]
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o;
    }
}

/// 2×2 matrix, used for the ABCD paraxial formalism and ray transfers.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix2<T> {
    /// Row 0, column 0.
    pub a: T,
    /// Row 0, column 1.
    pub b: T,
    /// Row 1, column 0.
    pub c: T,
    /// Row 1, column 1.
    pub d: T,
}

impl<T: Scalar> Matrix2<T> {
    /// Construct from row-major elements `[[a, b], [c, d]]`.
    #[inline]
    pub const fn new(a: T, b: T, c: T, d: T) -> Self {
        Matrix2 { a, b, c, d }
    }

    /// The 2×2 identity.
    #[inline]
    pub fn identity() -> Self {
        Matrix2::new(T::ONE, T::ZERO, T::ZERO, T::ONE)
    }

    /// Determinant `ad − bc`.
    #[inline]
    pub fn det(self) -> T {
        self.a * self.d - self.b * self.c
    }

    /// Matrix product.
    #[inline]
    pub fn mul_matrix(self, m: Self) -> Self {
        Matrix2::new(
            self.a * m.a + self.b * m.c,
            self.a * m.b + self.b * m.d,
            self.c * m.a + self.d * m.c,
            self.c * m.b + self.d * m.d,
        )
    }

    /// Apply to a ray vector `(y, u)` (height, angle).
    #[inline]
    pub fn apply_ray(self, y: T, u: T) -> (T, T) {
        (self.a * y + self.b * u, self.c * y + self.d * u)
    }

    /// Apply to a complex beam parameter: `q' = (a·q + b) / (c·q + d)`.
    ///
    /// This is the Gaussian-beam ABCD law.
    #[inline]
    pub fn apply_q(self, q: Complex<T>) -> Complex<T> {
        (Complex::real(self.a) * q + Complex::real(self.b))
            / (Complex::real(self.c) * q + Complex::real(self.d))
    }

    /// Fractional stability criterion for a round trip: `|trace| <= 2`.
    ///
    /// A resonator described by this round-trip matrix supports stable modes
    /// when `(a + d)/2` lies in `[-1, 1]`.
    #[inline]
    pub fn is_stable(self) -> bool {
        let half_trace = (self.a + self.d) * T::from_f64(0.5);
        half_trace.abs() <= T::ONE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_and_norm() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert!(x.cross(y).approx_eq(Vec3::new(0.0, 0.0, 1.0), 1e-15));
        assert!((x.cross(y).norm() - 1.0).abs() < 1e-15);
        assert!(x.dot(y).abs() < 1e-15);
    }

    #[test]
    fn matrix2_q_law_free_space() {
        // Free space: q -> q + d
        let m = Matrix2::<f64>::new(1.0, 2.5, 0.0, 1.0);
        let q = Complex::real(1.3);
        assert!((m.apply_q(q) - Complex::real(3.8)).approx_eq(Complex::zero(), 1e-14));
    }

    #[test]
    fn matrix2_thin_lens_q_law() {
        // Thin lens f: 1/q transforms with -1/f term; check against manual.
        let f = 0.1;
        let m = Matrix2::<f64>::new(1.0, 0.0, -1.0 / f, 1.0);
        let q = Complex::new(0.0, 1.0); // waist at lens
        let q2 = m.apply_q(q);
        // (i)/( -i/f * i + 1) = i/(1 - i/f)
        let expect = Complex::new(0.0, 1.0) / (Complex::new(1.0, 0.0) - Complex::new(0.0, 10.0));
        assert!(q2.approx_eq(expect, 1e-14));
    }

    #[test]
    fn stability() {
        assert!(Matrix2::<f64>::new(1.0, 0.1, -2.0, 1.0).is_stable());
        assert!(!Matrix2::<f64>::new(3.0, 0.0, 0.0, 1.0 / 3.0).is_stable());
    }
}

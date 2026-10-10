//! Complex numbers generic over [`Scalar`].

use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::num::Scalar;

/// A complex number `re + i·im`.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Complex<T> {
    /// Real part.
    pub re: T,
    /// Imaginary part.
    pub im: T,
}

impl<T: Scalar> Complex<T> {
    /// Construct from real and imaginary parts.
    #[inline]
    pub const fn new(re: T, im: T) -> Self {
        Complex { re, im }
    }

    /// The complex zero.
    #[inline]
    pub fn zero() -> Self {
        Complex::new(T::ZERO, T::ZERO)
    }

    /// The complex one.
    #[inline]
    pub fn one() -> Self {
        Complex::new(T::ONE, T::ZERO)
    }

    /// A purely real complex number.
    #[inline]
    pub fn real(re: T) -> Self {
        Complex::new(re, T::ZERO)
    }

    /// `e^{iθ}`.
    #[inline]
    pub fn cis(theta: T) -> Self {
        Complex::new(theta.cos(), theta.sin())
    }

    /// Complex conjugate.
    #[inline]
    pub fn conj(self) -> Self {
        Complex::new(self.re, -self.im)
    }

    /// Squared modulus `re² + im²` (differentiable where `abs` is not).
    #[inline]
    pub fn norm_sqr(self) -> T {
        self.re * self.re + self.im * self.im
    }

    /// Modulus.
    #[inline]
    pub fn abs(self) -> T {
        self.norm_sqr().sqrt()
    }

    /// Argument (phase angle).
    #[inline]
    pub fn arg(self) -> T {
        self.im.atan2(self.re)
    }

    /// Multiply by a real scalar.
    #[inline]
    pub fn scale(self, s: T) -> Self {
        Complex::new(self.re * s, self.im * s)
    }

    /// Principal square root with non-negative real part.
    #[inline]
    pub fn sqrt(self) -> Self {
        // Careful branch keeps `re >= 0` for the result.
        if self.re == T::ZERO && self.im == T::ZERO {
            return Self::zero();
        }
        let m = self.abs();
        let re = ((m + self.re) * T::from_f64(0.5)).sqrt();
        let im = ((m - self.re) * T::from_f64(0.5)).sqrt();
        if self.im < T::ZERO {
            Complex::new(re, -im)
        } else {
            Complex::new(re, im)
        }
    }

    /// Complex exponential.
    #[inline]
    pub fn exp(self) -> Self {
        let r = self.re.exp();
        Complex::new(r * self.im.cos(), r * self.im.sin())
    }

    /// Principal complex logarithm.
    #[inline]
    pub fn ln(self) -> Self {
        Complex::new(self.abs().ln(), self.arg())
    }

    /// Complex power `self^e` with real exponent.
    #[inline]
    pub fn powf(self, e: T) -> Self {
        (self.ln().scale(e)).exp()
    }

    /// Reciprocal `1/self`.
    #[inline]
    pub fn inv(self) -> Self {
        self.conj().scale(T::ONE / self.norm_sqr())
    }

    /// True if both parts are finite.
    #[inline]
    pub fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
}

impl Complex<f64> {
    /// Approximate equality with an absolute tolerance on both parts.
    pub fn approx_eq(self, other: Self, tol: f64) -> bool {
        (self.re - other.re).abs() < tol && (self.im - other.im).abs() < tol
    }
}

impl<T: Scalar> Add for Complex<T> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Complex::new(self.re + o.re, self.im + o.im)
    }
}

impl<T: Scalar> Sub for Complex<T> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Complex::new(self.re - o.re, self.im - o.im)
    }
}

impl<T: Scalar> Mul for Complex<T> {
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        Complex::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

impl<T: Scalar> Div for Complex<T> {
    type Output = Self;
    // `self * inv` is intentional: one multiplication instead of the
    // four-multiply/two-add division formula.
    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn div(self, o: Self) -> Self {
        self * o.inv()
    }
}

impl<T: Scalar> Neg for Complex<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Complex::new(-self.re, -self.im)
    }
}

impl<T: Scalar> AddAssign for Complex<T> {
    #[inline]
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl<T: Scalar> SubAssign for Complex<T> {
    #[inline]
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o;
    }
}

impl<T: Scalar> MulAssign for Complex<T> {
    #[inline]
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o;
    }
}

impl<T: Scalar> DivAssign for Complex<T> {
    #[inline]
    fn div_assign(&mut self, o: Self) {
        *self = *self / o;
    }
}

impl<T: Scalar> From<T> for Complex<T> {
    #[inline]
    fn from(re: T) -> Self {
        Complex::real(re)
    }
}

impl<T: Scalar + fmt::Display> fmt::Display for Complex<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}i)", self.re, self.im)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn arithmetic() {
        let a = Complex::new(3.0, 4.0);
        let b = Complex::new(1.0, -2.0);
        assert!(close((a * b).re, 11.0));
        assert!(close((a * b).im, -2.0));
        let q = a / b;
        assert!(close(q.re, -1.0));
        assert!(close(q.im, 2.0));
        assert!(close(a.abs(), 5.0));
        assert!(close(a.arg(), 4.0f64.atan2(3.0)));
    }

    #[test]
    fn sqrt_branch() {
        let m = Complex::new(-4.0, 0.0);
        assert!(close(m.sqrt().re, 0.0));
        assert!(close(m.sqrt().im, 2.0));
        let m2 = Complex::new(3.0, 4.0);
        assert!(close(m2.sqrt().re, 2.0));
        assert!(close(m2.sqrt().im, 1.0));
        // Principal root has non-negative real part everywhere.
        let z = Complex::new(-1.0, 1e-12);
        assert!(z.sqrt().re >= 0.0);
    }

    #[test]
    fn exp_log_roundtrip() {
        let z = Complex::new(0.3, 1.4);
        let w = z.exp().ln();
        assert!(w.approx_eq(z, 1e-12));
    }

    #[test]
    fn euler() {
        let c = Complex::cis(2.0_f64.atan2(1.0)); // unit circle point
        assert!(close(c.abs(), 1.0));
        assert!(close(c.re, 1.0 / 5.0f64.sqrt()));
    }
}

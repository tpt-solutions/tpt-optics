//! Generic scalar abstraction used to make optics code differentiable.
//!
//! Every geometric/wave routine that is interesting for inverse design is
//! generic over [`Scalar`]. Two implementations are provided:
//!
//! - `f64` — the fast path.
//! - [`crate::dual::Dual`] — forward-mode automatic differentiation, so the
//!   same trace code produces exact first derivatives.

use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A field-like scalar with the transcendental functions optics needs.
///
/// Implemented by `f64` and by [`crate::dual::Dual`] (forward-mode AD).
pub trait Scalar:
    Copy
    + fmt::Debug
    + fmt::Display
    + PartialEq
    + PartialOrd
    + Send
    + Sync
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
{
    /// The additive identity.
    const ZERO: Self;
    /// The multiplicative identity.
    const ONE: Self;

    /// Construct from an `f64` constant.
    fn from_f64(v: f64) -> Self;

    /// Lossy conversion to `f64` (identity for `f64`, value part for
    /// [`Dual`](crate::dual::Dual)).
    fn to_f64(self) -> f64;

    /// Square root.
    fn sqrt(self) -> Self;
    /// Absolute value.
    fn abs(self) -> Self;
    /// Sine.
    fn sin(self) -> Self;
    /// Cosine.
    fn cos(self) -> Self;
    /// Exponential.
    fn exp(self) -> Self;
    /// Natural logarithm.
    fn ln(self) -> Self;
    /// Power with a (locally constant) exponent.
    fn powf(self, e: Self) -> Self;
    /// Power with an integer exponent.
    fn powi(self, n: i32) -> Self;
    /// Four-quadrant arctangent of `self` (y) and `x`.
    fn atan2(self, x: Self) -> Self;
    /// True if the value is neither infinite nor NaN.
    fn is_finite(self) -> bool;

    /// `self * self`, in `f64` terms.
    fn squared(self) -> Self {
        self * self
    }
}

impl Scalar for f64 {
    const ZERO: f64 = 0.0;
    const ONE: f64 = 1.0;

    #[inline]
    fn from_f64(v: f64) -> Self {
        v
    }

    #[inline]
    fn to_f64(self) -> f64 {
        self
    }

    #[inline]
    fn sqrt(self) -> Self {
        libm::sqrt(self)
    }

    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }

    #[inline]
    fn sin(self) -> Self {
        libm::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        libm::cos(self)
    }

    #[inline]
    fn exp(self) -> Self {
        libm::exp(self)
    }

    #[inline]
    fn ln(self) -> Self {
        libm::log(self)
    }

    #[inline]
    fn powf(self, e: Self) -> Self {
        libm::pow(self, e)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        libm::pow(self, n as f64)
    }

    #[inline]
    fn atan2(self, x: Self) -> Self {
        libm::atan2(self, x)
    }

    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }
}

/// Linear interpolation: `a + t * (b - a)`.
pub fn lerp<T: Scalar>(a: T, b: T, t: T) -> T {
    a + t * (b - a)
}

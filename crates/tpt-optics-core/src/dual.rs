//! Forward-mode automatic differentiation via dual numbers.
//!
//! A [`Dual`] number carries a value `v` and a first derivative `d` with
//! respect to the single seeding variable. Arithmetic and the transcendental
//! functions propagate the derivative exactly (to floating-point precision),
//! which is what makes tpt-optics' generic ray tracing differentiable: trace
//! a system with one parameter seeded and the output derivative is the exact
//! Jacobian column, without finite-difference noise.

use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::num::Scalar;

/// A dual number `v + ε·d` with `ε² = 0`.
#[derive(Copy, Clone, Debug)]
pub struct Dual {
    /// Function value.
    pub v: f64,
    /// First derivative with respect to the seed variable.
    pub d: f64,
}

impl PartialEq for Dual {
    /// Equal by value only, consistent with [`PartialOrd`]: a seeded zero
    /// (`v = 0`, `d ≠ 0`) compares equal to a constant zero.
    fn eq(&self, other: &Self) -> bool {
        self.v == other.v
    }
}

impl PartialOrd for Dual {
    /// Ordered by the value part only (derivatives don't participate).
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        self.v.partial_cmp(&other.v)
    }
}

impl Dual {
    /// A dual number with zero derivative (a constant).
    #[inline]
    pub const fn new(v: f64) -> Self {
        Dual { v, d: 0.0 }
    }

    /// A dual number with unit derivative in the seed variable.
    #[inline]
    pub const fn var(v: f64) -> Self {
        Dual { v, d: 1.0 }
    }

    /// A constant combined with an explicit derivative.
    #[inline]
    pub const fn with_derivative(v: f64, d: f64) -> Self {
        Dual { v, d }
    }
}

impl From<f64> for Dual {
    #[inline]
    fn from(v: f64) -> Self {
        Dual::new(v)
    }
}

impl fmt::Display for Dual {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + ε·{}", self.v, self.d)
    }
}

impl Add for Dual {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Dual {
            v: self.v + rhs.v,
            d: self.d + rhs.d,
        }
    }
}

impl Sub for Dual {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Dual {
            v: self.v - rhs.v,
            d: self.d - rhs.d,
        }
    }
}

impl Mul for Dual {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Dual {
            v: self.v * rhs.v,
            d: self.d * rhs.v + self.v * rhs.d,
        }
    }
}

impl Div for Dual {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        let inv = 1.0 / (rhs.v * rhs.v);
        Dual {
            v: self.v / rhs.v,
            d: (self.d * rhs.v - self.v * rhs.d) * inv,
        }
    }
}

impl Neg for Dual {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Dual {
            v: -self.v,
            d: -self.d,
        }
    }
}

impl AddAssign for Dual {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Dual {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for Dual {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl DivAssign for Dual {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Scalar for Dual {
    const ZERO: Dual = Dual { v: 0.0, d: 0.0 };
    const ONE: Dual = Dual { v: 1.0, d: 0.0 };

    #[inline]
    fn from_f64(v: f64) -> Self {
        Dual::new(v)
    }

    #[inline]
    fn to_f64(self) -> f64 {
        self.v
    }

    #[inline]
    fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        Dual {
            v: s,
            d: if s == 0.0 { 0.0 } else { self.d / (2.0 * s) },
        }
    }

    #[inline]
    fn abs(self) -> Self {
        Dual {
            v: self.v.abs(),
            d: self.d * self.v.signum(),
        }
    }

    #[inline]
    fn sin(self) -> Self {
        Dual {
            v: libm::sin(self.v),
            d: self.d * libm::cos(self.v),
        }
    }

    #[inline]
    fn cos(self) -> Self {
        Dual {
            v: libm::cos(self.v),
            d: -self.d * libm::sin(self.v),
        }
    }

    #[inline]
    fn exp(self) -> Self {
        let e = libm::exp(self.v);
        Dual {
            v: e,
            d: self.d * e,
        }
    }

    #[inline]
    fn ln(self) -> Self {
        Dual {
            v: libm::log(self.v),
            d: self.d / self.v,
        }
    }

    /// Power with the exponent treated as a constant. (Differentiating with
    /// respect to the exponent is not used by any tpt-optics routine; the
    /// exponent's derivative component is ignored.)
    #[inline]
    fn powf(self, e: Self) -> Self {
        let p = libm::pow(self.v, e.v);
        Dual {
            v: p,
            d: self.d * e.v * libm::pow(self.v, e.v - 1.0),
        }
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        Dual {
            v: self.v.powi(n),
            d: n as f64 * libm::pow(self.v, n as f64 - 1.0) * self.d,
        }
    }

    /// atan2(y=self, x) derivative: `(x·dy − y·dx) / (x² + y²)`.
    #[inline]
    fn atan2(self, x: Self) -> Self {
        let denom = x.v * x.v + self.v * self.v;
        Dual {
            v: libm::atan2(self.v, x.v),
            d: (x.v * self.d - self.v * x.d) / denom,
        }
    }

    #[inline]
    fn is_finite(self) -> bool {
        self.v.is_finite() && self.d.is_finite()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn arithmetic_chain_rule() {
        let x = Dual::var(2.0);
        // f(x) = x^3 + 2x, f' = 3x^2 + 2
        let f = x * x * x + Dual::new(2.0) * x;
        assert!(close(f.v, 12.0, 1e-12));
        assert!(close(f.d, 14.0, 1e-12));
    }

    #[test]
    fn transcendentals() {
        let x = Dual::var(0.7);
        let f = x.exp() * x.sin();
        assert!(close(f.v, 0.7f64.exp() * 0.7f64.sin(), 1e-14));
        assert!(
            close(f.d, 0.7f64.exp() * (0.7f64.sin() + 0.7f64.cos()), 1e-14),
            "got {}",
            f.d
        );

        let g = (x * x + Dual::new(1.0)).ln();
        assert!(close(g.d, 2.0 * 0.7 / (0.49 + 1.0), 1e-14));
    }

    #[test]
    fn quotient_and_sqrt() {
        let x = Dual::var(3.0);
        let f = x.sqrt() / (x + Dual::new(1.0));
        let s = 3.0f64.sqrt();
        assert!(close(f.v, s / 4.0, 1e-14));
        // d/dx [sqrt(x)/(x+1)] = (1/(2√x)(x+1) − √x)/(x+1)²
        let expect = (4.0 / (2.0 * s) - s) / 16.0;
        assert!(close(f.d, expect, 1e-13));
    }

    #[test]
    fn constants_have_zero_derivative() {
        let x = Dual::new(1.5) + Dual::new(2.0);
        assert_eq!(x.v, 3.5);
        assert_eq!(x.d, 0.0);
        // A constant shifted by the variable keeps the variable's slope.
        let y = Dual::var(1.5) + Dual::new(2.0);
        assert_eq!(y.d, 1.0);
    }
}

#[cfg(test)]
mod equality_tests {
    use super::*;

    #[test]
    fn seeded_zero_equals_constant_zero() {
        // Equality must agree with the value-only ordering; a seeded zero
        // (e.g. a differentiated flat surface's curvature) is still zero.
        let seeded = Dual { v: 0.0, d: 1.0 };
        assert_eq!(seeded, Dual::ZERO);
        assert_eq!(
            seeded.partial_cmp(&Dual::ZERO),
            Some(core::cmp::Ordering::Equal)
        );
    }
}

//! Refractive index representations.

use core::fmt;

use crate::units::{Kelvin, Wavelength};

/// A real (non-absorbing) refractive index.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct RefractiveIndex(f64);

impl RefractiveIndex {
    /// Construct. Panics in debug builds for non-positive/non-finite values;
    /// silently stores them otherwise (callers check).
    #[inline]
    pub fn new(n: f64) -> Self {
        debug_assert!(
            n.is_finite() && n > 0.0,
            "refractive index must be positive"
        );
        RefractiveIndex(n)
    }

    /// The index value.
    #[inline]
    pub fn value(self) -> f64 {
        self.0
    }

    /// Relative permittivity ε = n² (non-magnetic medium).
    #[inline]
    pub fn permittivity(self) -> f64 {
        self.0 * self.0
    }

    /// Vacuum = 1.
    #[inline]
    pub fn vacuum() -> Self {
        RefractiveIndex(1.0)
    }
}

impl fmt::Display for RefractiveIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "n = {}", self.0)
    }
}

/// Complex refractive index `n + i·k`.
///
/// The imaginary part `k` (the *extinction coefficient*) is non-negative for
/// passive absorbing media; `k = 0` means lossless.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct ComplexRefractiveIndex {
    /// Real part (phase index).
    pub n: f64,
    /// Imaginary part (extinction coefficient, ≥ 0 for passive media).
    pub k: f64,
}

impl ComplexRefractiveIndex {
    /// Construct from `n` and `k`.
    #[inline]
    pub const fn new(n: f64, k: f64) -> Self {
        ComplexRefractiveIndex { n, k }
    }

    /// A lossless real index.
    #[inline]
    pub fn real(n: f64) -> Self {
        ComplexRefractiveIndex { n, k: 0.0 }
    }

    /// The real part.
    #[inline]
    pub fn re(self) -> f64 {
        self.n
    }

    /// The imaginary part.
    #[inline]
    pub fn im(self) -> f64 {
        self.k
    }

    /// True when the medium absorbs (`k > 0`).
    #[inline]
    pub fn is_absorbing(self) -> bool {
        self.k > 0.0
    }

    /// Relative complex permittivity `ε = (n + i·k)²`.
    #[inline]
    pub fn permittivity(self) -> crate::complex::Complex<f64> {
        let z = crate::complex::Complex::new(self.n, self.k);
        z * z
    }

    /// Intensity absorption coefficient `α = 4π·k / λ0` in 1/m
    /// (Beer–Lambert `I = I₀·exp(−α·d)`), for vacuum wavelength `λ0`.
    #[inline]
    pub fn absorption_coefficient(self, wl: Wavelength) -> f64 {
        4.0 * core::f64::consts::PI * self.k / wl.m()
    }

    /// Power reflectance at normal incidence from air/vacuum.
    #[inline]
    pub fn reflectance_normal(self) -> f64 {
        // |(ñ − 1)/(ñ + 1)|² with ñ = n + ik.
        let num = (self.n - 1.0) * (self.n - 1.0) + self.k * self.k;
        let den = (self.n + 1.0) * (self.n + 1.0) + self.k * self.k;
        num / den
    }
}

impl From<RefractiveIndex> for ComplexRefractiveIndex {
    #[inline]
    fn from(n: RefractiveIndex) -> Self {
        ComplexRefractiveIndex::real(n.value())
    }
}

impl From<f64> for ComplexRefractiveIndex {
    #[inline]
    fn from(n: f64) -> Self {
        ComplexRefractiveIndex::real(n)
    }
}

/// A medium that provides a (possibly wavelength/temperature dependent)
/// refractive index.
///
/// This is the interface the ray tracers consume; it keeps geometry crates
/// decoupled from specific dispersion models.
pub trait IndexAt: fmt::Debug + Send + Sync {
    /// Complex refractive index at the given vacuum wavelength and
    /// temperature.
    fn index_at(&self, wavelength: Wavelength, temperature: Kelvin) -> ComplexRefractiveIndex;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::constants;

    #[test]
    fn complex_index_relations() {
        let idx = ComplexRefractiveIndex::new(1.5, 0.1);
        let eps = idx.permittivity();
        // (1.5 + 0.1i)^2 = 2.24 + 0.3i
        assert!((eps.re - 2.24).abs() < 1e-12);
        assert!((eps.im - 0.3).abs() < 1e-12);

        // α = 4πk/λ
        let wl = Wavelength::from_nm(500.0);
        assert!(
            (idx.absorption_coefficient(wl) - 4.0 * core::f64::consts::PI * 0.1 / 5e-7).abs()
                < 1e-9
        );

        // Normal-incidence reflectance of glass: ((1.5-1)/(1.5+1))^2 = 0.04
        assert!((ComplexRefractiveIndex::real(1.5).reflectance_normal() - 0.04).abs() < 1e-15);
        let _ = constants::SPEED_OF_LIGHT;
    }
}

#[cfg(test)]
mod reflectance_tests {
    use super::*;

    #[test]
    fn metal_reflectance_includes_extinction() {
        // ñ = 0.2 + 3.0i: R = ((n−1)² + k²)/((n+1)² + k²) = 9.64/10.44.
        let r = ComplexRefractiveIndex::new(0.2, 3.0).reflectance_normal();
        assert!((r - 9.64 / 10.44).abs() < 1e-12, "R = {r}");
    }
}

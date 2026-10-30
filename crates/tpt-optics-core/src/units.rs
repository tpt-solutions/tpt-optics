//! Unit-safe optical physical quantities.
//!
//! Quantities are stored in SI base units and converted at the boundary.
//! Every newtype is `Copy` and cheap; arithmetic is restricted to
//! dimensionally meaningful operations (scaling, comparison), preventing
//! classic unit errors like adding a length to a wavelength or treating
//! degrees as radians.

use core::fmt;
use core::ops::Mul;

/// Physical constants used by the suite (CODATA-ish SI values).
pub mod constants {
    /// Speed of light in vacuum, m/s.
    pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;
    /// Planck constant, J·s.
    pub const PLANCK: f64 = 6.626_070_15e-34;
    /// Reduced Planck constant, J·s.
    pub const HBAR: f64 = 1.054_571_817e-34;
    /// Elementary charge, C.
    pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;
    /// Boltzmann constant, J/K.
    pub const BOLTZMANN: f64 = 1.380_649e-23;
    /// Vacuum permeability, H/m.
    pub const VACUUM_PERMEABILITY: f64 = 1.256_637_062_12e-6;
    /// Vacuum permittivity, F/m.
    pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_812_8e-12;
    /// Impedance of free space, Ω.
    pub const VACUUM_IMPEDANCE: f64 = 376.730_313_668;
}

macro_rules! unit_type {
    ($(#[$meta:meta])* $name:ident, $unit:expr) => {
        $(#[$meta])*
        #[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd)]
        pub struct $name(f64);

        impl $name {
            /// Raw value in SI base units.
            #[inline]
            pub const fn value(self) -> f64 {
                self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{} {}", self.0, $unit)
            }
        }

        impl Mul<f64> for $name {
            type Output = Self;
            #[inline]
            fn mul(self, s: f64) -> Self {
                Self(self.0 * s)
            }
        }
    };
}

unit_type!(
    /// Optical wavelength in vacuum, stored in meters.
    Wavelength,
    "m"
);
unit_type!(
    /// Temporal frequency, stored in hertz.
    Frequency,
    "Hz"
);
unit_type!(
    /// Optical path length (physical length × refractive index), meters.
    OpticalPathLength,
    "m"
);
unit_type!(
    /// Thermodynamic temperature, kelvin.
    Kelvin,
    "K"
);
unit_type!(
    /// Plane angle, radians.
    Angle,
    "rad"
);

impl Wavelength {
    /// Construct from nanometers.
    #[inline]
    pub fn from_nm(nm: f64) -> Self {
        Wavelength(nm * 1e-9)
    }

    /// Construct from micrometers.
    #[inline]
    pub fn from_um(um: f64) -> Self {
        Wavelength(um * 1e-6)
    }

    /// Construct from meters.
    #[inline]
    pub fn from_m(m: f64) -> Self {
        Wavelength(m)
    }

    /// Value in meters.
    #[inline]
    pub fn m(self) -> f64 {
        self.0
    }

    /// Value in micrometers.
    #[inline]
    pub fn um(self) -> f64 {
        self.0 * 1e6
    }

    /// Value in nanometers.
    #[inline]
    pub fn nm(self) -> f64 {
        self.0 * 1e9
    }

    /// Corresponding temporal frequency in vacuum.
    #[inline]
    pub fn frequency(self) -> Frequency {
        Frequency(constants::SPEED_OF_LIGHT / self.0)
    }

    /// Corresponding vacuum wavenumber in rad/m (`2π/λ`).
    #[inline]
    pub fn wavenumber(self) -> f64 {
        core::f64::consts::TAU / self.0
    }

    /// Photon energy in electronvolts.
    #[inline]
    pub fn energy_ev(self) -> f64 {
        constants::PLANCK * constants::SPEED_OF_LIGHT / self.0 / constants::ELEMENTARY_CHARGE
    }

    /// Wavelength inside a medium of refractive index `n` (λ₀/n).
    #[inline]
    pub fn in_medium(self, n: f64) -> Wavelength {
        Wavelength(self.0 / n)
    }
}

impl Frequency {
    /// Construct from hertz.
    #[inline]
    pub fn from_hz(hz: f64) -> Self {
        Frequency(hz)
    }

    /// Construct from terahertz.
    #[inline]
    pub fn from_thz(thz: f64) -> Self {
        Frequency(thz * 1e12)
    }

    /// Value in hertz.
    #[inline]
    pub fn hz(self) -> f64 {
        self.0
    }

    /// Value in terahertz.
    #[inline]
    pub fn thz(self) -> f64 {
        self.0 * 1e-12
    }

    /// Corresponding vacuum wavelength.
    #[inline]
    pub fn wavelength(self) -> Wavelength {
        Wavelength(constants::SPEED_OF_LIGHT / self.0)
    }
}

impl OpticalPathLength {
    /// Construct from meters of optical path.
    #[inline]
    pub fn from_m(m: f64) -> Self {
        OpticalPathLength(m)
    }

    /// Construct from millimeters of optical path.
    #[inline]
    pub fn from_mm(mm: f64) -> Self {
        OpticalPathLength(mm * 1e-3)
    }

    /// Value in meters of optical path.
    #[inline]
    pub fn m(self) -> f64 {
        self.0
    }

    /// Optical path expressed in waves at the given wavelength.
    #[inline]
    pub fn waves(self, wl: Wavelength) -> f64 {
        self.0 / wl.m()
    }
}

impl Kelvin {
    /// Construct from kelvin.
    #[inline]
    pub fn from_k(k: f64) -> Self {
        Kelvin(k)
    }

    /// Construct from degrees Celsius.
    #[inline]
    pub fn from_celsius(c: f64) -> Self {
        Kelvin(c + 273.15)
    }

    /// Value in kelvin.
    #[inline]
    pub fn k(self) -> f64 {
        self.0
    }

    /// Value in degrees Celsius.
    #[inline]
    pub fn celsius(self) -> f64 {
        self.0 - 273.15
    }
}

/// Reference temperature assumed by most catalog data (20 °C).
pub const REFERENCE_TEMPERATURE: Kelvin = Kelvin(293.15);

impl Angle {
    /// Construct from radians.
    #[inline]
    pub fn from_rad(rad: f64) -> Self {
        Angle(rad)
    }

    /// Construct from degrees.
    #[inline]
    pub fn from_degrees(deg: f64) -> Self {
        Angle(deg * core::f64::consts::PI / 180.0)
    }

    /// Value in radians.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
    }

    /// Value in degrees.
    #[inline]
    pub fn degrees(self) -> f64 {
        self.0 * 180.0 / core::f64::consts::PI
    }

    /// Sine of the angle.
    #[inline]
    pub fn sin(self) -> f64 {
        libm::sin(self.0)
    }

    /// Cosine of the angle.
    #[inline]
    pub fn cos(self) -> f64 {
        libm::cos(self.0)
    }

    /// Tangent of the angle.
    #[inline]
    pub fn tan(self) -> f64 {
        libm::tan(self.0)
    }

    /// `arcsin`, as an [`Angle`].
    #[inline]
    pub fn asin(x: f64) -> Option<Self> {
        if (-1.0..=1.0).contains(&x) {
            Some(Angle(libm::asin(x.clamp(-1.0, 1.0))))
        } else {
            None
        }
    }

    /// `arccos`, as an [`Angle`].
    #[inline]
    pub fn acos(x: f64) -> Option<Self> {
        if (-1.0..=1.0).contains(&x) {
            Some(Angle(libm::acos(x.clamp(-1.0, 1.0))))
        } else {
            None
        }
    }

    /// `arctan`, as an [`Angle`].
    #[inline]
    pub fn atan(x: f64) -> Self {
        Angle(libm::atan(x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wavelength_conversions() {
        let wl = Wavelength::from_nm(632.8);
        assert!((wl.um() - 0.6328).abs() < 1e-12);
        assert!((wl.m() - 6.328e-7).abs() < 1e-20);
        // f = c/λ ≈ 473.7 THz
        assert!((wl.frequency().thz() - 473.74).abs() < 0.05);
        // Photon energy ≈ 1.959 eV
        assert!((wl.energy_ev() - 1.9593).abs() < 1e-3);
    }

    #[test]
    fn path_in_waves() {
        let opl = OpticalPathLength::from_mm(1.0);
        let waves = opl.waves(Wavelength::from_nm(1000.0));
        assert!((waves - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn angles() {
        let a = Angle::from_degrees(90.0);
        assert!((a.rad() - core::f64::consts::FRAC_PI_2).abs() < 1e-15);
        assert!((Angle::from_degrees(30.0).sin() - 0.5).abs() < 1e-15);
        assert!(Angle::asin(1.5).is_none());
        assert!((Angle::asin(0.5).unwrap().degrees() - 30.0).abs() < 1e-12);
    }

    #[test]
    fn celsius_roundtrip() {
        let t = Kelvin::from_celsius(20.0);
        assert!((t.k() - 293.15).abs() < 1e-12);
        assert!((t.celsius() - 20.0).abs() < 1e-12);
    }

    #[test]
    fn in_medium() {
        let wl = Wavelength::from_nm(600.0);
        assert!((wl.in_medium(1.5).nm() - 400.0).abs() < 1e-9);
    }
}

/// Geometric length, stored in meters. `Display` picks a scale that keeps
/// the number readable: metres, millimetres, micrometres or nanometres.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Length(f64);

impl Length {
    /// From metres.
    #[inline]
    pub const fn from_m(m: f64) -> Self {
        Length(m)
    }

    /// From millimetres.
    #[inline]
    pub fn from_mm(mm: f64) -> Self {
        Length(mm * 1e-3)
    }

    /// From micrometres.
    #[inline]
    pub fn from_um(um: f64) -> Self {
        Length(um * 1e-6)
    }

    /// From nanometres.
    #[inline]
    pub fn from_nm(nm: f64) -> Self {
        Length(nm * 1e-9)
    }

    /// In metres.
    #[inline]
    pub const fn m(self) -> f64 {
        self.0
    }

    /// In millimetres.
    #[inline]
    pub fn mm(self) -> f64 {
        self.0 * 1e3
    }

    /// In micrometres.
    #[inline]
    pub fn um(self) -> f64 {
        self.0 * 1e6
    }

    /// In nanometres.
    #[inline]
    pub fn nm(self) -> f64 {
        self.0 * 1e9
    }
}

impl fmt::Display for Length {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let a = libm::fabs(self.0);
        if a == 0.0 || a >= 1.0 {
            write!(f, "{} m", self.0)
        } else if a >= 1e-3 {
            write!(f, "{} mm", self.mm())
        } else if a >= 1e-6 {
            write!(f, "{} µm", self.um())
        } else {
            write!(f, "{} nm", self.nm())
        }
    }
}

/// Surface curvature `c = 1/R`, stored in inverse metres. A radius of zero
/// is a flat surface, which has zero curvature.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Curvature(f64);

impl Curvature {
    /// From a curvature in inverse metres.
    #[inline]
    pub const fn from_inverse_m(c: f64) -> Self {
        Curvature(c)
    }

    /// From a radius of curvature in millimetres. Zero gives a flat surface.
    #[inline]
    pub fn from_radius_mm(r: f64) -> Self {
        if r == 0.0 {
            Curvature(0.0)
        } else {
            Curvature(1e3 / r)
        }
    }

    /// In inverse metres.
    #[inline]
    pub const fn inverse_m(self) -> f64 {
        self.0
    }

    /// Radius of curvature in millimetres, or `None` for a flat surface.
    #[inline]
    pub fn radius_mm(self) -> Option<f64> {
        if self.0 == 0.0 {
            None
        } else {
            Some(1e3 / self.0)
        }
    }
}

impl fmt::Display for Curvature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.radius_mm() {
            Some(r) => write!(f, "R = {r} mm"),
            None => write!(f, "flat"),
        }
    }
}

/// Optical power in dioptres (inverse metres), `P = 1/f`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Power(f64);

impl Power {
    /// From dioptres (1/m).
    #[inline]
    pub const fn from_dioptres(d: f64) -> Self {
        Power(d)
    }

    /// From a focal length in millimetres. Infinite focal length is zero power.
    #[inline]
    pub fn from_focal_length_mm(f: f64) -> Self {
        if f.is_infinite() {
            Power(0.0)
        } else {
            Power(1e3 / f)
        }
    }

    /// In dioptres.
    #[inline]
    pub const fn dioptres(self) -> f64 {
        self.0
    }

    /// Focal length in millimetres, or `None` for zero power.
    #[inline]
    pub fn focal_length_mm(self) -> Option<f64> {
        if self.0 == 0.0 {
            None
        } else {
            Some(1e3 / self.0)
        }
    }
}

impl fmt::Display for Power {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} D", self.0)
    }
}

#[cfg(test)]
mod quantity_tests {
    use super::*;

    #[test]
    fn length_converts_between_scales() {
        assert!((Length::from_mm(12.5).m() - 12.5e-3).abs() < 1e-15);
        assert!((Length::from_um(0.5).nm() - 500.0).abs() < 1e-9);
    }

    #[test]
    fn length_display_picks_a_readable_scale() {
        assert_eq!(alloc::format!("{}", Length::from_mm(12.5)), "12.5 mm");
        assert_eq!(alloc::format!("{}", Length::from_nm(632.8)), "632.8 nm");
        assert_eq!(alloc::format!("{}", Length::from_um(2.0)), "2 µm");
        assert_eq!(alloc::format!("{}", Length::from_m(1.5)), "1.5 m");
    }

    #[test]
    fn curvature_round_trips_through_radius() {
        let c = Curvature::from_radius_mm(100.0);
        assert!((c.inverse_m() - 10.0).abs() < 1e-12);
        assert!((c.radius_mm().unwrap() - 100.0).abs() < 1e-12);
        assert_eq!(Curvature::from_radius_mm(0.0).radius_mm(), None);
    }

    #[test]
    fn power_is_reciprocal_focal_length() {
        let p = Power::from_focal_length_mm(100.0);
        assert!((p.dioptres() - 10.0).abs() < 1e-12);
        assert_eq!(
            Power::from_focal_length_mm(f64::INFINITY).focal_length_mm(),
            None
        );
    }
}

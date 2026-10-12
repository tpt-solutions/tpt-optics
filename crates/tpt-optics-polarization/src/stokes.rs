//! Stokes vectors: the real 4-vector description of (partially) polarized
//! light: `S = (I, Q, U, V)` with `Q = I_x − I_y`, `U = I_45 − I_135`.
//!
//! Handedness convention (used throughout this crate): `V = −2·Im(E_x·E_y*)`,
//! so the right-circular state `E_y = −i·E_x` ([`JonesVector::circular_right`])
//! has `V = −1`. Under this convention `V = I_LCP − I_RCP`.

use crate::jones::JonesVector;

/// A Stokes vector `(s0, s1, s2, s3)` (intensity units).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct StokesVector {
    /// Total intensity.
    pub s0: f64,
    /// Linear 0°/90° component.
    pub s1: f64,
    /// Linear ±45° component.
    pub s2: f64,
    /// Circular component.
    pub s3: f64,
}

impl StokesVector {
    /// Construct from components.
    pub const fn new(s0: f64, s1: f64, s2: f64, s3: f64) -> Self {
        StokesVector { s0, s1, s2, s3 }
    }

    /// Fully unpolarized light of intensity `i`.
    pub fn unpolarized(i: f64) -> Self {
        StokesVector::new(i, 0.0, 0.0, 0.0)
    }

    /// Fully linearly polarized light at azimuth `psi` (rad), intensity `i`.
    pub fn linear(psi: f64, i: f64) -> Self {
        StokesVector::new(i, i * libm::cos(2.0 * psi), i * libm::sin(2.0 * psi), 0.0)
    }

    /// Fully circular polarization: `sign > 0` right (`E_y = −i·E_x`), `sign < 0`
    /// left. Matches [`StokesVector::from_jones`] of the same Jones state.
    pub fn circular(sign: f64, i: f64) -> Self {
        StokesVector::new(i, 0.0, 0.0, -i * sign.signum())
    }

    /// Jones vector → Stokes (deterministic, fully polarized).
    #[must_use]
    pub fn from_jones(j: &JonesVector) -> Self {
        let ix = j.ex.norm_sqr();
        let iy = j.ey.norm_sqr();
        let exy = j.ex * j.ey.conj();
        StokesVector::new(ix + iy, ix - iy, 2.0 * exy.re, -2.0 * exy.im)
    }

    /// Degree of polarization `sqrt(Q²+U²+V²)/I` in [0, 1].
    #[must_use]
    pub fn degree_of_polarization(&self) -> f64 {
        if self.s0 <= 0.0 {
            return 0.0;
        }
        (libm::sqrt(self.s1 * self.s1 + self.s2 * self.s2 + self.s3 * self.s3) / self.s0)
            .clamp(0.0, 1.0 + 1e-12)
    }

    /// Normalized form `(1, Q/I, U/I, V/I)`.
    #[must_use]
    pub fn normalized(&self) -> StokesVector {
        if self.s0 <= 0.0 {
            return *self;
        }
        StokesVector::new(1.0, self.s1 / self.s0, self.s2 / self.s0, self.s3 / self.s0)
    }

    /// Azimuth of the polarization ellipse (rad), if polarized.
    #[must_use]
    pub fn azimuth(&self) -> Option<f64> {
        if self.s1 == 0.0 && self.s2 == 0.0 {
            None
        } else {
            Some(0.5 * libm::atan2(self.s2, self.s1))
        }
    }

    /// Ellipticity angle (rad): `χ = asin(V / sqrt(Q²+U²+V²))`.
    #[must_use]
    pub fn ellipticity(&self) -> Option<f64> {
        let p = libm::sqrt(self.s1 * self.s1 + self.s2 * self.s2 + self.s3 * self.s3);
        if p <= 0.0 {
            None
        } else {
            Some(libm::asin(self.s3 / p))
        }
    }

    /// Incoherent sum of two beams.
    #[must_use]
    pub fn add(&self, other: &StokesVector) -> StokesVector {
        StokesVector::new(
            self.s0 + other.s0,
            self.s1 + other.s1,
            self.s2 + other.s2,
            self.s3 + other.s3,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn jones_to_stokes_basic_states() {
        let sx = StokesVector::from_jones(&JonesVector::linear_x());
        assert!(close(sx.s0, 1.0, 1e-15));
        assert!(close(sx.s1, 1.0, 1e-15));
        assert!(close(sx.s2, 0.0, 1e-15));
        assert!(close(sx.s3, 0.0, 1e-15));

        let s45 = StokesVector::from_jones(&JonesVector::linear_at(core::f64::consts::FRAC_PI_4));
        assert!(close(s45.s1, 0.0, 1e-15));
        assert!(close(s45.s2, 1.0, 1e-15));

        let sr = StokesVector::from_jones(&JonesVector::circular_right());
        assert!(close(sr.s3, -1.0, 1e-15), "s3 = {}", sr.s3);
        assert!(close(sr.s1, 0.0, 1e-15));
    }

    #[test]
    fn stokes_constructors_consistent() {
        let sl = StokesVector::linear(0.3, 5.0);
        assert!(close(sl.degree_of_polarization(), 1.0, 1e-12));
        assert!(close(sl.azimuth().unwrap(), 0.3, 1e-12));
        // Right-circular: ellipticity is −π/2 under V = −1 for RCP.
        let sc = StokesVector::circular(1.0, 2.0);
        assert!(close(
            sc.ellipticity().unwrap(),
            -core::f64::consts::FRAC_PI_2,
            1e-12
        ));
        let sj = StokesVector::from_jones(&JonesVector::circular_right());
        assert!(close(sc.s3, sj.s3 * 2.0, 1e-12), "{} vs {}", sc.s3, sj.s3);
        let su = StokesVector::unpolarized(3.0);
        assert!(close(su.degree_of_polarization(), 0.0, 1e-15));
        assert!(su.azimuth().is_none());
    }

    #[test]
    fn mixed_beam_dop() {
        // Half x-polarized + half unpolarized: DOP = 0.5.
        let a = StokesVector::linear(0.0, 1.0);
        let b = StokesVector::unpolarized(1.0);
        let m = a.add(&b);
        assert!(close(m.degree_of_polarization(), 0.5, 1e-12));
        assert!(close(m.normalized().s1, 0.5, 1e-12));
    }
}

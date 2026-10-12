//! Gaussian beams via the complex beam parameter (q).
//!
//! `q` encodes the wavefront curvature and spot size:
//! `1/q = 1/R − i·λ/(π·w²)`. Free space adds distance; ABCD optics transform
//! `q → (Aq+B)/(Cq+D)` (the ABCD law), implemented through the geometric
//! crate's [`Matrix2`](tpt_optics_core::Matrix2).

use tpt_optics_core::{Complex, Matrix2, Wavelength};

/// A Gaussian beam at a plane, described by its complex parameter.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GaussianBeam {
    /// Complex beam parameter q (meters).
    pub q: Complex<f64>,
    /// Vacuum wavelength (meters).
    pub wavelength: Wavelength,
    /// Refractive index of the host medium.
    pub index: f64,
}

impl GaussianBeam {
    /// A beam with waist `w0` at the current plane.
    #[must_use]
    pub fn at_waist(w0: f64, wavelength: Wavelength, index: f64) -> Self {
        let zr = rayleigh_range(w0, wavelength, index);
        GaussianBeam {
            q: Complex::new(0.0, zr),
            wavelength,
            index,
        }
    }

    /// A beam with given spot size `w` and wavefront radius `R` at the
    /// current plane (`R` may be infinite/`None` for a waist).
    #[must_use]
    pub fn from_w_and_r(w: f64, r: Option<f64>, wavelength: Wavelength, index: f64) -> Self {
        let inv_q_re = r.map_or(0.0, |r| 1.0 / r);
        let inv_q_im = -(wavelength.m() / index) / (core::f64::consts::PI * w * w);
        // 1/q = a + ib → q = 1/(a+ib)
        let inv_q = Complex::new(inv_q_re, inv_q_im);
        GaussianBeam {
            q: inv_q.inv(),
            wavelength,
            index,
        }
    }

    /// Rayleigh range z_R = π n w₀²/λ (waist-parameterised).
    #[must_use]
    pub fn rayleigh_range(&self) -> f64 {
        self.q.im
    }

    /// Spot size at this plane.
    #[must_use]
    pub fn w(&self) -> f64 {
        let inv = self.q.inv();
        libm::sqrt((self.wavelength.m() / self.index) / ((-inv.im).abs() * core::f64::consts::PI))
    }

    /// Wavefront curvature radius at this plane (None at a waist).
    #[must_use]
    pub fn r(&self) -> Option<f64> {
        let inv = self.q.inv();
        if inv.re.abs() < 1e-15 {
            None
        } else {
            Some(1.0 / inv.re)
        }
    }

    /// Propagate over free-space distance `d`.
    #[must_use]
    pub fn propagated(&self, d: f64) -> Self {
        self.through(&Matrix2::new(1.0, d, 0.0, 1.0))
    }

    /// Apply an arbitrary ABCD system (in the beam's index convention).
    #[must_use]
    pub fn through(&self, m: &Matrix2<f64>) -> Self {
        GaussianBeam {
            q: m.apply_q(self.q),
            wavelength: self.wavelength,
            index: self.index,
        }
    }

    /// Gouy phase at this plane relative to the nearest waist.
    #[must_use]
    pub fn gouy_phase(&self) -> f64 {
        libm::atan(self.q.re / self.q.im)
    }

    /// Divergence half-angle (rad): θ = atan(λ/πw₀).
    #[must_use]
    pub fn divergence(&self) -> f64 {
        let w0 = self.w0_equivalent();
        libm::atan(self.wavelength.m() / self.index / (core::f64::consts::PI * w0))
    }

    fn w0_equivalent(&self) -> f64 {
        let zr = self.rayleigh_range();
        // z_R = π·n·w0²/λ₀  ⇒  w0 = √(z_R·λ₀/(π·n)).
        libm::sqrt(zr * self.wavelength.m() / (core::f64::consts::PI * self.index))
    }
}

/// Rayleigh range for a waist `w0`: `z_R = π·n·w₀²/λ₀`.
#[must_use]
pub fn rayleigh_range(w0: f64, wavelength: Wavelength, index: f64) -> f64 {
    core::f64::consts::PI * index * w0 * w0 / wavelength.m()
}

/// Find a thin-lens focal length that images an input beam (given at the
/// lens plane) into a target waist `w0_out` located `d_out` after the lens.
///
/// One-dimensional root search over `f` on a logarithmic bracket; robust
/// for the usual mode-matching situation. Returns `None` if no focal length
/// in the scanned range produces the target waist.
#[must_use]
pub fn mode_match_thin_lens(
    q_in: Complex<f64>,
    d_out: f64,
    w0_out: f64,
    wavelength: Wavelength,
    index: f64,
) -> Option<f64> {
    const N_SCAN: usize = 400;
    let objective = |f: f64| -> f64 {
        let lens = Matrix2::new(1.0, 0.0, -1.0 / f, 1.0);
        let space = Matrix2::new(1.0, d_out, 0.0, 1.0);
        let q_out = space.apply_q(lens.apply_q(q_in));
        let beam = GaussianBeam {
            q: q_out,
            wavelength,
            index,
        };
        beam.w() - w0_out
    };
    // Dense logarithmic scan for a sign change, then bisection.
    let (f_min, f_max) = (1e-5_f64, 5.0_f64);
    let log_min = libm::log(f_min);
    let log_step = (libm::log(f_max) - log_min) / (N_SCAN - 1) as f64;
    let mut f_prev = f_min;
    let mut v_prev = objective(f_prev);
    for k in 1..N_SCAN {
        let f = libm::exp(log_min + log_step * k as f64);
        let v = objective(f);
        if v.is_finite() && v_prev.is_finite() && v.signum() != v_prev.signum() {
            let (mut a, mut b) = (f_prev, f);
            for _ in 0..80 {
                let mid = 0.5 * (a + b);
                if objective(mid).signum() == objective(a).signum() {
                    a = mid;
                } else {
                    b = mid;
                }
            }
            return Some(0.5 * (a + b));
        }
        f_prev = f;
        v_prev = v;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn waist_roundtrip() {
        let wl = Wavelength::from_nm(1064.0);
        let w0 = 200e-6;
        let beam = GaussianBeam::at_waist(w0, wl, 1.0);
        assert!(close(beam.w(), w0, 1e-12));
        assert!(beam.r().is_none());
        // From w and R=None again:
        let b2 = GaussianBeam::from_w_and_r(w0, None, wl, 1.0);
        assert!(b2.q.approx_eq(beam.q, 1e-12));
    }

    #[test]
    fn free_space_evolution_matches_analytic() {
        let wl = Wavelength::from_nm(632.8);
        let w0 = 100e-6;
        let zr = rayleigh_range(w0, wl, 1.0);
        let beam = GaussianBeam::at_waist(w0, wl, 1.0);
        for frac in [0.5, 1.0, 3.0] {
            let z = frac * zr;
            let b = beam.propagated(z);
            let w_analytic = w0 * (1.0 + frac * frac).sqrt();
            assert!(close(b.w(), w_analytic, 1e-9), "w at {frac} z_R");
            let r_analytic = z * (1.0 + 1.0 / (frac * frac));
            assert!(close(b.r().unwrap(), r_analytic, 1e-6), "R at {frac} z_R");
        }
        // Round trip: forward z_r then back z_r returns the waist.
        let b = beam.propagated(zr).propagated(-zr);
        assert!(close(b.w(), w0, 1e-9));
        assert!(b.q.approx_eq(beam.q, 1e-9));
    }

    #[test]
    fn thin_lens_abcd_law_focuses() {
        // Waist → thin lens at distance d: focus where the outgoing beam is
        // smallest; check the ABCD-law against the thin-lens formula for a
        // collimated beam: w_f ≈ λ f/(π w0).
        let wl = Wavelength::from_nm(1064.0);
        let w0 = 500e-6;
        let f = 0.05;
        // Collimated: waist far away (q = i z_r with z_r >> f).
        let beam = GaussianBeam::at_waist(w0, wl, 1.0);
        let big = beam.propagated(-1e3); // pull the waist far back
        let lens = Matrix2::new(1.0, 0.0, -1.0 / f, 1.0);
        let after = big.through(&lens);
        // Scan for the focus.
        let mut z_focus = 0.0;
        let mut w_min = f64::INFINITY;
        for i in 0..2000 {
            let z = f * i as f64 / 1000.0;
            let w = after.propagated(z).w();
            if w < w_min {
                w_min = w;
                z_focus = z;
            }
        }
        assert!(
            (z_focus - f).abs() < 2e-3,
            "focus at {z_focus}, expected ≈ {f}"
        );
        let w_expect = wl.m() * f / (core::f64::consts::PI * w0);
        assert!(
            (w_min - w_expect).abs() < 0.05 * w_expect,
            "w_f {w_min} vs {w_expect}"
        );
    }

    #[test]
    fn mode_match_recovers_known_lens() {
        // Image a waist to a target waist; the finder should reproduce the
        // analytic solution for symmetric 4f-like relay: |f| such that
        // magnification = 1 → f = d/4? Verify via self-consistency instead:
        // the found lens must actually produce the target waist.
        let wl = Wavelength::from_nm(1064.0);
        let q_in = GaussianBeam::at_waist(300e-6, wl, 1.0).propagated(0.1).q;
        let d_out = 0.2;
        let w_target = 250e-6; // min achievable at this relay is ~211 µm
        let f = mode_match_thin_lens(q_in, d_out, w_target, wl, 1.0).expect("solution exists");
        let lens = Matrix2::new(1.0, 0.0, -1.0 / f, 1.0);
        let q_out = Matrix2::<f64>::new(1.0, d_out, 0.0, 1.0).apply_q(lens.apply_q(q_in));
        let b = GaussianBeam {
            q: q_out,
            wavelength: wl,
            index: 1.0,
        };
        assert!(
            (b.w() - w_target).abs() < 1e-6,
            "matched w = {} vs {}",
            b.w(),
            w_target
        );
    }

    #[test]
    fn gouy_phase() {
        let wl = Wavelength::from_nm(632.8);
        let beam = GaussianBeam::at_waist(100e-6, wl, 1.0);
        assert!(beam.gouy_phase().abs() < 1e-12);
        let b = beam.propagated(beam.rayleigh_range());
        assert!((b.gouy_phase() - core::f64::consts::FRAC_PI_4).abs() < 1e-6);
    }

    #[test]
    fn divergence_is_lambertian_free() {
        // Divergence θ ≈ λ/(π w0) for a Gaussian.
        let wl = Wavelength::from_nm(1064.0);
        let w0 = 500e-6;
        let beam = GaussianBeam::at_waist(w0, wl, 1.0);
        let theta = wl.m() / (core::f64::consts::PI * w0);
        assert!((beam.divergence() - theta).abs() < 1e-6);
    }
}

#[cfg(test)]
mod divergence_in_medium_tests {
    use super::*;

    #[test]
    fn divergence_matches_waist_in_medium() {
        // θ = λ/(π·n·w0): the wavelength in the medium over π times the waist.
        let (w0, n) = (1.0e-3, 1.5);
        let wl = Wavelength::from_nm(632.8);
        let beam = GaussianBeam::at_waist(w0, wl, n);
        let expect = libm::atan(wl.m() / (core::f64::consts::PI * n * w0));
        assert!(
            (beam.divergence() - expect).abs() < 1e-12,
            "divergence {} vs {expect}",
            beam.divergence()
        );
    }
}

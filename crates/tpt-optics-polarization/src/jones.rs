//! Jones calculus: coherent, fully polarized light as 2-vectors of complex
//! amplitudes manipulated by 2×2 complex matrices.
//!
//! Basis convention: `(ê_x, ê_y)` in the local transverse plane, propagation
//! toward +z. Angles are measured counterclockwise from +x looking into the
//! beam.

use tpt_optics_core::Complex;

/// A Jones vector: transverse complex field amplitudes `(E_x, E_y)`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct JonesVector {
    /// x component.
    pub ex: Complex<f64>,
    /// y component.
    pub ey: Complex<f64>,
}

impl JonesVector {
    /// Construct from components.
    pub const fn new(ex: Complex<f64>, ey: Complex<f64>) -> Self {
        JonesVector { ex, ey }
    }

    /// x-linear unit polarization, unit intensity.
    pub fn linear_x() -> Self {
        JonesVector::new(Complex::one(), Complex::zero())
    }

    /// y-linear unit polarization.
    pub fn linear_y() -> Self {
        JonesVector::new(Complex::zero(), Complex::one())
    }

    /// Linear polarization at `angle` (rad from +x), normalized.
    pub fn linear_at(angle: f64) -> Self {
        JonesVector::new(
            Complex::real(libm::cos(angle)),
            Complex::real(libm::sin(angle)),
        )
    }

    /// Right-circular polarization (E_y = −i·E_x, IEEE/SRP convention).
    pub fn circular_right() -> Self {
        JonesVector::new(
            Complex::real(core::f64::consts::FRAC_1_SQRT_2),
            Complex::new(0.0, -core::f64::consts::FRAC_1_SQRT_2),
        )
    }

    /// Left-circular polarization.
    pub fn circular_left() -> Self {
        JonesVector::new(
            Complex::real(core::f64::consts::FRAC_1_SQRT_2),
            Complex::new(0.0, core::f64::consts::FRAC_1_SQRT_2),
        )
    }

    /// Intensity `|E_x|² + |E_y|²`.
    pub fn intensity(&self) -> f64 {
        self.ex.norm_sqr() + self.ey.norm_sqr()
    }

    /// Apply a Jones matrix.
    #[must_use]
    pub fn apply(&self, m: &JonesMatrix) -> JonesVector {
        JonesVector::new(m.a * self.ex + m.b * self.ey, m.c * self.ex + m.d * self.ey)
    }
}

/// A 2×2 complex Jones matrix (a polarization operator).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct JonesMatrix {
    /// Element (0,0).
    pub a: Complex<f64>,
    /// Element (0,1).
    pub b: Complex<f64>,
    /// Element (1,0).
    pub c: Complex<f64>,
    /// Element (1,1).
    pub d: Complex<f64>,
}

impl JonesMatrix {
    /// Construct from row-major elements.
    pub const fn new(a: Complex<f64>, b: Complex<f64>, c: Complex<f64>, d: Complex<f64>) -> Self {
        JonesMatrix { a, b, c, d }
    }

    /// Identity.
    pub fn identity() -> Self {
        JonesMatrix::new(
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::one(),
        )
    }

    /// Matrix product `self · other` (apply `other` first).
    #[must_use]
    pub fn mul_matrix(&self, other: &JonesMatrix) -> JonesMatrix {
        JonesMatrix::new(
            self.a * other.a + self.b * other.c,
            self.a * other.b + self.b * other.d,
            self.c * other.a + self.d * other.c,
            self.c * other.b + self.d * other.d,
        )
    }

    /// Rotate the operator's axes by `theta` (rad): `J(θ) = R(θ)·J·R(−θ)`,
    /// where `R` is the field rotation. An x-axis element becomes one with
    /// its axis at +θ, matching [`JonesMatrix::polarizer`].
    #[must_use]
    pub fn rotated(&self, theta: f64) -> JonesMatrix {
        let r = rotation(theta);
        let r_inv = rotation(-theta);
        r.mul_matrix(&self.mul_matrix(&r_inv))
    }

    /// Linear ideal polarizer with transmission axis at `theta` (rad).
    #[must_use]
    pub fn polarizer(theta: f64) -> JonesMatrix {
        let ct = libm::cos(theta);
        let st = libm::sin(theta);
        JonesMatrix::new(
            Complex::real(ct * ct),
            Complex::real(ct * st),
            Complex::real(ct * st),
            Complex::real(st * st),
        )
    }

    /// Retarder with phase delay `delta` (rad, slow minus fast axis) and
    /// fast axis at `theta` (rad): `J = R(θ)·diag(e^{iδ/2}?, ...)` in the
    /// common half-delay convention `diag(1, e^{iδ})` rotated.
    #[must_use]
    pub fn waveplate(delta: f64, theta: f64) -> JonesMatrix {
        // In the fast-axis frame: [[1, 0], [0, e^{iδ}]].
        let retard = Complex::cis(delta);
        JonesMatrix::new(Complex::one(), Complex::zero(), Complex::zero(), retard).rotated(theta)
    }

    /// Quarter-wave plate with fast axis at `theta`.
    #[must_use]
    pub fn quarter_wave(theta: f64) -> JonesMatrix {
        JonesMatrix::waveplate(core::f64::consts::FRAC_PI_2, theta)
    }

    /// Half-wave plate with fast axis at `theta`.
    #[must_use]
    pub fn half_wave(theta: f64) -> JonesMatrix {
        JonesMatrix::waveplate(core::f64::consts::PI, theta)
    }

    /// Optical rotator: rotates the polarization by `theta`.
    #[must_use]
    pub fn rotator(theta: f64) -> JonesMatrix {
        rotation(theta)
    }

    /// Linear diattenuator with amplitude transmissions `t_fast`, `t_slow`
    /// along the axes at `theta` (fast axis at `theta`).
    #[must_use]
    pub fn diattenuator(t_fast: f64, t_slow: f64, theta: f64) -> JonesMatrix {
        JonesMatrix::new(
            Complex::real(t_fast),
            Complex::zero(),
            Complex::zero(),
            Complex::real(t_slow),
        )
        .rotated(theta)
    }
}

/// 2-D rotation matrix `R(θ) = [[cosθ, −sinθ], [sinθ, cosθ]]`: rotates the
/// polarization azimuth by +θ.
#[must_use]
pub fn rotation(theta: f64) -> JonesMatrix {
    let c = libm::cos(theta);
    let s = libm::sin(theta);
    JonesMatrix::new(
        Complex::real(c),
        Complex::real(-s),
        Complex::real(s),
        Complex::real(c),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn malus_law() {
        // Crossed polarizers block; at 45° transmits half.
        let p0 = JonesMatrix::polarizer(0.0);
        let p90 = JonesMatrix::polarizer(core::f64::consts::FRAC_PI_2);
        let p45 = JonesMatrix::polarizer(core::f64::consts::FRAC_PI_4);
        let out_cross = JonesVector::linear_x().apply(&p90.mul_matrix(&p0));
        assert!(out_cross.intensity() < 1e-15);
        let out45 = JonesVector::linear_x().apply(&p45.mul_matrix(&p0));
        assert!(close(out45.intensity(), 0.5, 1e-15));
        // Parallel: full transmission.
        let out_par = JonesVector::linear_x().apply(&p0.mul_matrix(&p0));
        assert!(close(out_par.intensity(), 1.0, 1e-15));
    }

    #[test]
    fn half_wave_flips_linear_polarization() {
        // HWP with fast axis at θ maps linear angle φ → 2θ − φ.
        let hwp = JonesMatrix::half_wave(0.0); // axis along x: flips y sign
        let out = JonesVector::linear_at(core::f64::consts::FRAC_PI_3).apply(&hwp);
        let expect = JonesVector::linear_at(-core::f64::consts::FRAC_PI_3);
        let overlap = out.ex * expect.ex.conj() + out.ey * expect.ey.conj();
        assert!(close(overlap.abs(), out.intensity().sqrt(), 1e-12));
        // HWP at 45° turns x into y.
        let hwp45 = JonesMatrix::half_wave(core::f64::consts::FRAC_PI_4);
        let out2 = JonesVector::linear_x().apply(&hwp45);
        assert!(out2.ey.norm_sqr().sqrt() > 0.999);
        assert!(out2.ex.abs() < 1e-12);
    }

    #[test]
    fn quarter_wave_makes_circular() {
        // Linear at 45° through QWP (axis x) → circular.
        let qwp = JonesMatrix::quarter_wave(0.0);
        let out = JonesVector::linear_at(core::f64::consts::FRAC_PI_4).apply(&qwp);
        assert!(
            close(out.ex.abs(), out.ey.abs(), 1e-12),
            "|ex|={} |ey|={}",
            out.ex.abs(),
            out.ey.abs()
        );
        let phase = (out.ey * out.ex.conj()).arg();
        assert!(close(phase.abs(), core::f64::consts::FRAC_PI_2, 1e-12));
    }

    #[test]
    fn rotator_rotation() {
        let r = JonesMatrix::rotator(0.7);
        let out = JonesVector::linear_x().apply(&r);
        let expect = JonesVector::linear_at(0.7);
        let overlap = out.ex * expect.ex.conj() + out.ey * expect.ey.conj();
        assert!(close(overlap.re / out.intensity(), 1.0, 1e-12));
    }

    #[test]
    fn rotated_waveplate_equivalence() {
        // A waveplate at angle θ equals the rotated canonical retarder.
        let canonical = JonesMatrix::new(
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::cis(1.234),
        );
        let theta = 0.6;
        let rotated = canonical.rotated(theta);
        let direct = JonesMatrix::waveplate(1.234, theta);
        for (v, w) in [
            (rotated.a, direct.a),
            (rotated.b, direct.b),
            (rotated.c, direct.c),
            (rotated.d, direct.d),
        ] {
            assert!(v.approx_eq(w, 1e-14), "{v} vs {w}");
        }
    }

    #[test]
    fn diattenuator_transmissions() {
        let d = JonesMatrix::diattenuator(1.0, 0.5, 0.0);
        let out_fast = JonesVector::linear_x().apply(&d);
        assert!(close(out_fast.intensity(), 1.0, 1e-15));
        let out_slow = JonesVector::linear_y().apply(&d);
        assert!(close(out_slow.intensity(), 0.25, 1e-15));
    }
}

#[cfg(test)]
mod rotation_convention_tests {
    use super::*;
    use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    #[test]
    fn rotated_polarizer_transmits_along_plus_theta() {
        // A polarizer with its axis at +θ must pass light polarized at +θ
        // and block light at +θ + 90°.
        let theta = 0.4;
        let x_axis = JonesMatrix::new(
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::zero(),
        );
        let p = x_axis.rotated(theta);
        let pass = JonesVector::linear_at(theta).apply(&p);
        assert!(
            (pass.intensity() - 1.0).abs() < 1e-12,
            "T = {}",
            pass.intensity()
        );
        let block = JonesVector::linear_at(theta + FRAC_PI_2).apply(&p);
        assert!(block.intensity() < 1e-12);
    }

    #[test]
    fn rotated_polarizer_matches_constructor() {
        let theta = FRAC_PI_4 * 0.7;
        let x_axis = JonesMatrix::new(
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::zero(),
        );
        let a = x_axis.rotated(theta);
        let b = JonesMatrix::polarizer(theta);
        for (u, v) in [(a.a, b.a), (a.b, b.b), (a.c, b.c), (a.d, b.d)] {
            assert!((u - v).abs() < 1e-12, "rotated {u:?} vs polarizer {v:?}");
        }
    }
}

//! Mueller calculus: 4×4 real matrices acting on Stokes vectors, including
//! depolarizing elements that Jones calculus cannot express.

use crate::jones::{JonesMatrix, JonesVector};
use crate::stokes::StokesVector;
use tpt_optics_core::Complex;

/// A 4×4 Mueller matrix, row-major.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MuellerMatrix {
    /// Row-major 16 elements.
    pub m: [f64; 16],
}

impl MuellerMatrix {
    /// Construct from row-major elements.
    pub const fn new(m: [f64; 16]) -> Self {
        MuellerMatrix { m }
    }

    /// The identity (transparent, non-depolarizing).
    pub fn identity() -> Self {
        let mut m = [0.0; 16];
        for i in 0..4 {
            m[i * 4 + i] = 1.0;
        }
        MuellerMatrix::new(m)
    }

    /// Apply to a Stokes vector.
    #[must_use]
    pub fn apply(&self, s: &StokesVector) -> StokesVector {
        let v = [s.s0, s.s1, s.s2, s.s3];
        let mut out = [0.0; 4];
        for (r, o) in out.iter_mut().enumerate() {
            for (c, m_rc) in self.m[r * 4..r * 4 + 4].iter().enumerate() {
                *o += m_rc * v[c];
            }
        }
        StokesVector::new(out[0], out[1], out[2], out[3])
    }

    /// Matrix product `self · other` (apply `other` first).
    #[must_use]
    pub fn mul_matrix(&self, other: &MuellerMatrix) -> MuellerMatrix {
        let mut out = [0.0; 16];
        for r in 0..4 {
            for c in 0..4 {
                let mut acc = 0.0;
                for k in 0..4 {
                    acc += self.m[r * 4 + k] * other.m[k * 4 + c];
                }
                out[r * 4 + c] = acc;
            }
        }
        MuellerMatrix::new(out)
    }

    /// Rotate the element's axes by `theta` (rad):
    /// `M(θ) = R(θ)·M·R(−θ)` — conjugating an axis-aligned canonical
    /// element with the azimuth rotation `R` (this maps an x-axis element
    /// to one whose axis sits at azimuth +θ).
    #[must_use]
    pub fn rotated(&self, theta: f64) -> MuellerMatrix {
        rotation(theta).mul_matrix(&self.mul_matrix(&rotation(-theta)))
    }

    /// Ideal linear polarizer with axis at `theta` (rad).
    #[must_use]
    pub fn polarizer(theta: f64) -> Self {
        let mut m = [0.0; 16];
        m[0] = 0.5;
        m[1] = 0.5;
        m[4] = 0.5;
        m[5] = 0.5;
        MuellerMatrix::new(m).rotated(theta)
    }

    /// Retarder of phase delay `delta` (rad) with fast axis at `theta` (rad).
    #[must_use]
    pub fn waveplate(delta: f64, theta: f64) -> Self {
        // Canonical retarder (fast axis x).
        let c = libm::cos(delta);
        let s = libm::sin(delta);
        let m = [
            1.0, 0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, 0.0, //
            0.0, 0.0, c, s, //
            0.0, 0.0, -s, c, //
        ];
        MuellerMatrix::new(m).rotated(theta)
    }

    /// Quarter-wave plate at `theta`.
    #[must_use]
    pub fn quarter_wave(theta: f64) -> Self {
        MuellerMatrix::waveplate(core::f64::consts::FRAC_PI_2, theta)
    }

    /// Half-wave plate at `theta`.
    #[must_use]
    pub fn half_wave(theta: f64) -> Self {
        MuellerMatrix::waveplate(core::f64::consts::PI, theta)
    }

    /// Linear diattenuator with intensity transmissions `t_fast`, `t_slow`.
    #[must_use]
    pub fn diattenuator(t_fast: f64, t_slow: f64, theta: f64) -> Self {
        let a = (t_fast + t_slow) / 2.0;
        let b = (t_fast - t_slow) / 2.0;
        let m = [
            a, b, 0.0, 0.0, //
            b, a, 0.0, 0.0, //
            0.0, 0.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, 0.0, //
        ];
        MuellerMatrix::new(m).rotated(theta)
    }

    /// Deterministic Mueller matrix of a Jones matrix (no depolarization),
    /// derived through the 2×2 coherency (density) matrix so that the
    /// `(I, Q, U, V)` conventions match [`StokesVector::from_jones`]
    /// exactly by construction.
    #[must_use]
    pub fn from_jones(j: &JonesMatrix) -> Self {
        // Jones action on the coherency matrix: C' = J C J† (congruence,
        // correct for mixed states too — the linear extension that the
        // Mueller matrix columns need).
        let act = |c: [[Complex<f64>; 2]; 2]| -> [[Complex<f64>; 2]; 2] {
            let jj = [[j.a, j.b], [j.c, j.d]];
            let mut out = [[Complex::zero(); 2]; 2];
            for (i, row) in out.iter_mut().enumerate() {
                for (col, v) in row.iter_mut().enumerate() {
                    let mut acc = Complex::zero();
                    for (k, jk) in jj[i].iter().enumerate() {
                        for (l, cl) in c[k].iter().enumerate() {
                            acc += *jk * *cl * jj[col][l].conj();
                        }
                    }
                    *v = acc;
                }
            }
            out
        };
        // Stokes ↔ coherency: C = ½[[s0+s1, s2−i s3], [s2+i s3, s0−s1]].
        let stokes_to_c = |v: [f64; 4]| -> [[Complex<f64>; 2]; 2] {
            [
                [
                    Complex::real(0.5 * (v[0] + v[1])),
                    Complex::new(0.5 * v[2], -0.5 * v[3]),
                ],
                [
                    Complex::new(0.5 * v[2], 0.5 * v[3]),
                    Complex::real(0.5 * (v[0] - v[1])),
                ],
            ]
        };
        let c_to_stokes = |c: [[Complex<f64>; 2]; 2]| -> [f64; 4] {
            [
                (c[0][0] + c[1][1]).re,
                (c[0][0] - c[1][1]).re,
                2.0 * c[0][1].re,
                -2.0 * c[0][1].im,
            ]
        };
        // Column k of the Mueller matrix = image of the k-th Stokes basis.
        let mut m = [0.0; 16];
        for (k, basis) in [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
        .iter()
        .enumerate()
        {
            let out = c_to_stokes(act(stokes_to_c(*basis)));
            for r in 0..4 {
                m[r * 4 + k] = out[r];
            }
        }
        MuellerMatrix::new(m)
    }

    /// Ideal depolarizer: scales S₀ by `dop` retention... constructs the
    /// depolarizing matrix that maps any input to unpolarized light of the
    /// same total intensity.
    #[must_use]
    pub fn ideal_depolarizer() -> Self {
        let mut m = [0.0; 16];
        m[0] = 1.0;
        MuellerMatrix::new(m)
    }

    /// Partial depolarizer retaining degree of polarization factor `f`
    /// relative to the input DOP: `DOP_out = f · DOP_in` (s0 unchanged).
    #[must_use]
    pub fn partial_depolarizer(f: f64) -> Self {
        let mut m = [0.0; 16];
        m[0] = 1.0;
        m[5] = f;
        m[10] = f;
        m[15] = f;
        MuellerMatrix::new(m)
    }
}

/// Stokes-space azimuth rotation `R(θ)`: maps a state at azimuth ψ to
/// ψ + θ (mixes Q and U with the opposite sign of a passive frame change).
#[must_use]
pub fn rotation(theta: f64) -> MuellerMatrix {
    let c = libm::cos(2.0 * theta);
    let s = libm::sin(2.0 * theta);
    let mut m = [0.0; 16];
    m[0] = 1.0;
    m[5] = c;
    m[6] = -s;
    m[9] = s;
    m[10] = c;
    m[15] = 1.0;
    MuellerMatrix::new(m)
}

/// Convenience: Jones vector through a Mueller element (via Stokes).
#[must_use]
pub fn apply_to_jones(m: &MuellerMatrix, j: &JonesVector) -> (f64, f64, f64, f64) {
    let s = m.apply(&StokesVector::from_jones(j));
    (s.s0, s.s1, s.s2, s.s3)
}

#[cfg(test)]
mod tests {
    use super::*;

    use alloc::vec;
    use alloc::vec::Vec;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn malus_law_mueller() {
        let p0 = MuellerMatrix::polarizer(0.0);
        let out = p0.apply(&StokesVector::unpolarized(1.0));
        assert!(close(out.s0, 0.5, 1e-15));
        assert!(close(out.s1, 0.5, 1e-15));
        // Crossed: zero.
        let p90 = MuellerMatrix::polarizer(core::f64::consts::FRAC_PI_2);
        let dark = p90.apply(&out);
        assert!(dark.s0.abs() < 1e-15);
    }

    #[test]
    fn hwp_rotates_azimuth() {
        // HWP at θ maps azimuth φ → 2θ − φ (Stokes U flips sign at θ=0).
        let hwp = MuellerMatrix::half_wave(0.0);
        let s = StokesVector::linear(0.4, 1.0);
        let out = hwp.apply(&s);
        assert!(close(out.s1, s.s1, 1e-12));
        assert!(close(out.s2, -s.s2, 1e-12));
        // QWP turns linear ±45 into circular.
        let qwp = MuellerMatrix::quarter_wave(0.0);
        let s45 = StokesVector::linear(core::f64::consts::FRAC_PI_4, 1.0);
        let circ = qwp.apply(&s45);
        assert!(close(circ.s3.abs(), 1.0, 1e-12));
    }

    #[test]
    fn jones_mueller_consistency() {
        // For deterministic elements, Jones- and Mueller-calculus results
        // must agree on a battery of inputs.
        let elements: Vec<JonesMatrix> = vec![
            JonesMatrix::polarizer(0.3),
            JonesMatrix::waveplate(1.1, -0.4),
            JonesMatrix::diattenuator(1.0, 0.3, 0.9),
            JonesMatrix::rotator(0.77),
        ];
        let inputs = vec![
            JonesVector::linear_x(),
            JonesVector::linear_at(1.2),
            JonesVector::circular_left(),
            JonesVector::new(Complex::real(0.6), Complex::new(0.3, 0.2)),
        ];
        for je in &elements {
            let me = MuellerMatrix::from_jones(je);
            for j_in in &inputs {
                let s_in = StokesVector::from_jones(j_in);
                let via_jones = StokesVector::from_jones(&j_in.apply(je));
                let via_mueller = me.apply(&s_in);
                for k in 0..4 {
                    let a = [via_jones.s0, via_jones.s1, via_jones.s2, via_jones.s3][k];
                    let b = [
                        via_mueller.s0,
                        via_mueller.s1,
                        via_mueller.s2,
                        via_mueller.s3,
                    ][k];
                    assert!(close(a, b, 1e-12), "element/input {k}: {a} vs {b}");
                }
            }
        }
    }

    #[test]
    fn depolarizers() {
        let s = StokesVector::linear(0.2, 1.0);
        let out = MuellerMatrix::ideal_depolarizer().apply(&s);
        assert!(close(out.s0, 1.0, 1e-15));
        assert!(close(out.degree_of_polarization(), 0.0, 1e-15));
        let out5 = MuellerMatrix::partial_depolarizer(0.5).apply(&s);
        assert!(close(out5.degree_of_polarization(), 0.5, 1e-12));
    }

    #[test]
    fn rotation_matrix_acting() {
        // R(θ) maps a state at azimuth ψ to ψ + θ.
        let s = StokesVector::linear(0.3, 1.0);
        let a = rotation(0.2).apply(&s);
        let expect = StokesVector::linear(0.5, 1.0);
        assert!(close(a.s1, expect.s1, 1e-12));
        assert!(close(a.s2, expect.s2, 1e-12));
    }
}

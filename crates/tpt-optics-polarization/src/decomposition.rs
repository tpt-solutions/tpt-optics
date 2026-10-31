//! Lu–Chipman decomposition of a Mueller matrix into a depolarizer, a
//! retarder and a diattenuator, `M = M_Δ · M_R · M_D` (Lu & Chipman, J. Opt.
//! Soc. Am. A 13, 1996).
//!
//! The diattenuator is taken from the first row of `M`. The depolarizer and
//! retarder come from the polar decomposition of the 3 × 3 block that remains
//! once `M_D` is removed, `m = m_Δ · R`, with `m_Δ = (m mᵀ)^{1/2}`.

use crate::mueller::MuellerMatrix;

/// The three factors of a Mueller matrix.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LuChipman {
    /// Diattenuation vector `D` (the first row of `M`, divided by `m₀₀`).
    pub diattenuation: [f64; 3],
    /// Depolarizer block `m_Δ`: symmetric and positive semidefinite. It is the
    /// identity for a lossless system, and `f·I` for an ideal partial depolarizer.
    pub depolarizer: [[f64; 3]; 3],
    /// Retarder block `R`, an orthogonal 3 × 3 matrix.
    pub retarder: [[f64; 3]; 3],
}

/// Decompose `m` into depolarizer, retarder and diattenuator.
///
/// Returns `None` when `m₀₀ ≤ 0`, when the diattenuation reaches one, or when
/// the remaining block is singular. In those cases the decomposition is not
/// defined.
#[must_use]
pub fn lu_chipman(m: &MuellerMatrix) -> Option<LuChipman> {
    let a = |r: usize, c: usize| m.m[r * 4 + c];
    let m00 = a(0, 0);
    if !(m00 > 0.0 && m00.is_finite()) {
        return None;
    }
    let d = [a(0, 1) / m00, a(0, 2) / m00, a(0, 3) / m00];
    let d2: f64 = d.iter().map(|v| v * v).sum();
    if d2 >= 1.0 {
        return None;
    }
    // Pure diattenuator with unit m₀₀: [[1, Dᵀ], [D, m_D]], m_D = s·I + b·D·Dᵀ.
    let s = libm::sqrt(1.0 - d2);
    let b = if d2 > 0.0 { (1.0 - s) / d2 } else { 0.0 };
    let mut md = [0.0; 16];
    md[0] = 1.0;
    for i in 0..3 {
        md[i + 1] = d[i];
        md[(i + 1) * 4] = d[i];
        for j in 0..3 {
            let delta = if i == j { s } else { 0.0 };
            md[(i + 1) * 4 + (j + 1)] = delta + b * d[i] * d[j];
        }
    }
    let md_inv = invert4(&md)?;
    // X = (M / m₀₀) · M_D⁻¹ has no diattenuation, so its first row is (x, 0, 0, 0).
    let normalised: [f64; 16] = core::array::from_fn(|i| m.m[i] / m00);
    let x = mul4(&normalised, &md_inv);
    let mut block = [[0.0; 3]; 3];
    for (i, row) in block.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = x[(i + 1) * 4 + (j + 1)];
        }
    }
    // m_Δ = (m mᵀ)^{1/2}, R = m_Δ⁻¹ m.
    let g = mul3(&block, &transpose3(&block));
    let (lambda, vecs) = jacobi3(g);
    if lambda.iter().any(|l| l.is_nan() || *l <= 0.0) {
        return None;
    }
    let sqrt_g = map_sym(&vecs, &lambda, libm::sqrt);
    let inv_sqrt_g = map_sym(&vecs, &lambda, |l| 1.0 / libm::sqrt(l));
    let retarder = mul3(&inv_sqrt_g, &block);
    Some(LuChipman {
        diattenuation: d,
        depolarizer: sqrt_g,
        retarder,
    })
}

impl LuChipman {
    /// Diattenuation magnitude `|D|`.
    #[must_use]
    pub fn diattenuation_magnitude(&self) -> f64 {
        libm::sqrt(self.diattenuation.iter().map(|v| v * v).sum())
    }

    /// Retardance (rad) from the trace of the retarder: `tr R = 1 + 2 cos δ`.
    /// The cosine is clamped to `[−1, 1]` before the arccosine.
    #[must_use]
    pub fn retardance(&self) -> f64 {
        let tr = self.retarder[0][0] + self.retarder[1][1] + self.retarder[2][2];
        libm::acos(((tr - 1.0) / 2.0).clamp(-1.0, 1.0))
    }
}

/// `f(M) = Σ f(λᵢ) vᵢ vᵢᵀ` for a symmetric matrix with eigenpairs `(λᵢ, vᵢ)`.
fn map_sym(vecs: &[[f64; 3]; 3], lambda: &[f64; 3], f: impl Fn(f64) -> f64) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for (k, &lk) in lambda.iter().enumerate() {
        let w = f(lk);
        for (i, row) in out.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v += w * vecs[i][k] * vecs[j][k];
            }
        }
    }
    out
}

/// Eigen-decomposition of a symmetric 3 × 3 matrix by cyclic Jacobi rotations.
/// Returns the eigenvalues and the eigenvectors as the columns of a matrix.
fn jacobi3(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..100 {
        let total: f64 = a.iter().flatten().map(|x| x * x).sum();
        let off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if off <= 1e-30 * total || total == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (libm::fabs(theta) + libm::sqrt(theta * theta + 1.0));
            let c = 1.0 / libm::sqrt(t * t + 1.0);
            let s = t * c;
            // A ← Jᵀ A J: columns p and q, then rows p and q; V ← V J.
            for row in &mut a {
                let (akp, akq) = (row[p], row[q]);
                row[p] = c * akp - s * akq;
                row[q] = s * akp + c * akq;
            }
            let (row_p, row_q) = (a[p], a[q]);
            for k in 0..3 {
                a[p][k] = c * row_p[k] - s * row_q[k];
                a[q][k] = s * row_p[k] + c * row_q[k];
            }
            for row in &mut v {
                let (vkp, vkq) = (row[p], row[q]);
                row[p] = c * vkp - s * vkq;
                row[q] = s * vkp + c * vkq;
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

fn mul3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    out
}

fn transpose3(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    core::array::from_fn(|i| core::array::from_fn(|j| a[j][i]))
}

fn mul4(a: &[f64; 16], b: &[f64; 16]) -> [f64; 16] {
    let mut out = [0.0; 16];
    for r in 0..4 {
        for c in 0..4 {
            out[r * 4 + c] = (0..4).map(|k| a[r * 4 + k] * b[k * 4 + c]).sum();
        }
    }
    out
}

/// Inverse of a 4 × 4 row-major matrix by Gauss–Jordan elimination with
/// partial pivoting. Returns `None` when the matrix is singular.
fn invert4(m: &[f64; 16]) -> Option<[f64; 16]> {
    let mut a = *m;
    let mut inv = [0.0; 16];
    for i in 0..4 {
        inv[i * 4 + i] = 1.0;
    }
    for col in 0..4 {
        let pivot = (col..4)
            .max_by(|&i, &j| libm::fabs(a[i * 4 + col]).total_cmp(&libm::fabs(a[j * 4 + col])))?;
        if libm::fabs(a[pivot * 4 + col]) < 1e-300 {
            return None;
        }
        for k in 0..4 {
            a.swap(col * 4 + k, pivot * 4 + k);
            inv.swap(col * 4 + k, pivot * 4 + k);
        }
        let p = a[col * 4 + col];
        for k in 0..4 {
            a[col * 4 + k] /= p;
            inv[col * 4 + k] /= p;
        }
        for row in 0..4 {
            if row != col {
                let f = a[row * 4 + col];
                for k in 0..4 {
                    a[row * 4 + k] -= f * a[col * 4 + k];
                    inv[row * 4 + k] -= f * inv[col * 4 + k];
                }
            }
        }
    }
    Some(inv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_the_factors_of_a_known_product() {
        // M = M_Δ · M_R · M_D with a partial depolarizer, a waveplate and a
        // diattenuator, so the three factors are known in advance.
        let f = 0.8;
        let dep = MuellerMatrix::partial_depolarizer(f);
        let ret = MuellerMatrix::waveplate(0.9, 0.3);
        let dia = MuellerMatrix::diattenuator(0.9, 0.5, 0.0);
        let m = dep.mul_matrix(&ret).mul_matrix(&dia);
        let lc = lu_chipman(&m).expect("decomposable");

        let expect_d = [
            dia.m[1] / dia.m[0],
            dia.m[2] / dia.m[0],
            dia.m[3] / dia.m[0],
        ];
        for (i, (got, want)) in lc.diattenuation.iter().zip(expect_d).enumerate() {
            assert!((got - want).abs() < 1e-10, "D[{i}]");
        }
        for (i, row) in lc.depolarizer.iter().enumerate() {
            for (j, got) in row.iter().enumerate() {
                let want = if i == j { f } else { 0.0 };
                assert!((got - want).abs() < 1e-9, "m_Δ[{i}][{j}]");
            }
        }
        let tr_ret = ret.m[5] + ret.m[10] + ret.m[15];
        let tr_rec = lc.retarder[0][0] + lc.retarder[1][1] + lc.retarder[2][2];
        assert!((tr_rec - tr_ret).abs() < 1e-9, "trace {tr_rec} vs {tr_ret}");
        // The recovered retarder is orthogonal.
        let rrt = mul3(&lc.retarder, &transpose3(&lc.retarder));
        for (i, row) in rrt.iter().enumerate() {
            for (j, got) in row.iter().enumerate() {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((got - want).abs() < 1e-9, "R Rᵀ [{i}][{j}]");
            }
        }
    }

    #[test]
    fn retardance_matches_the_waveplate() {
        let ret = MuellerMatrix::waveplate(1.2, 0.0);
        let lc = lu_chipman(&ret).unwrap();
        assert!((lc.retardance() - 1.2).abs() < 1e-9, "{}", lc.retardance());
        assert!(lc.diattenuation_magnitude() < 1e-12);
    }

    #[test]
    fn non_positive_m00_is_refused() {
        assert!(lu_chipman(&MuellerMatrix::new([0.0; 16])).is_none());
    }
}

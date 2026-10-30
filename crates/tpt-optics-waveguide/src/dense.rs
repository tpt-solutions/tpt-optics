//! Dense linear-algebra helpers shared by the mode solvers.
//!
//! The solvers assemble dense systems of a few hundred unknowns. They factor
//! once and then solve for several right-hand sides.

use alloc::vec::Vec;

/// Pivots below this magnitude are treated as singular.
const SINGULAR_PIVOT: f64 = 1e-300;

/// LU factorisation with partial pivoting of a dense `n × n` matrix
/// (row-major), stored for repeated solves.
pub(crate) struct LuFactor {
    lu: Vec<f64>,
    piv: Vec<usize>,
    n: usize,
}

impl LuFactor {
    /// Factorises `a`. Returns `None` if a pivot is zero, non-finite, or
    /// smaller than `1e-300` in magnitude. Each case means the matrix is
    /// singular to working precision.
    pub(crate) fn new(mut a: Vec<f64>, n: usize) -> Option<Self> {
        debug_assert_eq!(a.len(), n * n);
        let mut piv: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let mut best = k;
            for i in (k + 1)..n {
                if libm::fabs(a[i * n + k]) > libm::fabs(a[best * n + k]) {
                    best = i;
                }
            }
            if best != k {
                for j in 0..n {
                    a.swap(k * n + j, best * n + j);
                }
                piv.swap(k, best);
            }
            let pivot = a[k * n + k];
            let magnitude = libm::fabs(pivot);
            if magnitude.is_nan() || magnitude < SINGULAR_PIVOT {
                return None;
            }
            for i in (k + 1)..n {
                let f = a[i * n + k] / pivot;
                a[i * n + k] = f;
                for j in (k + 1)..n {
                    a[i * n + j] -= f * a[k * n + j];
                }
            }
        }
        Some(Self { lu: a, piv, n })
    }

    /// Solves `A·x = b` using the stored factors.
    pub(crate) fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x: Vec<f64> = self.piv.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            for k in 0..i {
                x[i] -= self.lu[i * n + k] * x[k];
            }
        }
        for i in (0..n).rev() {
            for k in (i + 1)..n {
                x[i] -= self.lu[i * n + k] * x[k];
            }
            x[i] /= self.lu[i * n + i];
        }
        x
    }
}

/// Product `A·x` for a dense row-major `n × n` matrix.
pub(crate) fn mat_vec(a: &[f64], n: usize, x: &[f64]) -> Vec<f64> {
    (0..n)
        .map(|i| {
            a[i * n..(i + 1) * n]
                .iter()
                .zip(x)
                .map(|(m, v)| m * v)
                .sum()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn solves_a_pivoting_system() {
        // First pivot is zero, so the solver must swap rows.
        let a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 2.0, 0.0, 3.0];
        let x_true = [1.0, -2.0, 0.5];
        let b = mat_vec(&a, 3, &x_true);
        let x = LuFactor::new(a, 3).unwrap().solve(&b);
        for (got, want) in x.iter().zip(x_true) {
            assert!((got - want).abs() < 1e-12, "{x:?}");
        }
    }

    #[test]
    fn singular_matrix_is_rejected() {
        let a = vec![1.0, 2.0, 2.0, 4.0];
        assert!(LuFactor::new(a, 2).is_none());
    }

    #[test]
    fn nan_pivot_is_rejected() {
        let a = vec![f64::NAN, 1.0, 1.0, 1.0];
        assert!(LuFactor::new(a, 2).is_none());
    }
}

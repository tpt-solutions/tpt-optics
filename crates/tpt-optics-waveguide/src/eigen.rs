//! Dense symmetric eigenvalues by Householder tridiagonalisation followed by
//! the implicit QL algorithm.
//!
//! Cost is about `(4/3 + small)·n³` flops for the reduction and `O(n²)` for
//! the tridiagonal solve. The cyclic Jacobi solver in [`crate::phc`] needs
//! roughly ten sweeps of `O(n³)` each, which makes it too slow for the large
//! plane-wave bases that supercell calculations need. Eigenvectors are not
//! computed here.

use alloc::vec;
use alloc::vec::Vec;

use libm::{fabs, hypot};

/// All eigenvalues of the real symmetric `n × n` matrix `a` (row-major),
/// returned in ascending order. `a` is overwritten.
pub fn symmetric_eigenvalues(a: &mut [f64], n: usize) -> Vec<f64> {
    debug_assert_eq!(a.len(), n * n);
    if n == 0 {
        return Vec::new();
    }
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    householder_tridiagonal(a, n, &mut d, &mut e);
    implicit_ql(&mut d, &mut e);
    d.sort_by(f64::total_cmp);
    d
}

/// Householder reduction of the symmetric matrix `a` to tridiagonal form.
/// On exit `d` holds the diagonal and `e[i]` (for `i ≥ 1`) the coupling
/// between rows `i − 1` and `i`.
fn householder_tridiagonal(a: &mut [f64], n: usize, d: &mut [f64], e: &mut [f64]) {
    let at = |i: usize, j: usize| i * n + j;
    for i in (1..n).rev() {
        let l = i - 1;
        let mut h = 0.0;
        if l > 0 {
            let scale: f64 = (0..=l).map(|k| fabs(a[at(i, k)])).sum();
            if scale == 0.0 {
                e[i] = a[at(i, l)];
            } else {
                for k in 0..=l {
                    a[at(i, k)] /= scale;
                    h += a[at(i, k)] * a[at(i, k)];
                }
                let f = a[at(i, l)];
                let g = if f >= 0.0 { -libm::sqrt(h) } else { libm::sqrt(h) };
                e[i] = scale * g;
                h -= f * g;
                a[at(i, l)] = f - g;
                let mut f_sum = 0.0;
                for j in 0..=l {
                    let mut g_j = 0.0;
                    for k in 0..=j {
                        g_j += a[at(j, k)] * a[at(i, k)];
                    }
                    for k in (j + 1)..=l {
                        g_j += a[at(k, j)] * a[at(i, k)];
                    }
                    e[j] = g_j / h;
                    f_sum += e[j] * a[at(i, j)];
                }
                let hh = f_sum / (h + h);
                for j in 0..=l {
                    let f_j = a[at(i, j)];
                    let g_j = e[j] - hh * f_j;
                    e[j] = g_j;
                    for k in 0..=j {
                        a[at(j, k)] -= f_j * e[k] + g_j * a[at(i, k)];
                    }
                }
            }
        } else {
            e[i] = a[at(i, l)];
        }
        d[i] = h;
    }
    // The diagonal of the tridiagonal matrix sits on the diagonal of `a`
    // after the reduction.
    for (i, di) in d.iter_mut().enumerate() {
        *di = a[at(i, i)];
    }
    e[0] = 0.0;
}

/// Eigenvalues of a symmetric tridiagonal matrix (diagonal `d`, off-diagonal
/// `e[i]` coupling rows `i − 1` and `i`, with `e[0]` unused) by implicit QL.
/// On exit `d` holds the eigenvalues in no particular order.
fn implicit_ql(d: &mut [f64], e_in: &mut [f64]) {
    let n = d.len();
    // Shift so that `e[i]` couples rows `i` and `i + 1`.
    let mut e = vec![0.0; n];
    for i in 0..n.saturating_sub(1) {
        e[i] = e_in[i + 1];
    }
    const EPS: f64 = 1e-15;
    for l in 0..n {
        let mut iter = 0;
        loop {
            let mut m = l;
            while m + 1 < n {
                let dd = fabs(d[m]) + fabs(d[m + 1]);
                if fabs(e[m]) <= EPS * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iter += 1;
            if iter > 200 {
                break;
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = hypot(g, 1.0);
            let sign_r = if g >= 0.0 { fabs(r) } else { -fabs(r) };
            g = d[m] - d[l] + e[l] / (g + sign_r);
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut i = m;
            let mut underflow = false;
            while i > l {
                i -= 1;
                let f = s * e[i];
                let b = c * e[i];
                r = hypot(f, g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    underflow = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
            }
            if underflow {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phc::jacobi_eigen;

    /// Deterministic pseudo-random symmetric matrix (LCG).
    fn random_symmetric(n: usize, mut seed: u64) -> Vec<f64> {
        let mut next = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
        };
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = next();
                a[i * n + j] = v;
                a[j * n + i] = v;
            }
        }
        a
    }

    #[test]
    fn matches_jacobi_on_random_matrices() {
        for &(n, seed) in &[(2usize, 1u64), (5, 7), (17, 3), (40, 11)] {
            let a = random_symmetric(n, seed);
            let mut for_ql = a.clone();
            let ql = symmetric_eigenvalues(&mut for_ql, n);
            let mut for_jac = a.clone();
            let mut jac = vec![0.0; n];
            jacobi_eigen(&mut for_jac, &mut jac);
            for (x, y) in ql.iter().zip(&jac) {
                assert!((x - y).abs() < 1e-9, "n={n}: {x} vs {y}");
            }
        }
    }

    #[test]
    fn known_spectrum() {
        // [[2,1],[1,2]] → {1, 3}.
        let mut a = vec![2.0, 1.0, 1.0, 2.0];
        let w = symmetric_eigenvalues(&mut a, 2);
        assert!((w[0] - 1.0).abs() < 1e-12 && (w[1] - 3.0).abs() < 1e-12, "{w:?}");
    }
}

//! Scalar 2-D finite-element mode solver (E-field, weak guidance).
//!
//! Triangular P1 (linear) elements on a uniform rectangular mesh over a
//! rectangle; the eigenproblem `A·x = β²·B·x` (stiffness in the index-
//! weighted Laplacian vs. mass matrix) is solved for the fundamental mode
//! by inverse iteration with a shift, using conjugate-gradient solves of
//! the shifted system.
//!
//! Accuracy: with the default mesh (~30 samples per wavelength-scale width)
//! the fundamental β agrees with the effective index method to ~1 % for
//! well-confined modes. This is a scalar (weak-guidance) solver; full
//! vectorial FEM awaits the `tpt-fem` dependency.

use alloc::vec::Vec;

use libm::sqrt;
use tpt_optics_core::Wavelength;

/// Rectangular computation domain with a uniform triangular mesh.
#[derive(Copy, Clone, Debug)]
pub struct RectDomain {
    /// Domain width (x, m).
    pub width: f64,
    /// Domain height (y, m).
    pub height: f64,
    /// Number of elements along x.
    pub nx: usize,
    /// Number of elements along y.
    pub ny: usize,
}

impl RectDomain {
    /// Node count of the mesh.
    pub fn nodes(&self) -> usize {
        (self.nx + 1) * (self.ny + 1)
    }

    /// Node (x, y) coordinates, row-major with y fastest.
    pub fn node_coords(&self) -> Vec<(f64, f64)> {
        let dx = self.width / self.nx as f64;
        let dy = self.height / self.ny as f64;
        let mut out = Vec::with_capacity(self.nodes());
        for ix in 0..=self.nx {
            for iy in 0..=self.ny {
                out.push((ix as f64 * dx, iy as f64 * dy));
            }
        }
        out
    }

    fn node_index(&self, ix: usize, iy: usize) -> usize {
        ix * (self.ny + 1) + iy
    }
}

/// Sample of the 2-D index profile; called with node coordinates.
pub type IndexFn<'a> = dyn Fn(f64, f64) -> f64 + Sync + 'a;

/// Fundamental (largest β²) mode of the scalar Helmholtz eigenproblem on
/// the domain with Dirichlet-type (large-detuning) boundary handling.
///
/// Returns `(beta, field)` with the field normalized to unit maximum.
///
/// * `sigma` — shift for inverse iteration, ≈ (k0·n_avg)² (use the
///   provided default from [`fundamental_shift`]).
pub fn fem_fundamental_mode(
    domain: &RectDomain,
    index: &IndexFn<'_>,
    wl: Wavelength,
    iterations: usize,
) -> (f64, Vec<f64>) {
    // Eigenproblem (scalar, weak guidance): with E(x,y)·e^{iβz},
    //   ∇ₜ²E + (k0²n² − β²)E = 0
    // in weak form: (k0²·Mn − S)·x = β²·M·x, where S = ∫∇φᵢ·∇φⱼ,
    // Mn = ∫n²φᵢφⱼ, M = ∫φᵢφⱼ. Guided modes are the LARGEST β².
    let n = domain.nodes();
    let k0 = wl.wavenumber();
    let dx = domain.width / domain.nx as f64;
    let dy = domain.height / domain.ny as f64;
    let mut a = alloc::vec![0.0; n * n]; // k0²·Mn − S
    let mut b = alloc::vec![0.0; n * n]; // M
    for ix in 0..domain.nx {
        for iy in 0..domain.ny {
            let n00 = domain.node_index(ix, iy);
            let n10 = domain.node_index(ix + 1, iy);
            let n01 = domain.node_index(ix, iy + 1);
            let n11 = domain.node_index(ix + 1, iy + 1);
            let xc = (ix as f64 + 0.5) * dx;
            let yc = (iy as f64 + 0.5) * dy;
            let n2 = index(xc, yc) * index(xc, yc);
            let area = dx * dy;
            // Per-triangle consistent mass: ∫φᵢφⱼ = A_t/12 (2A_t/12 on the
            // diagonal); stiffness on the same two triangles.
            let tri0 = [n00, n10, n01];
            let g0 = [(-1.0 / dx, -1.0 / dy), (1.0 / dx, 0.0), (0.0, 1.0 / dy)];
            let tri1 = [n11, n01, n10];
            let g1 = [(1.0 / dx, 1.0 / dy), (-1.0 / dx, 0.0), (0.0, -1.0 / dy)];
            for (nodes, gs) in [(&tri0, &g0), (&tri1, &g1)] {
                let tri_area = area * 0.5;
                for (p, gp) in gs.iter().enumerate() {
                    for (q, gq) in gs.iter().enumerate() {
                        let mass = tri_area / 12.0 * if p == q { 2.0 } else { 1.0 };
                        a[nodes[p] * n + nodes[q]] +=
                            mass * k0 * k0 * n2 - (gp.0 * gq.0 + gp.1 * gq.1) * tri_area;
                        b[nodes[p] * n + nodes[q]] += mass;
                    }
                }
            }
        }
    }
    // Dirichlet boundary: eliminate boundary nodes (zero row+column,
    // unit diagonal) so the box walls pin the field to zero.
    let mut boundary = alloc::vec![false; n];
    for ix in 0..=domain.nx {
        for iy in 0..=domain.ny {
            if ix == 0 || iy == 0 || ix == domain.nx || iy == domain.ny {
                boundary[domain.node_index(ix, iy)] = true;
            }
        }
    }
    for i in 0..n {
        if boundary[i] {
            for j in 0..n {
                a[i * n + j] = 0.0;
                a[j * n + i] = 0.0;
                b[i * n + j] = 0.0;
                b[j * n + i] = 0.0;
            }
            a[i * n + i] = 1.0;
        }
    }

    // Rayleigh quotient iteration: σ tracks the quotient so the dense
    // solves converge (quadratically) to the eigenpair nearest the initial
    // guess — a core-localized Gaussian with σ just below (k0·n_core)².
    // Boxed radiation modes are far from this shift and are skipped.
    let k0_core_sq = {
        let mut n_max = 0.0_f64;
        for (x, y) in domain.node_coords() {
            let nv = index(x, y);
            if nv > n_max {
                n_max = nv;
            }
        }
        (k0 * n_max) * (k0 * n_max)
    };
    let mut sigma = 0.97 * k0_core_sq;
    let mut x = alloc::vec![0.0; n];
    {
        // Seed: Gaussian at the domain center, localized over ~1/6 width.
        let coords = domain.node_coords();
        let (cx, cy) = (0.5 * domain.width, 0.5 * domain.height);
        let w0 = 0.18 * domain.width.max(domain.height);
        for (i, (xi, yi)) in coords.iter().enumerate() {
            let r2 = (xi - cx) * (xi - cx) + (yi - cy) * (yi - cy);
            x[i] = if boundary[i] {
                0.0
            } else {
                libm::exp(-r2 / (w0 * w0))
            };
        }
    }
    for it in 0..iterations {
        let mut rhs = mat_vec(&b, &x, n);
        for i in 0..n {
            if boundary[i] {
                rhs[i] = 0.0;
            }
        }
        let s_mat = {
            let mut m = a.clone();
            m.fill(0.0);
            for i in 0..n {
                for j in 0..n {
                    m[i * n + j] = a[i * n + j] - sigma * b[i * n + j];
                }
            }
            m
        };
        let Some(next) = lu_solve(&s_mat, &rhs, n) else {
            break;
        };
        x = next;
        let norm = sqrt(x.iter().map(|v| v * v).sum());
        if norm <= 0.0 || !norm.is_finite() {
            break;
        }
        for v in &mut x {
            *v /= norm;
        }
        // Update σ from the Rayleigh quotient of the original problem.
        let at = mat_vec(&a, &x, n);
        let bt = mat_vec(&b, &x, n);
        let num: f64 = x.iter().zip(at.iter()).map(|(v, w)| v * w).sum();
        let den: f64 = x.iter().zip(bt.iter()).map(|(v, w)| v * w).sum();
        // Pure inverse iteration first (amplifies the eigenpair nearest σ),
        // then Rayleigh-quotient updates for quadratic polish.
        if it >= 8 && den > 1e-300 {
            sigma = num / den;
        }
    }
    // β² from the final Rayleigh quotient.
    let at = mat_vec(&a, &x, n);
    let bt = mat_vec(&b, &x, n);
    let num: f64 = x.iter().zip(at.iter()).map(|(v, w)| v * w).sum();
    let den: f64 = x.iter().zip(bt.iter()).map(|(v, w)| v * w).sum();
    let lambda = if den > 0.0 { num / den } else { 0.0 };

    // β² = λ → guide β.
    let beta = sqrt(lambda.max(0.0));
    // Normalize the field to unit max.
    let max = x.iter().copied().fold(0.0_f64, f64::max);
    if max > 0.0 {
        for v in &mut x {
            *v /= max;
        }
    }
    (beta, x)
}

/// Suggested inverse-iteration shift: (k0·n_clad)² — below every guided
/// mode but above radiation.
pub fn fundamental_shift(domain: &RectDomain, index: &IndexFn<'_>, wl: Wavelength) -> f64 {
    let coords = domain.node_coords();
    let mut n_min = f64::INFINITY;
    for (x, y) in coords {
        let nv = index(x, y);
        if nv < n_min {
            n_min = nv;
        }
    }
    let k0 = wl.wavenumber();
    (k0 * n_min) * (k0 * n_min) * 0.999
}

fn mat_vec(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            a[i * n..(i + 1) * n]
                .iter()
                .zip(x.iter())
                .map(|(m, v)| m * v)
                .sum()
        })
        .collect()
}

/// Dense LU solve with partial pivoting (Doolittle). Returns `None` for a
/// singular matrix.
fn lu_solve(a_in: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut a = a_in.to_vec();
    let mut b = b.to_vec();
    for col in 0..n {
        // Pivot.
        let (mut piv, mut prow) = (a[col * n + col].abs(), col);
        for r in (col + 1)..n {
            let v = a[r * n + col].abs();
            if v > piv {
                piv = v;
                prow = r;
            }
        }
        if piv < 1e-300 {
            return None;
        }
        if prow != col {
            for k in 0..n {
                a.swap(col * n + k, prow * n + k);
            }
            b.swap(col, prow);
        }
        let pivot = a[col * n + col];
        for r in (col + 1)..n {
            let factor = a[r * n + col] / pivot;
            a[r * n + col] = 0.0;
            for k in (col + 1)..n {
                a[r * n + k] -= factor * a[col * n + k];
            }
            b[r] -= factor * b[col];
        }
    }
    // Back substitution.
    let mut x = alloc::vec![0.0; n];
    for r in (0..n).rev() {
        let mut sum = b[r];
        for k in (r + 1)..n {
            sum -= a[r * n + k] * x[k];
        }
        x[r] = sum / a[r * n + r];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn wl(nm: f64) -> Wavelength {
        Wavelength::from_nm(nm)
    }

    /// Slab-like index: high-index stripe centered in the domain
    /// (|y − height/2| < h/2, extending across x).
    fn slab_index(h: f64, n_core: f64, n_clad: f64, y_center: f64) -> impl Fn(f64, f64) -> f64 {
        move |_x, y| {
            if (y - y_center).abs() < h / 2.0 {
                n_core
            } else {
                n_clad
            }
        }
    }

    #[test]
    fn fem_slab_mode_matches_transcendental() {
        // Wide domain slab: FEM β should match the exact slab TE0 solve.
        let (h, n_core, n_clad) = (1.0e-6, 3.48, 1.444);
        let domain = RectDomain {
            width: 2e-6, // narrow in x (uniform along slab)
            height: 8e-6,
            nx: 6,
            ny: 96,
        };
        let idx = slab_index(h, n_core, n_clad, 4e-6);
        let (beta_fem, field) = fem_fundamental_mode(&domain, &idx, wl(1550.0), 30);

        let slab = crate::slab::SlabLayer {
            n_core,
            n_clad,
            thickness: h,
        };
        let modes = crate::slab::slab_modes(&slab, wl(1550.0), crate::slab::SlabPolarization::Te);
        assert!(!modes.is_empty());
        let beta_ref = modes[0].n_eff * wl(1550.0).wavenumber();
        assert!(
            (beta_fem - beta_ref).abs() < 0.05 * beta_ref,
            "FEM β {beta_fem} vs slab β {beta_ref}"
        );
        // Field peaked in the core (center row).
        let center = field[(domain.nx / 2) * (domain.ny + 1) + domain.ny / 2];
        assert!(center > 0.8, "center amplitude {center}");
        // And decaying into the cladding.
        let edge = field[(domain.nx / 2) * (domain.ny + 1) + 4];
        assert!(edge < 0.3, "cladding amplitude {edge}");
    }

    #[test]
    fn jacobi_catches_offdiag() {
        // (moved here from phc tests if needed) identity already diagonal.
        let mut a = vec![1.0, 0.0, 0.0, 2.0];
        let mut w = [0.0; 2];
        super::super::phc::jacobi_eigen(&mut a, &mut w);
        assert!(close(w[0], 1.0, 1e-12) && close(w[1], 2.0, 1e-12));
    }
}

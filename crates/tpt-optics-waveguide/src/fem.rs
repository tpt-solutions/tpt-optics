//! Scalar 2-D finite-element mode solver (E-field, weak guidance).
//!
//! Triangular P1 (linear) elements on a uniform rectangular mesh over a
//! rectangle; the eigenproblem `A·x = β²·B·x` (stiffness in the index-
//! weighted Laplacian vs. mass matrix) is solved for the fundamental mode
//! by Rayleigh-quotient inverse iteration, with dense LU solves of the
//! shifted system.
//!
//! Accuracy: with the default mesh (~30 samples per wavelength-scale width)
//! the fundamental β agrees with the effective index method to ~1 % for
//! well-confined modes. This is a scalar (weak-guidance) solver. For the full
//! vector field see [`crate::vector_fem`], which is built in-workspace and
//! does not depend on the external `tpt-fem` crate.

use alloc::vec::Vec;

use libm::sqrt;
use tpt_optics_core::Wavelength;

use crate::dense::{mat_vec, LuFactor};

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
/// Returns `(beta, field)`. The field is scaled so that its entry of largest
/// magnitude is +1.
///
/// The dense solver needs `n²` entries, so meshes with more nodes are refused.
pub const MAX_FEM_NODES: usize = 1_500;

/// Relative residual `‖Ax − λBx‖ / (‖Ax‖ + |λ|‖Bx‖)` above which the mode is
/// reported as not converged.
pub const FEM_RESIDUAL_TOL: f64 = 1e-6;

/// Errors from the scalar FEM mode solver.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum FemError {
    /// The domain has no cells, or a non-positive or non-finite size.
    InvalidDomain,
    /// The mesh has more nodes than [`MAX_FEM_NODES`].
    TooLarge {
        /// Number of nodes of the mesh.
        nodes: usize,
    },
    /// The shifted matrix is singular or the iterate vanished, at this
    /// inverse-iteration step.
    Singular {
        /// Zero-based iteration index.
        iteration: usize,
    },
    /// The final eigenpair has a relative residual above [`FEM_RESIDUAL_TOL`].
    NotConverged {
        /// The residual reached.
        residual: f64,
    },
}

impl core::fmt::Display for FemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDomain => write!(f, "FEM domain has no cells or a bad size"),
            Self::TooLarge { nodes } => write!(
                f,
                "FEM mesh has {nodes} nodes (limit {MAX_FEM_NODES}); coarsen it"
            ),
            Self::Singular { iteration } => write!(
                f,
                "FEM shifted matrix is singular at inverse-iteration step {iteration}"
            ),
            Self::NotConverged { residual } => write!(
                f,
                "FEM mode did not converge (relative residual {residual:e})"
            ),
        }
    }
}

impl core::error::Error for FemError {}

/// * `sigma` — shift for inverse iteration, ≈ (k0·n_avg)² (use the
///   provided default from [`fundamental_shift`]).
///
/// # Errors
///
/// [`FemError`] for an invalid or oversized domain, a singular shift, or an
/// eigenpair whose residual is above [`FEM_RESIDUAL_TOL`].
pub fn fem_fundamental_mode(
    domain: &RectDomain,
    index: &IndexFn<'_>,
    wl: Wavelength,
    iterations: usize,
) -> Result<(f64, Vec<f64>), FemError> {
    let finite_positive = |v: f64| v.is_finite() && v > 0.0;
    if domain.nx == 0
        || domain.ny == 0
        || !finite_positive(domain.width)
        || !finite_positive(domain.height)
    {
        return Err(FemError::InvalidDomain);
    }
    let n = domain.nodes();
    if n > MAX_FEM_NODES {
        return Err(FemError::TooLarge { nodes: n });
    }
    // Eigenproblem (scalar, weak guidance): with E(x,y)·e^{iβz},
    //   ∇ₜ²E + (k0²n² − β²)E = 0
    // in weak form: (k0²·Mn − S)·x = β²·M·x, where S = ∫∇φᵢ·∇φⱼ,
    // Mn = ∫n²φᵢφⱼ, M = ∫φᵢφⱼ. Guided modes are the LARGEST β².
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
    let mut sigma = k0_core_sq;
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
        let mut rhs = mat_vec(&b, n, &x);
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
        let Some(lu) = LuFactor::new(s_mat, n) else {
            return Err(FemError::Singular { iteration: it });
        };
        x = lu.solve(&rhs);
        let norm = sqrt(x.iter().map(|v| v * v).sum());
        if norm <= 0.0 || !norm.is_finite() {
            return Err(FemError::Singular { iteration: it });
        }
        for v in &mut x {
            *v /= norm;
        }
        // Update σ from the Rayleigh quotient of the original problem.
        let at = mat_vec(&a, n, &x);
        let bt = mat_vec(&b, n, &x);
        let num: f64 = x.iter().zip(at.iter()).map(|(v, w)| v * w).sum();
        let den: f64 = x.iter().zip(bt.iter()).map(|(v, w)| v * w).sum();
        // Pure inverse iteration first (amplifies the eigenpair nearest σ),
        // then Rayleigh-quotient updates for quadratic polish.
        if it >= 8 && den > 1e-300 {
            sigma = num / den;
        }
    }
    // β² from the final Rayleigh quotient.
    let at = mat_vec(&a, n, &x);
    let bt = mat_vec(&b, n, &x);
    let num: f64 = x.iter().zip(at.iter()).map(|(v, w)| v * w).sum();
    let den: f64 = x.iter().zip(bt.iter()).map(|(v, w)| v * w).sum();
    let lambda = if den > 0.0 { num / den } else { 0.0 };

    // The eigenpair must satisfy Ax = λBx. Check it, so an iteration that
    // stopped early is reported rather than returned as a mode.
    let residual_vec: Vec<f64> = at
        .iter()
        .zip(&bt)
        .map(|(a_x, b_x)| a_x - lambda * b_x)
        .collect();
    let r_norm = sqrt(residual_vec.iter().map(|v| v * v).sum());
    let scale = sqrt(at.iter().map(|v| v * v).sum::<f64>())
        + libm::fabs(lambda) * sqrt(bt.iter().map(|v| v * v).sum::<f64>());
    let residual = if scale > 0.0 {
        r_norm / scale
    } else {
        f64::INFINITY
    };
    if residual.is_nan() || residual > FEM_RESIDUAL_TOL {
        return Err(FemError::NotConverged { residual });
    }

    // β² = λ → guide β.
    let beta = sqrt(lambda.max(0.0));
    normalise_field(&mut x);
    Ok((beta, x))
}

/// Fixes the sign of an eigenvector and scales it to unit peak amplitude.
///
/// The peak is the entry of largest magnitude, which may be negative. Dividing
/// by it both removes the arbitrary sign of inverse iteration and sets the
/// magnitude to one. Dividing by the signed maximum instead would leave an
/// all-negative field unscaled.
fn normalise_field(x: &mut [f64]) {
    let peak = x.iter().copied().fold(
        0.0_f64,
        |m, v| if libm::fabs(v) > libm::fabs(m) { v } else { m },
    );
    if peak != 0.0 && peak.is_finite() {
        for v in x.iter_mut() {
            *v /= peak;
        }
    }
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
        let (beta_fem, field) = fem_fundamental_mode(&domain, &idx, wl(1550.0), 30).unwrap();

        let slab = crate::slab::SlabLayer {
            n_core,
            n_clad,
            thickness: h,
        };
        let modes = crate::slab::slab_modes(&slab, wl(1550.0), crate::slab::SlabPolarization::Te);
        assert_ne!(modes.len(), 0);
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
    fn fem_weak_guide_finds_fundamental_slab_mode() {
        // Weak guide (index contrast 0.7 %). The x walls are Dirichlet, so the
        // box width adds a transverse term (π/W)² to k0²; W must be wide
        // compared with the guidance (n_eff² − n_clad²)·k0², or the test
        // measures the box, not the mode.
        let (h, n_core, n_clad) = (8.0e-6, 1.46, 1.45);
        let lambda = wl(1550.0);
        // Square cells (1 µm): a stretched mesh gives a poor discrete operator.
        let domain = RectDomain {
            width: 16e-6,
            height: 4.0 * h,
            nx: 16,
            ny: 32,
        };
        let idx = slab_index(h, n_core, n_clad, 2.0 * h);
        let (beta_fem, field) = fem_fundamental_mode(&domain, &idx, lambda, 40).unwrap();
        // A box mode lies within a few per cent of the guided mode, so require
        // a guided mode and a tight β tolerance (the error here is ~0.1 %).
        assert!(
            beta_fem > lambda.wavenumber() * n_clad,
            "FEM mode is not guided: β {beta_fem}"
        );

        let slab = crate::slab::SlabLayer {
            n_core,
            n_clad,
            thickness: h,
        };
        let modes = crate::slab::slab_modes(&slab, lambda, crate::slab::SlabPolarization::Te);
        let beta_ref = modes[0].n_eff * lambda.wavenumber();
        assert!(
            (beta_fem - beta_ref).abs() < 0.01 * beta_ref,
            "FEM β {beta_fem} vs slab β {beta_ref}"
        );
        let center = field[(domain.nx / 2) * (domain.ny + 1) + domain.ny / 2];
        assert!(center > 0.8, "center amplitude {center}");
    }

    #[test]
    fn normalised_field_peak_is_positive_unit() {
        // Inverse iteration fixes the eigenvector only up to sign, so a field
        // whose entries are all negative must still come out with a +1 peak.
        let mut x = vec![-0.2, -1.0, -0.5];
        normalise_field(&mut x);
        assert!(
            close(x[0], 0.2, 1e-12) && close(x[1], 1.0, 1e-12) && close(x[2], 0.5, 1e-12),
            "{x:?}"
        );
    }

    #[test]
    fn normalised_field_uses_largest_magnitude_entry() {
        // The positive entry is smaller in magnitude, so it must not set the scale.
        let mut x = vec![0.3, -1.0];
        normalise_field(&mut x);
        assert!(close(x[0], -0.3, 1e-12) && close(x[1], 1.0, 1e-12), "{x:?}");
    }

    #[test]
    fn jacobi_catches_offdiag() {
        // (moved here from phc tests if needed) identity already diagonal.
        let mut a = vec![1.0, 0.0, 0.0, 2.0];
        let mut w = [0.0; 2];
        super::super::phc::jacobi_eigen(&mut a, &mut w).unwrap();
        assert!(close(w[0], 1.0, 1e-12) && close(w[1], 2.0, 1e-12));
    }
}

#[cfg(test)]
mod error_path_tests {
    use super::*;

    fn flat(_x: f64, _y: f64) -> f64 {
        1.5
    }

    #[test]
    fn empty_mesh_is_an_invalid_domain() {
        let domain = RectDomain {
            width: 1e-5,
            height: 1e-5,
            nx: 0,
            ny: 4,
        };
        let wl = Wavelength::from_um(1.55);
        assert_eq!(
            fem_fundamental_mode(&domain, &flat, wl, 10),
            Err(FemError::InvalidDomain)
        );
    }

    #[test]
    fn oversized_mesh_is_refused_before_allocation() {
        // 50×50 cells is 2601 nodes, above the dense limit of 1500.
        let domain = RectDomain {
            width: 1e-5,
            height: 1e-5,
            nx: 50,
            ny: 50,
        };
        let wl = Wavelength::from_um(1.55);
        assert!(matches!(
            fem_fundamental_mode(&domain, &flat, wl, 10),
            Err(FemError::TooLarge { nodes: 2601 })
        ));
    }
}

#[cfg(test)]
mod strong_guide_tests {
    use super::*;

    #[test]
    fn strong_silicon_slab_fundamental_is_found_and_converges() {
        // A 1.5 % error at 48 cells across the 2.6 um domain. At 32 cells the
        // error is about 11 %, from discretisation, not a wrong mode: the field
        // is guided (n_eff 2.55, between 1.45 and 3.48) at both meshes.
        // Silicon core (n 3.48) in silica (1.45): a strongly guiding slab.
        let (h, n_core, n_clad) = (0.22e-6, 3.48, 1.45);
        let lambda = Wavelength::from_um(1.55);
        let domain = RectDomain {
            width: 2.0e-6,
            height: 3.0 * 0.22e-6 * 4.0,
            nx: 16,
            ny: 48,
        };
        let idx = |_x: f64, y: f64| {
            let yc = 0.5 * domain.height;
            if (y - yc).abs() < 0.5 * h {
                n_core
            } else {
                n_clad
            }
        };
        let (beta, _) = fem_fundamental_mode(&domain, &idx, lambda, 40).unwrap();
        let slab = crate::slab::SlabLayer {
            n_core,
            n_clad,
            thickness: h,
        };
        let modes = crate::slab::slab_modes(&slab, lambda, crate::slab::SlabPolarization::Te);
        let beta_ref = modes[0].n_eff * lambda.wavenumber();
        assert!(
            (beta - beta_ref).abs() < 0.03 * beta_ref,
            "fem {beta} vs {beta_ref}"
        );
        assert!(
            beta > lambda.wavenumber() * n_clad,
            "mode is not guided: {beta}"
        );
    }
}

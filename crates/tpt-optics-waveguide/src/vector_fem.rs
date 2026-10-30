//! Full-vector mode solver: hybrid edge–nodal finite elements on a
//! triangulated rectangle.
//!
//! The transverse field `E_t` uses first-order Nédélec (edge) elements and
//! the longitudinal field `E_z` uses P1 nodal elements. Fields are
//! `E(x, y)·e^{-iβz}` with `ε(x, y) = n²`, and the PEC walls set
//! `n × E = 0` on the boundary.
//!
//! Writing `ψ = i·E_z` makes every operator real. The weak form is
//!
//! `(P − k₀²Q)·e + β²·M·e − β·G·ψ = 0`,  `−β·Gᵀe + (S − k₀²T)·ψ = 0`,
//!
//! with edge curl–curl `P`, ε-weighted edge mass `Q`, edge mass `M`, the
//! edge–node coupling `G`, nodal stiffness `S` and ε-weighted nodal mass `T`.
//! Eliminating `ψ = β·K⁻¹Gᵀe` with `K = S − k₀²T` gives the generalised
//! eigenproblem `A·e = β²·B·e` with `A = P − k₀²Q` and `B = G·K⁻¹·Gᵀ − M`.
//!
//! The eigenvalue closest to a shift is found by inverse iteration on
//! `(A − σB)`. The solver is dense, so keep the mesh modest (a few hundred
//! edges) to stay practical.
//!
//! Accuracy limits: piecewise-constant ε per triangle (sampled at the
//! centroid), first-order elements, and a dense solve whose cost grows as
//! `n³` in the edge count. Edge elements avoid the spurious modes of nodal
//! vector formulations, but the solver has not been checked on tensor-weighted
//! or high-contrast interface meshes.

use alloc::collections::BTreeMap;
use alloc::vec;
use alloc::vec::Vec;

use libm::{fabs, sqrt};
use tpt_optics_core::Wavelength;

use crate::dense::{mat_vec, LuFactor};
use crate::fem::{IndexFn, RectDomain};

/// Largest number of edges, and of nodes, that the dense solver accepts. The
/// matrices are `n²` entries each, so a finer mesh is refused before any
/// allocation rather than exhausting memory.
pub const MAX_DENSE_UNKNOWNS: usize = 1_500;

/// Errors from the vector mode solver.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum VectorFemError {
    /// An input is non-finite or out of range. The message names it.
    InvalidParameter(&'static str),
    /// The mesh needs more edges or nodes than [`MAX_DENSE_UNKNOWNS`].
    TooLarge {
        /// Number of edges the mesh would have.
        edges: usize,
        /// Number of nodes the mesh would have.
        nodes: usize,
    },
    /// A dense factorisation hit a zero or non-finite pivot. The shift is
    /// probably exactly an eigenvalue, or the matrix is degenerate.
    Singular,
    /// The converged `β²` is not positive, so the mode is evanescent (below
    /// cutoff) and has no real effective index.
    NotGuided {
        /// The converged `β²` (1/m²), not positive.
        beta_sq: f64,
    },
    /// The iteration produced a non-finite value.
    NonFinite,
}

impl core::fmt::Display for VectorFemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidParameter(msg) => write!(f, "invalid vector-FEM parameter: {msg}"),
            Self::TooLarge { edges, nodes } => write!(
                f,
                "mesh has {edges} edges and {nodes} nodes (limit {MAX_DENSE_UNKNOWNS}); coarsen it"
            ),
            Self::Singular => write!(f, "dense factorisation is singular"),
            Self::NotGuided { beta_sq } => {
                write!(
                    f,
                    "mode is evanescent (β² = {beta_sq}); no real effective index"
                )
            }
            Self::NonFinite => write!(f, "vector-FEM iteration produced a non-finite value"),
        }
    }
}

impl core::error::Error for VectorFemError {}

/// A converged vector mode.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct VectorMode {
    /// Effective index `√β² / k₀`.
    pub n_eff: f64,
    /// Relative residual `‖A·x − β²·B·x‖ / ‖A·x‖` of the returned mode.
    pub residual: f64,
    /// True when the Rayleigh quotient changed by less than `1e-10` relative
    /// over the last inverse-iteration step.
    pub converged: bool,
}

/// Effective index of the mode closest to `(k₀·n_guess)²` in the full-vector
/// hybrid edge–nodal formulation on `domain` with PEC walls.
///
/// `iterations` inverse-iteration steps are used; 30–60 is usually enough
/// when `n_guess` is within a few percent of the answer. Check
/// [`VectorMode::converged`] and [`VectorMode::residual`] before relying on
/// the result.
///
/// # Errors
/// See [`VectorFemError`]. An evanescent solution is an error, not an index.
pub fn vector_mode_neff(
    domain: &RectDomain,
    index: &IndexFn<'_>,
    wl: Wavelength,
    n_guess: f64,
    iterations: usize,
) -> Result<VectorMode, VectorFemError> {
    let bad = VectorFemError::InvalidParameter;
    if domain.nx == 0 || domain.ny == 0 {
        return Err(bad("mesh needs at least one cell in each direction"));
    }
    if !(domain.width.is_finite() && domain.height.is_finite())
        || domain.width <= 0.0
        || domain.height <= 0.0
    {
        return Err(bad("domain width and height must be finite and positive"));
    }
    let k0 = wl.wavenumber();
    if !(k0.is_finite() && k0 > 0.0) {
        return Err(bad("wavelength must be finite and positive"));
    }
    if !(n_guess.is_finite() && n_guess > 0.0) {
        return Err(bad("n_guess must be finite and positive"));
    }
    if iterations == 0 {
        return Err(bad("at least one inverse-iteration step is needed"));
    }
    let (nx, ny) = (domain.nx, domain.ny);
    let n_edges_expected = nx * (ny + 1) + (nx + 1) * ny + nx * ny;
    let n_nodes_expected = (nx + 1) * (ny + 1);
    if n_edges_expected > MAX_DENSE_UNKNOWNS || n_nodes_expected > MAX_DENSE_UNKNOWNS {
        return Err(VectorFemError::TooLarge {
            edges: n_edges_expected,
            nodes: n_nodes_expected,
        });
    }
    let k0sq = k0 * k0;
    let nodes = domain.node_coords();
    let n_nodes = nodes.len();
    let node_at = |ix: usize, iy: usize| ix * (ny + 1) + iy;
    let node_ij = |n: usize| (n / (ny + 1), n % (ny + 1));

    // Triangles and global edges (keyed by sorted node pair).
    let mut tris: Vec<[usize; 3]> = Vec::with_capacity(2 * nx * ny);
    for ix in 0..nx {
        for iy in 0..ny {
            let n00 = node_at(ix, iy);
            let n10 = node_at(ix + 1, iy);
            let n01 = node_at(ix, iy + 1);
            let n11 = node_at(ix + 1, iy + 1);
            tris.push([n00, n10, n01]);
            tris.push([n11, n01, n10]);
        }
    }
    let mut edge_id: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut edge_nodes: Vec<(usize, usize)> = Vec::new();
    // Per triangle: global edge id and orientation sign for local edges
    // (0,1), (1,2), (2,0).
    let local_pairs = [(0usize, 1usize), (1, 2), (2, 0)];
    let mut tri_edges: Vec<[(usize, f64); 3]> = Vec::with_capacity(tris.len());
    for t in &tris {
        let mut row = [(0usize, 0.0f64); 3];
        for (slot, &(a, b)) in local_pairs.iter().enumerate() {
            let (p, q) = (t[a], t[b]);
            let key = if p < q { (p, q) } else { (q, p) };
            let next = edge_nodes.len();
            let id = *edge_id.entry(key).or_insert_with(|| {
                edge_nodes.push(key);
                next
            });
            let sign = if p < q { 1.0 } else { -1.0 };
            row[slot] = (id, sign);
        }
        tri_edges.push(row);
    }
    let n_edges = edge_nodes.len();

    // Assemble dense global matrices (row-major).
    let mut p_mat = vec![0.0; n_edges * n_edges];
    let mut m_mat = vec![0.0; n_edges * n_edges];
    let mut q_mat = vec![0.0; n_edges * n_edges];
    let mut g_mat = vec![0.0; n_edges * n_nodes]; // G[edge][node]
    let mut s_mat = vec![0.0; n_nodes * n_nodes];
    let mut t_mat = vec![0.0; n_nodes * n_nodes];

    for (ti, t) in tris.iter().enumerate() {
        let (grad, area) = triangle_gradients([nodes[t[0]], nodes[t[1]], nodes[t[2]]]);
        let cx = (nodes[t[0]].0 + nodes[t[1]].0 + nodes[t[2]].0) / 3.0;
        let cy = (nodes[t[0]].1 + nodes[t[1]].1 + nodes[t[2]].1) / 3.0;
        let eps = {
            let n = index(cx, cy);
            n * n
        };
        if !(eps.is_finite() && eps > 0.0) {
            return Err(bad("index must be positive and finite at every triangle"));
        }
        // Local edge basis N_(i,j) = λ_i ∇λ_j − λ_j ∇λ_i for the three local
        // edges, dotted against each other.
        let mut le = [[0.0f64; 3]; 3];
        let mut le_q = [[0.0f64; 3]; 3];
        let mut le_p = [[0.0f64; 3]; 3];
        for (a, &(i, j)) in local_pairs.iter().enumerate() {
            for (b, &(k, l)) in local_pairs.iter().enumerate() {
                le[a][b] = edge_dot(i, j, k, l, &grad, area);
                le_q[a][b] = eps * le[a][b];
                le_p[a][b] = 4.0 * area * cross(grad[i], grad[j]) * cross(grad[k], grad[l]);
            }
        }
        let ed = tri_edges[ti];
        for a in 0..3 {
            for b in 0..3 {
                let (ia, sa) = ed[a];
                let (ib, sb) = ed[b];
                let s = sa * sb;
                m_mat[ia * n_edges + ib] += s * le[a][b];
                q_mat[ia * n_edges + ib] += s * le_q[a][b];
                p_mat[ia * n_edges + ib] += s * le_p[a][b];
            }
            // Edge–node coupling: ∫ N_(i,j)·∇λ_c = (A/3)(∇λ_j − ∇λ_i)·∇λ_c.
            let (i, j) = local_pairs[a];
            let (ia, sa) = ed[a];
            for c in 0..3 {
                let val = area / 3.0 * (dot(grad[j], grad[c]) - dot(grad[i], grad[c]));
                let gn = t[c];
                g_mat[ia * n_nodes + gn] += sa * val;
            }
        }
        for c in 0..3 {
            for d in 0..3 {
                let gn_c = t[c];
                let gn_d = t[d];
                s_mat[gn_c * n_nodes + gn_d] += area * dot(grad[c], grad[d]);
                let mass = eps * area / 12.0 * if c == d { 2.0 } else { 1.0 };
                t_mat[gn_c * n_nodes + gn_d] += mass;
            }
        }
    }

    // PEC: drop boundary edges and nodes.
    let on_boundary_node = |n: usize| {
        let (ix, iy) = node_ij(n);
        ix == 0 || iy == 0 || ix == nx || iy == ny
    };
    let interior_edge: Vec<usize> = (0..n_edges)
        .filter(|&e| {
            let (p, q) = edge_nodes[e];
            let (ip, jp) = node_ij(p);
            let (iq, jq) = node_ij(q);
            let side_x = ip == iq && (ip == 0 || ip == nx);
            let side_y = jp == jq && (jp == 0 || jp == ny);
            !(side_x || side_y)
        })
        .collect();
    let interior_node: Vec<usize> = (0..n_nodes).filter(|&n| !on_boundary_node(n)).collect();
    let ne = interior_edge.len();
    let nn = interior_node.len();

    let sub = |full: &[f64], rows: &[usize], cols: &[usize], stride: usize| -> Vec<f64> {
        let mut out = vec![0.0; rows.len() * cols.len()];
        for (r, &ri) in rows.iter().enumerate() {
            for (c, &ci) in cols.iter().enumerate() {
                out[r * cols.len() + c] = full[ri * stride + ci];
            }
        }
        out
    };
    let p_int = sub(&p_mat, &interior_edge, &interior_edge, n_edges);
    let m_int = sub(&m_mat, &interior_edge, &interior_edge, n_edges);
    let q_int = sub(&q_mat, &interior_edge, &interior_edge, n_edges);
    let s_int = sub(&s_mat, &interior_node, &interior_node, n_nodes);
    let t_int = sub(&t_mat, &interior_node, &interior_node, n_nodes);
    // G_int[e][n] and its transpose.
    let g_int = sub(&g_mat, &interior_edge, &interior_node, n_nodes);

    // K = S − k₀²T; W = K⁻¹·Gᵀ (nn × ne).
    let mut k_mat = vec![0.0; nn * nn];
    for i in 0..nn * nn {
        k_mat[i] = s_int[i] - k0sq * t_int[i];
    }
    let mut gt = vec![0.0; nn * ne];
    for e in 0..ne {
        for n in 0..nn {
            gt[n * ne + e] = g_int[e * nn + n];
        }
    }
    let w = solve_dense(&k_mat, nn, &gt, ne)?;

    // A = P − k₀²Q;  B = G·W − M.
    let mut a_mat = vec![0.0; ne * ne];
    let mut b_mat = vec![0.0; ne * ne];
    for i in 0..ne * ne {
        a_mat[i] = p_int[i] - k0sq * q_int[i];
    }
    for i in 0..ne {
        for j in 0..ne {
            let mut s = 0.0;
            for n in 0..nn {
                s += g_int[i * nn + n] * w[n * ne + j];
            }
            b_mat[i * ne + j] = s - m_int[i * ne + j];
        }
    }

    // Inverse iteration on (A − σB)⁻¹B.
    let sigma = (k0 * n_guess) * (k0 * n_guess);
    let mut shifted = vec![0.0; ne * ne];
    for i in 0..ne * ne {
        shifted[i] = a_mat[i] - sigma * b_mat[i];
    }
    let lu = LuFactor::new(shifted, ne).ok_or(VectorFemError::Singular)?;
    let mut x = vec![1.0; ne];
    let mut prev = f64::NAN;
    let mut converged = false;
    for _ in 0..iterations {
        let y = mat_vec(&b_mat, ne, &x);
        let z = lu.solve(&y);
        let norm = sqrt(z.iter().map(|v| v * v).sum::<f64>());
        if !(norm.is_finite() && norm > 0.0) {
            return Err(VectorFemError::NonFinite);
        }
        x = z.iter().map(|v| v / norm).collect();
        let lambda = rayleigh_quotient(&a_mat, &b_mat, &x, ne)?;
        converged = libm::fabs(lambda - prev) <= 1e-10 * libm::fabs(lambda);
        prev = lambda;
    }
    let ax = mat_vec(&a_mat, ne, &x);
    let bx = mat_vec(&b_mat, ne, &x);
    let beta_sq = rayleigh_quotient(&a_mat, &b_mat, &x, ne)?;
    if beta_sq <= 0.0 {
        return Err(VectorFemError::NotGuided { beta_sq });
    }
    let residual_sq: f64 = ax
        .iter()
        .zip(&bx)
        .map(|(p, q)| (p - beta_sq * q) * (p - beta_sq * q))
        .sum();
    let ax_norm = sqrt(ax.iter().map(|p| p * p).sum::<f64>());
    Ok(VectorMode {
        n_eff: sqrt(beta_sq) / k0,
        residual: sqrt(residual_sq) / ax_norm,
        converged,
    })
}

/// Rayleigh quotient `xᵀA x / xᵀB x`, which estimates the eigenvalue β².
fn rayleigh_quotient(a: &[f64], b: &[f64], x: &[f64], n: usize) -> Result<f64, VectorFemError> {
    let ax = mat_vec(a, n, x);
    let bx = mat_vec(b, n, x);
    let num: f64 = x.iter().zip(&ax).map(|(p, q)| p * q).sum();
    let den: f64 = x.iter().zip(&bx).map(|(p, q)| p * q).sum();
    let q = num / den;
    if q.is_finite() {
        Ok(q)
    } else {
        Err(VectorFemError::NonFinite)
    }
}

/// Gradients `∇λ_i` and area of a triangle with vertices `v`.
fn triangle_gradients(v: [(f64, f64); 3]) -> ([(f64, f64); 3], f64) {
    let (x0, y0) = v[0];
    let (x1, y1) = v[1];
    let (x2, y2) = v[2];
    let det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    let grad = [
        ((y1 - y2) / det, (x2 - x1) / det),
        ((y2 - y0) / det, (x0 - x2) / det),
        ((y0 - y1) / det, (x1 - x0) / det),
    ];
    (grad, fabs(det) * 0.5)
}

/// `∫ N_(i,j)·N_(k,l) dA` for the local edge basis on one triangle.
fn edge_dot(i: usize, j: usize, k: usize, l: usize, grad: &[(f64, f64); 3], area: f64) -> f64 {
    // N_(i,j)·N_(k,l) = λ_iλ_k(∇λ_j·∇λ_l) − λ_iλ_l(∇λ_j·∇λ_k)
    //                 − λ_jλ_k(∇λ_i·∇λ_l) + λ_jλ_l(∇λ_i·∇λ_k).
    let lam = |p: usize, q: usize| area * if p == q { 2.0 } else { 1.0 } / 12.0;
    lam(i, k) * dot(grad[j], grad[l])
        - lam(i, l) * dot(grad[j], grad[k])
        - lam(j, k) * dot(grad[i], grad[l])
        + lam(j, l) * dot(grad[i], grad[k])
}

fn dot(a: (f64, f64), b: (f64, f64)) -> f64 {
    a.0 * b.0 + a.1 * b.1
}

fn cross(a: (f64, f64), b: (f64, f64)) -> f64 {
    a.0 * b.1 - a.1 * b.0
}

/// Dense solve `A·X = B` with `B` of shape `n × m` (row-major), partial
/// pivoting. Returns `X` in the same shape.
fn solve_dense(a: &[f64], n: usize, b: &[f64], m: usize) -> Result<Vec<f64>, VectorFemError> {
    let lu = LuFactor::new(a.to_vec(), n).ok_or(VectorFemError::Singular)?;
    let mut out = vec![0.0; n * m];
    for c in 0..m {
        let col: Vec<f64> = (0..n).map(|r| b[r * m + c]).collect();
        let sol = lu.solve(&col);
        for r in 0..n {
            out[r * m + c] = sol[r];
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    #[test]
    fn rectangular_pec_guide_te10_cutoff() {
        // Uniform ε = 1 inside a PEC box of width a = 2 (b = 1). The lowest
        // mode is TE10 with β² = k₀² − (π/a)².
        let domain = RectDomain {
            width: 2.0,
            height: 1.0,
            nx: 16,
            ny: 8,
        };
        let k0: f64 = 3.0;
        let wl = Wavelength::from_m(2.0 * PI / k0);
        let exact_beta_sq = k0 * k0 - (PI / 2.0) * (PI / 2.0);
        let exact_neff = sqrt(exact_beta_sq) / k0;
        let index = |_x: f64, _y: f64| 1.0;
        let neff = vector_mode_neff(&domain, &index, wl, exact_neff, 60)
            .unwrap()
            .n_eff;
        let rel = fabs(neff - exact_neff) / exact_neff;
        assert!(rel < 0.01, "FEM {neff} vs exact {exact_neff} (rel {rel})");
    }

    #[test]
    fn dielectric_slab_te0_matches_exact() {
        // Wide-slab limit: the TE0 field is uniform across the width, so the
        // PEC side walls do not perturb it. The exact index comes from the
        // slab solver. Core y ∈ (3.5, 4.5) in a box of height 8.
        use crate::slab::{slab_modes, SlabLayer, SlabPolarization};
        let (n1, n2) = (1.5, 1.45);
        let wl = Wavelength::from_m(1.0);
        let slab = SlabLayer {
            n_core: n1,
            n_clad: n2,
            thickness: 1.0,
        };
        let exact = slab_modes(&slab, wl, SlabPolarization::Te)[0].n_eff;
        let domain = RectDomain {
            width: 2.0,
            height: 8.0,
            nx: 4,
            ny: 32,
        };
        let index = move |_x: f64, y: f64| if y > 3.5 && y < 4.5 { n1 } else { n2 };
        let neff = vector_mode_neff(&domain, &index, wl, exact, 60)
            .unwrap()
            .n_eff;
        let rel = fabs(neff - exact) / exact;
        assert!(rel < 1e-3, "FEM {neff} vs slab TE0 {exact} (rel {rel})");
    }

    #[test]
    fn strip_agrees_with_scalar_solver_in_weak_guidance() {
        // Weak guidance (Δ ≈ 0.7 %): the vector correction to the scalar
        // index is small, so the two solvers should agree closely.
        use crate::fem::{fem_fundamental_mode, IndexFn};
        let (n1, n2) = (1.5, 1.49);
        let wl = Wavelength::from_m(1.0);
        let index = move |x: f64, y: f64| {
            if (x - 1.0).abs() < 0.5 && (y - 4.0).abs() < 0.5 {
                n1
            } else {
                n2
            }
        };
        let index_ref: &IndexFn = &index;
        let domain = RectDomain {
            width: 2.0,
            height: 8.0,
            nx: 8,
            ny: 16,
        };
        let (beta, _) = fem_fundamental_mode(&domain, index_ref, wl, 80);
        let scalar = beta / wl.wavenumber();
        let vector = vector_mode_neff(&domain, index_ref, wl, scalar, 60)
            .unwrap()
            .n_eff;
        assert!(
            fabs(vector - scalar) < 2e-3,
            "vector {vector} vs scalar {scalar}"
        );
    }

    fn uniform_box() -> RectDomain {
        RectDomain {
            width: 2.0,
            height: 1.0,
            nx: 8,
            ny: 4,
        }
    }

    #[test]
    fn rejects_bad_inputs() {
        let one = |_x: f64, _y: f64| 1.0;
        let wl = Wavelength::from_m(2.0 * PI / 3.0);
        let mut d = uniform_box();
        d.nx = 0;
        assert!(matches!(
            vector_mode_neff(&d, &one, wl, 0.5, 20),
            Err(VectorFemError::InvalidParameter(_))
        ));
        let d = uniform_box();
        assert!(matches!(
            vector_mode_neff(&d, &one, wl, 0.5, 0),
            Err(VectorFemError::InvalidParameter(_))
        ));
        assert!(matches!(
            vector_mode_neff(&d, &one, wl, f64::NAN, 20),
            Err(VectorFemError::InvalidParameter(_))
        ));
        // A non-positive index at a centroid used to pass straight into the matrices.
        let bad = |_x: f64, _y: f64| 0.0;
        assert!(matches!(
            vector_mode_neff(&d, &bad, wl, 0.5, 20),
            Err(VectorFemError::InvalidParameter(_))
        ));
    }

    #[test]
    fn refuses_a_mesh_too_large_for_dense_storage() {
        // 100 × 100 cells give about 30 000 edges. The check must run before
        // any dense matrix is allocated.
        let big = RectDomain {
            width: 1.0,
            height: 1.0,
            nx: 100,
            ny: 100,
        };
        let one = |_x: f64, _y: f64| 1.0;
        let wl = Wavelength::from_m(2.0 * PI / 3.0);
        assert!(matches!(
            vector_mode_neff(&big, &one, wl, 0.5, 5),
            Err(VectorFemError::TooLarge { .. })
        ));
    }

    #[test]
    fn reports_residual_and_convergence() {
        let one = |_x: f64, _y: f64| 1.0;
        let k0: f64 = 3.0;
        let wl = Wavelength::from_m(2.0 * PI / k0);
        let d = RectDomain {
            width: 2.0,
            height: 1.0,
            nx: 16,
            ny: 8,
        };
        let mode = vector_mode_neff(&d, &one, wl, 0.9, 60).unwrap();
        assert!(mode.converged, "{mode:?}");
        assert!(mode.residual < 1e-6, "{mode:?}");
    }

    #[test]
    fn evanescent_solution_is_an_error_not_a_real_index() {
        // Below the TE10 cutoff (k0 = 1 < π/2) the mode decays, so β² < 0.
        // The old code returned √|β²|/k0 here, which looked like a guided index.
        let one = |_x: f64, _y: f64| 1.0;
        let wl = Wavelength::from_m(2.0 * PI / 1.0);
        let d = RectDomain {
            width: 2.0,
            height: 1.0,
            nx: 16,
            ny: 8,
        };
        let r = vector_mode_neff(&d, &one, wl, 0.5, 60);
        assert!(matches!(r, Err(VectorFemError::NotGuided { .. })), "{r:?}");
    }
}

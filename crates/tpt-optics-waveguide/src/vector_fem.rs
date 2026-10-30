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

use crate::fem::{IndexFn, RectDomain};

/// Effective index of the mode closest to `(k₀·n_guess)²` in the full-vector
/// hybrid edge–nodal formulation on `domain` with PEC walls.
///
/// `iterations` inverse-iteration steps are used; 30–60 is usually enough
/// when `n_guess` is within a few percent of the answer.
pub fn vector_mode_neff(
    domain: &RectDomain,
    index: &IndexFn<'_>,
    wl: Wavelength,
    n_guess: f64,
    iterations: usize,
) -> f64 {
    let k0 = wl.wavenumber();
    let k0sq = k0 * k0;
    let nodes = domain.node_coords();
    let n_nodes = nodes.len();
    let nx = domain.nx;
    let ny = domain.ny;
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
    let w = solve_dense(&k_mat, nn, &gt, ne);

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
    let lu = LuFactor::new(shifted, ne);
    let mut x = vec![1.0; ne];
    for _ in 0..iterations {
        let y = mat_vec(&b_mat, ne, &x);
        let z = lu.solve(&y);
        let norm = sqrt(z.iter().map(|v| v * v).sum::<f64>());
        if norm == 0.0 {
            break;
        }
        x = z.iter().map(|v| v / norm).collect();
    }
    let ax = mat_vec(&a_mat, ne, &x);
    let bx = mat_vec(&b_mat, ne, &x);
    let num: f64 = x.iter().zip(&ax).map(|(p, q)| p * q).sum();
    let den: f64 = x.iter().zip(&bx).map(|(p, q)| p * q).sum();
    let beta_sq = num / den;
    sqrt(fabs(beta_sq)) / k0
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

fn mat_vec(a: &[f64], n: usize, x: &[f64]) -> Vec<f64> {
    (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
        .collect()
}

/// Dense solve `A·X = B` with `B` of shape `n × m` (row-major), partial
/// pivoting. Returns `X` in the same shape.
fn solve_dense(a: &[f64], n: usize, b: &[f64], m: usize) -> Vec<f64> {
    let lu = LuFactor::new(a.to_vec(), n);
    let mut out = vec![0.0; n * m];
    for c in 0..m {
        let col: Vec<f64> = (0..n).map(|r| b[r * m + c]).collect();
        let sol = lu.solve(&col);
        for r in 0..n {
            out[r * m + c] = sol[r];
        }
    }
    out
}

/// LU factorisation with partial pivoting, stored for repeated solves.
struct LuFactor {
    lu: Vec<f64>,
    piv: Vec<usize>,
    n: usize,
}

impl LuFactor {
    fn new(mut a: Vec<f64>, n: usize) -> Self {
        let mut piv: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let mut best = k;
            for i in (k + 1)..n {
                if fabs(a[i * n + k]) > fabs(a[best * n + k]) {
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
            if pivot == 0.0 {
                continue;
            }
            for i in (k + 1)..n {
                let f = a[i * n + k] / pivot;
                a[i * n + k] = f;
                for j in (k + 1)..n {
                    a[i * n + j] -= f * a[k * n + j];
                }
            }
        }
        Self { lu: a, piv, n }
    }

    fn solve(&self, b: &[f64]) -> Vec<f64> {
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
            let d = self.lu[i * n + i];
            x[i] = if d == 0.0 { 0.0 } else { x[i] / d };
        }
        x
    }
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
        let neff = vector_mode_neff(&domain, &index, wl, exact_neff, 60);
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
        let neff = vector_mode_neff(&domain, &index, wl, exact, 60);
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
        let vector = vector_mode_neff(&domain, index_ref, wl, scalar, 60);
        assert!(
            fabs(vector - scalar) < 2e-3,
            "vector {vector} vs scalar {scalar}"
        );
    }
}

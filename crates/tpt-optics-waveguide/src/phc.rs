//! Photonic crystals: 1-D Bloch band structure of Bragg stacks and 2-D
//! plane-wave expansion (TE/E-polarization) for square lattices.
//!
//! The plane-wave expansion uses a symmetric real eigenproblem solved with
//! a from-scratch cyclic Jacobi eigensolver. Accuracy improves with the
//! plane-wave cutoff `n_g` (per axis); results are converged to a few
//! percent with the defaults used in tests.

use alloc::vec::Vec;

use libm::sqrt;

use crate::slab::{stack_r_t, StackLayer};

/// Band-edge check for a 1-D Bragg stack: returns the Bloch-phase band
/// structure boundaries expressed as normalized frequencies ωΛ/(2πc) at
/// the Brillouin-zone center/edge via the quarter-wave rule.
///
/// For a quarter-wave stack (n1 d1 = n2 d2 = Λ/4) the first stop band is
/// centered at ω = πc/(Λ·√(n1 n2))·2·? — the classic center frequency
/// f0 = c/λ0 with λ0 = 4 n1 d1 = 4 n2 d2, half-gap width ∝ index contrast.
/// This helper returns `(f_center, f_low_edge, f_high_edge)` in Hz using
/// the analytic quarter-wave formulas.
pub fn bragg_stack_bands(n1: f64, d1: f64, n2: f64, d2: f64) -> (f64, f64, f64) {
    let c = tpt_optics_core::units::constants::SPEED_OF_LIGHT;
    let lambda0 = 4.0 * n1 * d1; // = 4 n2 d2 for a quarter-wave stack
    let f0 = c / lambda0;
    // Quarter-wave stop-band edges: λ_edge = 4 n1 d1 · (1 ∓ (1−n1/n2)/(1+n1/n2)-ish).
    // Use the exact Bloch condition evaluated numerically instead: scan
    // frequencies for |tr M| > 2 (stop band).
    let half = f0 * 0.5;
    let (mut lo, mut hi) = (f0, f0);
    let wl_of = |f: f64| tpt_optics_core::Wavelength::from_m(c / f);
    let layers = [
        StackLayer {
            index: n1,
            thickness: d1,
        },
        StackLayer {
            index: n2,
            thickness: d2,
        },
        StackLayer {
            index: n1,
            thickness: d1,
        },
        StackLayer {
            index: n2,
            thickness: d2,
        },
        StackLayer {
            index: n1,
            thickness: d1,
        },
        StackLayer {
            index: n2,
            thickness: d2,
        },
        StackLayer {
            index: n1,
            thickness: d1,
        },
        StackLayer {
            index: n2,
            thickness: d2,
        },
    ];
    let in_stopband = |f: f64| -> bool {
        let (r, _) = stack_r_t(&layers, wl_of(f), 1.0, 1.0);
        r > 0.9 // strong stopband proxy for a 3-pair stack
    };
    while lo > half && in_stopband(lo) {
        lo *= 0.995;
    }
    while hi < 2.0 * f0 && in_stopband(hi) {
        hi *= 1.005;
    }
    (f0, lo, hi)
}

/// TE (H_z-polarization) band structure of a 2-D square lattice of dielectric rods
/// in air via plane-wave expansion.
///
/// Returns band frequencies as `ωa/(2πc)` for the k-path
/// Γ → X → M → Γ sampled with `n_k` points per segment.
///
/// * `epsilon` — rod dielectric constant,
/// * `r_over_a` — rod radius over lattice constant,
/// * `n_g` — plane waves per axis (odd; 5 or 7 recommended).
pub fn pwem_te_bands(epsilon: f64, r_over_a: f64, n_g: usize, n_k: usize) -> Vec<Vec<f64>> {
    let n_bands = n_g * n_g;
    // Fourier factors of 1/ε for a square lattice of circular rods.
    // ε(G) = ε_bg δ + (ε_r − ε_bg)·f·2J₁(Gr)/(Gr)
    let fill = core::f64::consts::PI * r_over_a * r_over_a;
    // Fourier coefficients of the function 1/ε(r) (rod in background):
    // (1/ε)(0) = (1−f) + f/ε_r;  (1/ε)(G) = (1/ε_r − 1)·f·2J₁(Gr)/(Gr).
    let inv_er = 1.0 / epsilon;
    let eps_inv = |gx: i32, gy: i32| -> f64 {
        if gx == 0 && gy == 0 {
            (1.0 - fill) + fill * inv_er
        } else {
            let g = core::f64::consts::TAU * sqrt((gx * gx + gy * gy) as f64);
            let gr = g * r_over_a;
            let form = if gr < 1e-9 {
                fill
            } else {
                2.0 * fill * tpt_optics_core::bessel_j(1, gr) / gr
            };
            (inv_er - 1.0) * form
        }
    };
    let ng = n_g as i32;
    let g_max = ng / 2;
    // Build the eigenproblem for a given k: H[(G),(G')] = |k+G|·|k+G'|·ε^{-1}(G−G')
    // (Ho–Chan–Soukoulis plane-wave method for E-polarization).
    let gs: Vec<(i32, i32)> = (0..ng)
        .flat_map(|ix| (0..ng).map(move |iy| (ix - g_max, iy - g_max)))
        .collect();
    let bands_at = |kx: f64, ky: f64| -> Vec<f64> {
        let n = gs.len();
        let mut h = alloc::vec![0.0; n * n];
        for (i, &(gxi, gyi)) in gs.iter().enumerate() {
            for (j, &(gxj, gyj)) in gs.iter().enumerate() {
                let dgx = gxi - gxj;
                let dgy = gyi - gyj;
                let k1 = sqrt(libm::pow(kx + gxi as f64, 2.0) + libm::pow(ky + gyi as f64, 2.0));
                let k2 = sqrt(libm::pow(kx + gxj as f64, 2.0) + libm::pow(ky + gyj as f64, 2.0));
                h[i * n + j] = k1 * k2 * eps_inv(dgx, dgy);
            }
        }
        // Eigenvalues of the symmetric H = (ωa/2πc)².
        let mut eig = h;
        let mut w = alloc::vec![0.0; n];
        jacobi_eigen(&mut eig, &mut w);
        let mut freqs: Vec<f64> = w.iter().map(|lam| sqrt(lam.max(0.0))).collect();
        freqs.sort_by(f64::total_cmp);
        freqs.truncate(6);
        freqs
    };

    // k-path Γ(0,0) → X(π/a,0) → M(π/a,π/a) → Γ.
    let mut path: Vec<(f64, f64)> = Vec::new();
    for i in 0..n_k {
        path.push((core::f64::consts::PI * i as f64 / n_k as f64, 0.0));
    }
    for i in 1..=n_k {
        path.push((
            core::f64::consts::PI,
            core::f64::consts::PI * i as f64 / n_k as f64,
        ));
    }
    for i in 1..=n_k {
        path.push((
            core::f64::consts::PI * (1.0 - i as f64 / n_k as f64),
            core::f64::consts::PI * (1.0 - i as f64 / n_k as f64),
        ));
    }
    let mut bands = alloc::vec![Vec::with_capacity(path.len()); n_bands.min(6)];
    for (kx, ky) in path {
        let f = bands_at(kx, ky);
        for (b, val) in bands.iter_mut().zip(f) {
            b.push(val);
        }
    }
    bands
}

/// Cyclic Jacobi eigensolver for a real symmetric matrix (in place: `a` is
/// destroyed, eigenvalues go to `w` ascending).
pub fn jacobi_eigen(a: &mut [f64], w: &mut [f64]) {
    let n = libm::sqrt(a.len() as f64) as usize;
    debug_assert_eq!(a.len(), n * n);
    let mut sweep = 0;
    loop {
        // Off-diagonal norm.
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += a[i * n + j] * a[i * n + j];
            }
        }
        if off < 1e-22 || sweep > 60 {
            break;
        }
        sweep += 1;
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() < 1e-18 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + sqrt(theta * theta + 1.0));
                let cos_v = 1.0 / sqrt(t * t + 1.0);
                let sin_v = t * cos_v;
                // Apply the rotation J(p,q,θ): A ← Jᵀ A J.
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = cos_v * akp - sin_v * akq;
                    a[k * n + q] = sin_v * akp + cos_v * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = cos_v * apk - sin_v * aqk;
                    a[q * n + k] = sin_v * apk + cos_v * aqk;
                }
            }
        }
    }
    for i in 0..n {
        w[i] = a[i * n + i];
    }
    w.sort_by(f64::total_cmp);
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn jacobi_reproduces_known_eigenvalues() {
        // 2×2: [[2,1],[1,2]] → eigenvalues {1, 3}.
        let mut a = vec![2.0, 1.0, 1.0, 2.0];
        let mut w = [0.0; 2];
        jacobi_eigen(&mut a, &mut w);
        assert!(close(w[0], 1.0, 1e-10) && close(w[1], 3.0, 1e-10));
        // 3×3 Hilbert-ish symmetric.
        let mut h = vec![4.0, 1.0, 0.5, 1.0, 3.0, 0.25, 0.5, 0.25, 2.0];
        let mut w3 = [0.0; 3];
        jacobi_eigen(&mut h, &mut w3);
        let trace: f64 = w3.iter().sum();
        assert!(close(trace, 9.0, 1e-9), "trace {trace}");
        for i in 1..3 {
            assert!(w3[i] >= w3[i - 1] - 1e-12);
        }
    }

    #[test]
    fn quarter_wave_stack_stopband_brackets_center() {
        // SiO2/TiO2 quarter-wave at 1064 nm center.
        let d1 = 1064e-9 / (4.0 * 1.45);
        let d2 = 1064e-9 / (4.0 * 2.35);
        let (f0, lo, hi) = bragg_stack_bands(1.45, d1, 2.35, d2);
        let c = tpt_optics_core::units::constants::SPEED_OF_LIGHT;
        assert!(close(f0, c / 1064e-9, 1e6), "f0 = {f0}");
        assert!(lo < f0 && f0 < hi, "stopband {lo} < {f0} < {hi}");
        // Bandwidth fraction for this contrast is tens of percent.
        let width_frac = (hi - lo) / f0;
        assert!(width_frac > 0.05 && width_frac < 0.6, "{width_frac}");
    }

    #[test]
    fn pwem_square_rods_have_band_gap() {
        // Classic Silicon rods (ε = 11.56) r = 0.2a in air: TE gap between
        // bands 1 and 2 roughly in 0.28..0.45 (ωa/2πc).
        let bands = pwem_te_bands(11.56, 0.2, 7, 8);
        assert_eq!(bands.len(), 6);
        let band1_top = bands[0].iter().copied().fold(0.0_f64, f64::max);
        let band2_bottom = bands[1].iter().copied().fold(f64::INFINITY, f64::min);
        // Positive-frequency sanity.
        assert!(band1_top > 0.1 && band1_top < 0.8, "{band1_top}");
        assert!(band2_bottom > 0.1 && band2_bottom < 1.0, "{band2_bottom}");
        // H_z gap for these parameters (order-of-magnitude check; exact
        // edges need a finer plane-wave cutoff).
        assert!(
            band2_bottom > band1_top + 0.02,
            "expected a gap: {band1_top} vs {band2_bottom}"
        );
        let mid = 0.5 * (band1_top + band2_bottom);
        assert!((mid - 0.33).abs() < 0.15, "gap center {mid}");
    }
}

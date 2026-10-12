//! Bessel functions of the first kind `Jν` and modified second kind `Kν`.
//!
//! From-scratch implementations (series + asymptotic expansion) used by the
//! fiber mode solver and coherence theory. Accuracy is ≳1e-8 across the
//! argument ranges the optics kernels exercise; these are *not* general-purpose
//! special-function libraries.

use libm::sqrt;

/// Bessel function of the first kind, order 0.
pub fn bessel_j0(x: f64) -> f64 {
    bessel_j(0, x)
}

/// Bessel function of the first kind, order 1.
pub fn bessel_j1(x: f64) -> f64 {
    bessel_j(1, x)
}

/// Bessel function of the first kind `Jν(x)` for integer `ν ≥ 0`.
///
/// Series summation for small arguments (|x| < 16), asymptotic expansion
/// above. Uses upward recurrence from `J₀, J₁` (stable for `ν ≤ x`).
///
/// # Panics
/// Panics in debug builds for negative orders.
pub fn bessel_j(nu: i32, x: f64) -> f64 {
    debug_assert!(nu >= 0, "negative order not supported");
    let ax = x.abs();
    if ax < 1e-10 {
        return if nu == 0 {
            1.0
        } else {
            // Jν(x) ≈ (x/2)^ν / ν!
            let mut inv_fact = 1.0;
            for k in 2..=nu {
                inv_fact /= k as f64;
            }
            libm::pow(0.5 * x, nu as f64) * inv_fact
        };
    }
    // Large-x asymptotics are used for J₀ and J₁ only. Higher orders go
    // through Miller's recurrence, which is accurate at every x.
    if ax > 20.0 && nu < 2 {
        return j_asymptotic(nu, x);
    }
    // J_ν(−x) = (−1)^ν J_ν(x): compute at ax > 0 and restore the sign.
    let sign = if x < 0.0 && nu % 2 == 1 { -1.0 } else { 1.0 };
    let (j0, j1) = j0_j1_series(ax);
    let pos = match nu {
        0 => j0,
        1 => j1,
        // Upward recurrence J_{ν+1} = 2ν/x · J_ν − J_{ν−1} is unstable once
        // ν exceeds x, so use Miller's backward recurrence for ν ≥ 2.
        _ => miller_j(nu as usize, ax),
    };
    sign * pos
}

/// `J_n(x)` for `x > 0` by Miller's backward recurrence, normalized with
/// `J₀ + 2Σ J_{2k} = 1`. Stable for every order; the start index is chosen
/// well beyond `max(n, x)` so the seed error has decayed by the target.
fn miller_j(n: usize, x: f64) -> f64 {
    let base = n.max(libm::ceil(x) as usize);
    let m = 2 * ((base + 20 + libm::sqrt(40.0 * base as f64) as usize) / 2);
    let mut jp = 0.0; // j_{k+1}
    let mut jk = 1e-30; // j_k, starting at k = m
                        // Σ over the normalization identity, with j_m already included.
    let mut total = 2.0 * jk;
    let mut j_n = 0.0;
    for k in (1..=m).rev() {
        let jkm1 = 2.0 * k as f64 / x * jk - jp;
        jp = jk;
        jk = jkm1;
        let idx = k - 1; // jk now holds j_idx
        if idx == n {
            j_n = jk;
        }
        if idx == 0 {
            total += jk;
        } else if idx % 2 == 0 {
            total += 2.0 * jk;
        }
        if jk.abs() > 1e100 {
            // Rescale everything together; the ratio j_n / total is unchanged.
            jk *= 1e-100;
            jp *= 1e-100;
            total *= 1e-100;
            j_n *= 1e-100;
        }
    }
    j_n / total
}

/// J₀ and J₁ by series: `J₀ = Σ (−1)^k (x²/4)^k/(k!)²`,
/// `J₁ = (x/2)·Σ (−1)^k (x²/4)^k/(k!(k+1)!)`.
fn j0_j1_series(x: f64) -> (f64, f64) {
    let half_x2 = 0.25 * x * x;
    let mut s0 = 1.0;
    let mut term0 = 1.0;
    let mut s1 = 1.0;
    let mut term1 = 1.0;
    for k in 1..200u32 {
        term0 *= -half_x2 / (k as f64 * k as f64);
        s0 += term0;
        term1 *= -half_x2 / (k as f64 * (k as f64 + 1.0));
        s1 += term1;
        if term0.abs() < 1e-18 && term1.abs() < 1e-18 {
            break;
        }
    }
    (s0, 0.5 * x * s1)
}

/// Asymptotic form for large |x| (DLMF 10.17.3):
/// `Jν(x) ≈ √(2/(πx))·(P·cos ω − Q·sin ω)`, `ω = x − νπ/2 − π/4`, with the
/// standard P/Q polynomial corrections through 4th/3rd order.
fn j_asymptotic(nu: i32, x: f64) -> f64 {
    let ax = x.abs();
    let mu = 4.0 * nu as f64 * nu as f64;
    let chi = ax - nu as f64 * core::f64::consts::FRAC_PI_2 - core::f64::consts::FRAC_PI_4;
    let a = mu - 1.0;
    let e8 = 8.0 * ax;
    let p = 1.0 - a * (mu - 9.0) / (2.0 * e8 * e8)
        + a * (mu - 9.0) * (mu - 17.0) * (mu - 25.0) / (24.0 * e8 * e8 * e8 * e8);
    let q = a / e8 - a * (mu - 9.0) * (mu - 25.0) / (6.0 * e8 * e8 * e8);
    let amp = sqrt(2.0 / (core::f64::consts::PI * ax));
    let val = amp * (libm::cos(chi) * p - libm::sin(chi) * q);
    if x < 0.0 {
        // Jν(−x) = (−1)^ν Jν(x)
        if nu % 2 == 0 {
            val
        } else {
            -val
        }
    } else {
        val
    }
}

/// Modified Bessel function of the second kind, order 0.
pub fn bessel_k0(x: f64) -> f64 {
    bessel_k(0, x)
}

/// Modified Bessel function of the second kind, order 1.
pub fn bessel_k1(x: f64) -> f64 {
    bessel_k(1, x)
}

/// Modified Bessel function of the second kind `Kν(x)` for `x > 0`, integer
/// `ν ≥ 0`.
///
/// Builds `K₀, K₁` from `I₀, I₁` (Abramowitz & Stegun 9.8.x polynomials for
/// `x ≤ 2`, asymptotic above), then upward recurrence.
///
/// # Panics
/// Panics for `x ≤ 0` or negative orders.
pub fn bessel_k(nu: i32, x: f64) -> f64 {
    debug_assert!(nu >= 0, "negative order not supported");
    debug_assert!(x > 0.0, "Kν defined for x > 0");
    let (k0, k1) = k0_k1(x);
    match nu {
        0 => k0,
        1 => k1,
        _ => {
            let mut km = k0;
            let mut k = k1;
            for i in 1..nu {
                let next = 2.0 * i as f64 / x * k + km;
                km = core::mem::replace(&mut k, next);
            }
            k
        }
    }
}

/// Asymptotic expansion of `K_ν(x)` for large `x` (A&S 9.7.1), with
/// `mu = 4ν²`: `√(π/2x)·e^{−x}·Σ_k a_k/x^k`, where
/// `a_k = Π_{j≤k}(mu − (2j−1)²) / (k!·8^k)`. Summed until the terms stop
/// shrinking (optimal truncation; error ≈ e^{−2x}, about 1e−7 at x = 8).
fn k_asymptotic(mu: f64, x: f64) -> f64 {
    let mut sum = 1.0;
    let mut term = 1.0;
    let mut smallest = f64::INFINITY;
    for k in 1..60u32 {
        let odd = (2 * k - 1) as f64;
        term *= (mu - odd * odd) / (k as f64 * 8.0 * x);
        if term.abs() > smallest {
            break;
        }
        smallest = term.abs();
        sum += term;
        if term.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    sqrt(core::f64::consts::FRAC_PI_2 / x) * libm::exp(-x) * sum
}

fn k0_k1(x: f64) -> (f64, f64) {
    if x <= 8.0 {
        // Frobenius series with digamma terms (exact for integer order):
        //   K0 = −ln(x/2)·I0 + Σ (x²/4)^k/(k!)² · ψ(k+1)
        //   K1 = 1/x + Σ (x/2)^{2k+1}/(k!(k+1)!) · [ln(x/2) − ½ψ(k+1) − ½ψ(k+2)]
        let gamma = 0.577_215_664_901_532_9;
        // ψ(n) for integer n ≥ 1: −γ + H_{n−1}
        let psi = |n: u32| -> f64 {
            let mut h = 0.0;
            for j in 1..n {
                h += 1.0 / j as f64;
            }
            h - gamma
        };
        let ln_half_x = libm::log(x * 0.5);
        let y = 0.25 * x * x;

        let mut i0 = 1.0; // I0 partial: Σ y^k/(k!)²
        let mut k0 = psi(1); // Σ y^k/(k!)² · ψ(k+1)
        let mut term = 1.0;
        for k in 1..40u32 {
            term *= y / (k as f64 * k as f64);
            i0 += term;
            k0 += term * psi(k + 1);
            if term.abs() < 1e-19 {
                break;
            }
        }
        let k0 = -ln_half_x * i0 + k0;

        let half_x = 0.5 * x;
        let mut k1 = 1.0 / x;
        let mut term = half_x; // (x/2)^1/(0!·1!)
        for k in 0..40u32 {
            let bracket = ln_half_x - 0.5 * psi(k + 1) - 0.5 * psi(k + 2);
            k1 += term * bracket;
            let next = k + 1;
            term *= y / (next as f64 * (next as f64 + 1.0));
            if term.abs() < 1e-19 {
                break;
            }
        }
        (k0, k1)
    } else {
        (k_asymptotic(0.0, x), k_asymptotic(4.0, x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn j_reference_values() {
        // A&S tables / scipy reference values.
        assert!(close(bessel_j0(0.0), 1.0, 1e-12));
        assert!(close(bessel_j0(1.0), 0.765_197_686_557_966_6, 1e-12));
        assert!(close(bessel_j1(1.0), 0.440_050_585_744_933_5, 1e-12));
        assert!(close(bessel_j0(5.0), -0.177_596_771_314_338_3, 1e-10));
        assert!(close(bessel_j1(5.0), -0.327_579_137_591_465_2, 1e-10));
        assert!(close(bessel_j0(10.0), -0.245_935_764_451_348_4, 1e-9));
        assert!(close(bessel_j1(10.0), 0.043_472_746_168_861_4, 1e-9));
        // Series branch edge.
        assert!(close(bessel_j0(20.0), 0.167_024_664_262_906_2, 1e-8));
        // Asymptotic branch (documented ~1e-6 accuracy for large x).
        assert!(close(bessel_j1(30.0), -0.118_750_56, 1e-5));
        // Recurrence branch: J3(5)
        assert!(close(bessel_j(3, 5.0), 0.364_831_230_613_667_2, 1e-10));
        // Zeros: J0 first zero ≈ 2.4048255577
        let a = bessel_j0(2.404_5);
        let b = bessel_j0(2.405_1);
        assert!(
            a < 0.0 && b > 0.0 || a > 0.0 && b < 0.0,
            "sign change expected"
        );
    }

    #[test]
    fn k_reference_values() {
        assert!(close(bessel_k0(1.0), 0.421_024_438_240_708_3, 1e-10));
        assert!(close(bessel_k1(1.0), 0.601_907_230_197_234_6, 1e-10));
        assert!(close(bessel_k0(0.5), 0.924_419_071_227_666_3, 1e-9));
        assert!(close(bessel_k1(0.5), 1.656_441_120_003_300_5, 1e-9));
        // Series branch edge.
        assert!(close(bessel_k0(5.0), 0.003_691_098_334_042_568, 1e-12));
        assert!(close(bessel_k1(5.0), 0.004_044_613_445_478_336, 1e-12));
        // Recurrence branch.
        assert!(close(bessel_k(2, 2.0), 0.253_759_754_566_055_9, 1e-10));
    }

    #[test]
    fn k_diverges_at_zero() {
        assert!(bessel_k0(1e-6) > 10.0);
    }
}

#[cfg(test)]
mod modified_k_reference_tests {
    use super::*;

    #[test]
    fn k0_k1_asymptotic_branch_matches_tabulated_values() {
        // Reference: K0(10) = 1.7780062316e-5, K1(10) = 1.8648773454e-5.
        let k0 = bessel_k0(10.0);
        let k1 = bessel_k1(10.0);
        assert!(
            (k0 / 1.778_006_231_6e-5 - 1.0).abs() < 1e-6,
            "K0(10) = {k0}"
        );
        assert!(
            (k1 / 1.864_877_345_4e-5 - 1.0).abs() < 1e-6,
            "K1(10) = {k1}"
        );
    }

    #[test]
    fn series_and_asymptotic_branches_join_smoothly() {
        // x = 8 uses the series, just above 8 the asymptotic expansion.
        let (series, asym) = (bessel_k0(8.0), bessel_k0(8.0 + 1e-9));
        assert!((series / asym - 1.0).abs() < 1e-6, "{series} vs {asym}");
        let (series1, asym1) = (bessel_k1(8.0), bessel_k1(8.0 + 1e-9));
        assert!((series1 / asym1 - 1.0).abs() < 1e-6, "{series1} vs {asym1}");
    }
}

#[cfg(test)]
mod high_order_tests {
    use super::*;

    #[test]
    fn high_order_matches_published_values() {
        // Upward recurrence loses everything for ν > x; Miller's stays exact.
        let j10 = bessel_j(10, 1.0);
        assert!(
            (j10 / 2.630_615_123_687_4e-10 - 1.0).abs() < 1e-6,
            "J10(1) = {j10}"
        );
        let j2 = bessel_j(2, 3.0);
        assert!((j2 - 0.486_091_260_585_891).abs() < 1e-9, "J2(3) = {j2}");
    }

    #[test]
    fn negative_argument_follows_parity() {
        // J_ν(−x) = (−1)^ν J_ν(x).
        // J3(3) from the recurrence J3 = (4/3)J2 − J1 with J1(3) = 0.339058958.
        assert!((bessel_j(3, -3.0) + 0.309_062_722_25).abs() < 1e-9);
        assert!((bessel_j(2, -3.0) - 0.486_091_260_585_891).abs() < 1e-9);
    }
}

#[cfg(test)]
mod large_x_high_order_tests {
    use super::*;

    #[test]
    fn high_order_at_large_x_matches_stable_recurrence() {
        // Upward recurrence is stable for ν ≤ x; build J₅(30) from the
        // asymptotic J₀, J₁ and compare with the Miller-based value.
        let x = 30.0;
        let (mut jm, mut j) = (bessel_j(0, x), bessel_j(1, x));
        for k in 1..5 {
            let next = 2.0 * k as f64 / x * j - jm;
            jm = j;
            j = next;
        }
        let miller = bessel_j(5, x);
        assert!((miller - j).abs() < 1e-6, "J5(30): {miller} vs {j}");
    }
}

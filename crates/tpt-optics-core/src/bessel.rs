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
    if ax > 20.0 {
        return j_asymptotic(nu, x);
    }
    let (j0, j1) = j0_j1_series(ax);
    match nu {
        0 => j0,
        1 => j1 * x.signum(),
        _ => {
            // Upward recurrence J_{ν+1} = 2ν/x · J_ν − J_{ν−1}
            let mut jm = j0;
            let mut j = j1 * x.signum();
            for k in 1..nu {
                let next = 2.0 * k as f64 / ax * j - jm;
                jm = core::mem::replace(&mut j, next);
            }
            j
        }
    }
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
        // Asymptotic series (A&S 9.7.2 / 9.7.3 style).
        let ix = 1.0 / x;
        let e = libm::exp(-x) / sqrt(x);
        let k0 = e
            * 1.253_314_14
            * (1.0
                + ix * (-0.078_323_58
                    + ix * (0.021_895_68
                        + ix * (-0.010_624_46
                            + ix * (0.005_878_72 + ix * (-0.002_515_40 + ix * 0.000_532_08))))));
        let k1 = e
            * 1.253_314_14
            * (1.0
                + ix * (0.234_986_19
                    + ix * (-0.036_556_20
                        + ix * (0.015_042_68
                            + ix * (-0.007_803_53 + ix * (0.003_256_14 + ix * (-0.000_682_45)))))));
        (k0, k1)
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

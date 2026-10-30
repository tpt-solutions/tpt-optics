//! Monte-Carlo tolerance analysis with compensators.
//!
//! Perturb a design's parameters with random (normal or uniform) errors,
//! optionally re-optimize a *compensator* parameter (e.g., the final
//! air-space) to absorb the error, and report the resulting distribution of
//! a merit function.

use alloc::vec::Vec;

use tpt_optics_core::Rng;

/// One tolerance: parameter index, distribution type, and magnitude.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Tolerance {
    /// Normal error with the given standard deviation.
    Normal {
        /// Parameter index.
        index: usize,
        /// Standard deviation.
        sigma: f64,
    },
    /// Uniform error on [−half_width, +half_width].
    Uniform {
        /// Parameter index.
        index: usize,
        /// Half width.
        half_width: f64,
    },
}

impl Tolerance {
    /// Draw one error sample.
    fn sample(&self, rng: &mut Rng) -> (usize, f64) {
        match *self {
            Tolerance::Normal { index, sigma } => (index, rng.normal(0.0, sigma)),
            Tolerance::Uniform { index, half_width } => (index, rng.range(-half_width, half_width)),
        }
    }
}

/// Statistics of the merit distribution.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ToleranceStats {
    /// Mean merit.
    pub mean: f64,
    /// Merit standard deviation.
    pub sigma: f64,
    /// 90th percentile of the merit.
    pub p90: f64,
    /// Fraction of runs exceeding the pass limit.
    pub failure_fraction: f64,
}

/// Run the tolerance Monte Carlo.
///
/// * `base` — nominal parameter vector,
/// * `tolerances` — perturbations applied per run,
/// * `compensator` — optional `(index, target_residual)` compensator:
///   after perturbation, that parameter is solved so that the residual
///   function hits the target (1-D bisection on the merit),
/// * `merit` — scalar merit of a parameter vector (smaller is better),
/// * `passes` — number of Monte-Carlo runs,
/// * `pass_limit` — merit above which the unit "fails".
pub fn analyze<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    merit: &M,
    passes: usize,
    pass_limit: f64,
    seed: u64,
) -> ToleranceStats
where
    M: Fn(&[f64]) -> f64,
{
    let mut rng = Rng::new(seed);
    let values: Vec<f64> = (0..passes)
        .map(|_| one_run(base, tolerances, compensator, merit, &mut rng))
        .collect();
    summarize(values, pass_limit)
}

/// Like [`analyze`], but each run draws from its own generator, seeded from
/// `seed` and the run index. The result does not depend on how the runs are
/// scheduled, so [`analyze_parallel`] gives exactly the same statistics.
pub fn analyze_seeded<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    merit: &M,
    passes: usize,
    pass_limit: f64,
    seed: u64,
) -> ToleranceStats
where
    M: Fn(&[f64]) -> f64,
{
    let values: Vec<f64> = (0..passes)
        .map(|i| {
            let mut rng = Rng::new(run_seed(seed, i as u64));
            one_run(base, tolerances, compensator, merit, &mut rng)
        })
        .collect();
    summarize(values, pass_limit)
}

/// Parallel form of [`analyze_seeded`], on the rayon thread pool. The
/// statistics are identical to the serial result for the same seed.
#[cfg(feature = "parallel")]
pub fn analyze_parallel<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    merit: &M,
    passes: usize,
    pass_limit: f64,
    seed: u64,
) -> ToleranceStats
where
    M: Fn(&[f64]) -> f64 + Sync,
{
    use rayon::prelude::*;
    let values: Vec<f64> = (0..passes)
        .into_par_iter()
        .map(|i| {
            let mut rng = Rng::new(run_seed(seed, i as u64));
            one_run(base, tolerances, compensator, merit, &mut rng)
        })
        .collect();
    summarize(values, pass_limit)
}

/// Seed of run `index` derived from the base `seed` (SplitMix64 mixing).
fn run_seed(seed: u64, index: u64) -> u64 {
    fn mix(mut z: u64) -> u64 {
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    mix(seed ^ mix(index.wrapping_add(0x9E37_79B9_7F4A_7C15)))
}

/// One tolerance run: perturb the parameters, apply the compensator, and
/// return the merit.
fn one_run<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    merit: &M,
    rng: &mut Rng,
) -> f64
where
    M: Fn(&[f64]) -> f64,
{
    let mut x = base.to_vec();
    for t in tolerances {
        let (i, e) = t.sample(rng);
        x[i] += e;
    }
    if let Some((ci, target)) = compensator {
        x[ci] = solve_compensator(merit, &x, ci, target);
    }
    merit(&x)
}

/// Statistics of a list of merit values.
fn summarize(mut values: Vec<f64>, pass_limit: f64) -> ToleranceStats {
    values.sort_by(f64::total_cmp);
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
    let p90 = values[libm::round(0.9 * (values.len() as f64 - 1.0)) as usize];
    // A NaN merit (undefined or blocked trace) is a failure, not a pass.
    let failures = values
        .iter()
        .filter(|v| v.is_nan() || **v > pass_limit)
        .count();
    ToleranceStats {
        mean,
        sigma: libm::sqrt(var),
        p90,
        failure_fraction: failures as f64 / n,
    }
}

/// Solve one parameter so the merit hits `target` (bisection over
/// `x[ci] ± search_half_width`, merit assumed monotone in that interval).
fn solve_compensator<M>(merit: &M, x: &[f64], ci: usize, target: f64) -> f64
where
    M: Fn(&[f64]) -> f64,
{
    let half_width = 0.1 * libm::fabs(x[ci]).max(1e-3) * 10.0;
    let g = |c: f64| -> f64 {
        let mut xt = x.to_vec();
        xt[ci] = c;
        merit(&xt) - target
    };
    let (mut lo, mut hi) = (x[ci] - half_width, x[ci] + half_width);
    let (glo, ghi) = (g(lo), g(hi));
    if glo * ghi > 0.0 {
        return x[ci]; // no bracket — keep nominal
    }
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if g(lo) * g(mid) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn tolerancing_spreads_the_merit() {
        // Merit = x0 (a thickness-like parameter); tolerance sigma = 0.1.
        let base = [1.0];
        let tols = [Tolerance::Normal {
            index: 0,
            sigma: 0.1,
        }];
        let stats = analyze(&base, &tols, None, &|x| x[0], 20_000, 1.5, 99);
        assert!(
            close(stats.mean, 1.0, 0.01),
            "mean {mean}",
            mean = stats.mean
        );
        assert!(close(stats.sigma, 0.1, 0.01), "sigma {}", stats.sigma);
        assert!(
            close(stats.failure_fraction, 1.0 - 0.99997, 0.01),
            "fail {} (P(x>1.5) under N(1,0.1) ≈ 2.9e-6)",
            stats.failure_fraction
        );
    }

    #[test]
    fn compensator_absorbs_the_error() {
        // Two mirrors distance-like parameters; merit depends on the sum
        // x0 + x1. Tolerate x0 but compensate with x1: merit stays put.
        let base = [10.0, 20.0];
        let tols = [Tolerance::Normal {
            index: 0,
            sigma: 0.2,
        }];
        let stats = analyze(
            &base,
            &tols,
            Some((1, 30.0)), // solve x1 so merit = 30
            &|x| x[0] + x[1],
            5_000,
            31.0,
            7,
        );
        assert!(close(stats.mean, 30.0, 1e-6), "mean {}", stats.mean);
        assert!(stats.sigma < 1e-6, "sigma {}", stats.sigma);
        assert!(close(stats.failure_fraction, 0.0, 1e-9));
    }

    #[test]
    fn uniform_tolerances() {
        let base = [0.0];
        let tols = [Tolerance::Uniform {
            index: 0,
            half_width: 0.1,
        }];
        let stats = analyze(&base, &tols, None, &|x| x[0], 20_000, 0.05, 3);
        // Uniform(±0.1): P(x > 0.05) = 0.25.
        assert!(
            close(stats.failure_fraction, 0.25, 0.02),
            "fail {}",
            stats.failure_fraction
        );
        assert!(stats.mean.abs() < 0.005);
    }
}

#[cfg(test)]
mod nan_tolerance_tests {
    use super::*;

    #[test]
    fn nan_merit_counts_as_failure() {
        // A merit that is undefined everywhere must fail every sample, not pass.
        let stats = analyze(&[0.0], &[], None, &|_x: &[f64]| f64::NAN, 8, 1.0, 1);
        assert_eq!(stats.failure_fraction, 1.0);
    }
}

#[cfg(test)]
mod seeded_tests {
    use super::*;

    fn merit(x: &[f64]) -> f64 {
        libm::sqrt(x.iter().map(|v| v * v).sum::<f64>())
    }

    fn tols() -> [Tolerance; 2] {
        [
            Tolerance::Normal {
                index: 0,
                sigma: 0.1,
            },
            Tolerance::Uniform {
                index: 1,
                half_width: 0.2,
            },
        ]
    }

    #[test]
    fn seeded_analysis_is_reproducible() {
        let a = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 200, 0.5, 7);
        let b = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 200, 0.5, 7);
        assert_eq!(a, b);
        let c = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 200, 0.5, 8);
        assert_ne!(a.mean, c.mean);
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_analysis_matches_the_serial_result_exactly() {
        let serial = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 500, 0.5, 11);
        let parallel = analyze_parallel(&[1.0, 2.0], &tols(), None, &merit, 500, 0.5, 11);
        assert_eq!(serial, parallel);
    }
}

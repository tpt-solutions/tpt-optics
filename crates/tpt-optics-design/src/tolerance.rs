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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    /// Normal error with standard deviation `sigma`, limited to `[−limit, +limit]`.
    /// This is how a vendor's tolerance is usually stated: a sigma inside a
    /// hard bound. Draws outside the bound are redrawn, and a draw that stays
    /// outside after a fixed number of tries is clamped.
    TruncatedNormal {
        /// Parameter index.
        index: usize,
        /// Standard deviation of the untruncated normal.
        sigma: f64,
        /// Hard bound on the error (must be positive).
        limit: f64,
    },
}

/// Redraws allowed before a truncated normal is clamped to its bound.
const TRUNCATION_TRIES: usize = 1_000;

impl Tolerance {
    /// Parameter index this tolerance perturbs.
    fn index(&self) -> usize {
        match *self {
            Tolerance::Normal { index, .. }
            | Tolerance::Uniform { index, .. }
            | Tolerance::TruncatedNormal { index, .. } => index,
        }
    }

    /// Whether the tolerance's magnitudes are usable.
    fn is_valid(&self) -> bool {
        match *self {
            Tolerance::TruncatedNormal { sigma, limit, .. } => {
                sigma.is_finite() && limit.is_finite() && limit > 0.0
            }
            _ => true,
        }
    }

    /// Draw one error sample.
    fn sample(&self, rng: &mut Rng) -> (usize, f64) {
        match *self {
            Tolerance::Normal { index, sigma } => (index, rng.normal(0.0, sigma)),
            Tolerance::Uniform { index, half_width } => (index, rng.range(-half_width, half_width)),
            Tolerance::TruncatedNormal {
                index,
                sigma,
                limit,
            } => {
                let mut e = rng.normal(0.0, sigma);
                for _ in 0..TRUNCATION_TRIES {
                    if libm::fabs(e) <= limit {
                        break;
                    }
                    e = rng.normal(0.0, sigma);
                }
                (index, e.clamp(-limit, limit))
            }
        }
    }
}

/// Statistics of the merit distribution.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

/// Why a tolerance analysis was refused before any run.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ToleranceError {
    /// The analysis needs at least one Monte-Carlo run.
    NoRuns,
    /// A tolerance or the compensator names a parameter that the base
    /// vector does not have.
    IndexOutOfRange {
        /// The index that was named.
        index: usize,
        /// Length of the base parameter vector.
        len: usize,
    },
    /// A tolerance has a magnitude that cannot be sampled, such as a
    /// truncation limit that is not positive.
    InvalidTolerance {
        /// Parameter index of the tolerance.
        index: usize,
    },
}

impl core::fmt::Display for ToleranceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoRuns => write!(f, "tolerance analysis needs at least one run"),
            Self::IndexOutOfRange { index, len } => write!(
                f,
                "parameter index {index} is out of range for {len} parameters"
            ),
            Self::InvalidTolerance { index } => write!(
                f,
                "tolerance on parameter {index} has a magnitude that cannot be sampled"
            ),
        }
    }
}

impl core::error::Error for ToleranceError {}

/// Check the inputs shared by every analysis entry point.
fn validate(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    passes: usize,
) -> Result<(), ToleranceError> {
    if passes == 0 {
        return Err(ToleranceError::NoRuns);
    }
    if let Some(bad) = tolerances.iter().find(|t| !t.is_valid()) {
        return Err(ToleranceError::InvalidTolerance { index: bad.index() });
    }
    let indices = tolerances
        .iter()
        .map(Tolerance::index)
        .chain(compensator.map(|(ci, _)| ci));
    for index in indices {
        if index >= base.len() {
            return Err(ToleranceError::IndexOutOfRange {
                index,
                len: base.len(),
            });
        }
    }
    Ok(())
}

/// Run the tolerance Monte Carlo.
///
/// * `base` — nominal parameter vector,
/// * `tolerances` — perturbations applied per run,
/// * `compensator` — optional `(index, target_residual)` compensator:
///   after perturbation, that parameter is solved so that the residual
///   function hits the target (1-D bisection on the merit),
/// * `merit` — scalar merit of a parameter vector (smaller is better),
/// * `passes` — number of Monte-Carlo runs (at least one),
/// * `pass_limit` — merit above which the unit "fails".
///
/// # Errors
///
/// [`ToleranceError::NoRuns`] when `passes` is zero, and
/// [`ToleranceError::IndexOutOfRange`] when a tolerance or the compensator
/// names a parameter outside `base`.
pub fn analyze<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    merit: &M,
    passes: usize,
    pass_limit: f64,
    seed: u64,
) -> Result<ToleranceStats, ToleranceError>
where
    M: Fn(&[f64]) -> f64,
{
    validate(base, tolerances, compensator, passes)?;
    let mut rng = Rng::new(seed);
    let values: Vec<f64> = (0..passes)
        .map(|_| one_run(base, tolerances, compensator, merit, &mut rng))
        .collect();
    Ok(summarize(values, pass_limit))
}

/// Like [`analyze`], but each run draws from its own generator, seeded from
/// `seed` and the run index. The result does not depend on how the runs are
/// scheduled, so [`analyze_parallel`] gives exactly the same statistics.
///
/// # Errors
///
/// The same as [`analyze`].
pub fn analyze_seeded<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    merit: &M,
    passes: usize,
    pass_limit: f64,
    seed: u64,
) -> Result<ToleranceStats, ToleranceError>
where
    M: Fn(&[f64]) -> f64,
{
    validate(base, tolerances, compensator, passes)?;
    let values: Vec<f64> = (0..passes)
        .map(|i| {
            let mut rng = Rng::new(run_seed(seed, i as u64));
            one_run(base, tolerances, compensator, merit, &mut rng)
        })
        .collect();
    Ok(summarize(values, pass_limit))
}

/// Parallel form of [`analyze_seeded`], on the rayon thread pool. The
/// statistics are identical to the serial result for the same seed.
///
/// # Errors
///
/// The same as [`analyze`].
#[cfg(feature = "parallel")]
pub fn analyze_parallel<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    compensator: Option<(usize, f64)>,
    merit: &M,
    passes: usize,
    pass_limit: f64,
    seed: u64,
) -> Result<ToleranceStats, ToleranceError>
where
    M: Fn(&[f64]) -> f64 + Sync,
{
    use rayon::prelude::*;
    validate(base, tolerances, compensator, passes)?;
    let values: Vec<f64> = (0..passes)
        .into_par_iter()
        .map(|i| {
            let mut rng = Rng::new(run_seed(seed, i as u64));
            one_run(base, tolerances, compensator, merit, &mut rng)
        })
        .collect();
    Ok(summarize(values, pass_limit))
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
/// return the merit. The caller has validated every index.
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

/// Statistics of a non-empty list of merit values.
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
        let stats = analyze(&base, &tols, None, &|x| x[0], 20_000, 1.5, 99).unwrap();
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
        )
        .unwrap();
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
        let stats = analyze(&base, &tols, None, &|x| x[0], 20_000, 0.05, 3).unwrap();
        // Uniform(±0.1): P(x > 0.05) = 0.25.
        assert!(
            close(stats.failure_fraction, 0.25, 0.02),
            "fail {}",
            stats.failure_fraction
        );
        assert!(stats.mean.abs() < 0.005);
    }

    #[test]
    fn zero_runs_is_an_error_not_a_panic() {
        let result = analyze(&[0.0], &[], None, &|x: &[f64]| x[0], 0, 1.0, 1);
        assert_eq!(result, Err(ToleranceError::NoRuns));
    }

    #[test]
    fn out_of_range_tolerance_is_an_error_not_a_panic() {
        let tols = [Tolerance::Normal {
            index: 3,
            sigma: 0.1,
        }];
        let result = analyze(&[0.0, 1.0], &tols, None, &|x: &[f64]| x[0], 4, 1.0, 1);
        assert_eq!(
            result,
            Err(ToleranceError::IndexOutOfRange { index: 3, len: 2 })
        );
    }

    #[test]
    fn out_of_range_compensator_is_an_error_not_a_panic() {
        let result = analyze(&[0.0], &[], Some((1, 0.0)), &|x: &[f64]| x[0], 4, 1.0, 1);
        assert_eq!(
            result,
            Err(ToleranceError::IndexOutOfRange { index: 1, len: 1 })
        );
    }
}

#[cfg(test)]
mod nan_tolerance_tests {
    use super::*;

    #[test]
    fn nan_merit_counts_as_failure() {
        // A merit that is undefined everywhere must fail every sample, not pass.
        let stats = analyze(&[0.0], &[], None, &|_x: &[f64]| f64::NAN, 8, 1.0, 1).unwrap();
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
        let a = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 200, 0.5, 7).unwrap();
        let b = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 200, 0.5, 7).unwrap();
        assert_eq!(a, b);
        let c = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 200, 0.5, 8).unwrap();
        assert_ne!(a.mean, c.mean);
    }

    #[test]
    fn seeded_analysis_refuses_zero_runs() {
        let result = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 0, 0.5, 7);
        assert_eq!(result, Err(ToleranceError::NoRuns));
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_analysis_matches_the_serial_result_exactly() {
        let serial = analyze_seeded(&[1.0, 2.0], &tols(), None, &merit, 500, 0.5, 11).unwrap();
        let parallel = analyze_parallel(&[1.0, 2.0], &tols(), None, &merit, 500, 0.5, 11).unwrap();
        assert_eq!(serial, parallel);
    }
}

#[cfg(test)]
mod truncated_tests {
    use super::*;

    #[test]
    fn truncated_normal_stays_within_its_limit() {
        // Merit |x| with sigma 0.1 would exceed 0.2 now and then. Truncated at
        // 0.15, it never does.
        let tols = [Tolerance::TruncatedNormal {
            index: 0,
            sigma: 0.1,
            limit: 0.15,
        }];
        let stats =
            analyze_seeded(&[0.0], &tols, None, &|x: &[f64]| x[0].abs(), 5000, 0.2, 3).unwrap();
        assert_eq!(stats.failure_fraction, 0.0);
        assert!(stats.p90 <= 0.15 + 1e-12, "{}", stats.p90);
    }

    #[test]
    fn non_positive_truncation_limit_is_refused() {
        // A zero limit cannot be sampled, so it is refused up front rather
        // than looping.
        let tols = [Tolerance::TruncatedNormal {
            index: 0,
            sigma: 0.1,
            limit: 0.0,
        }];
        let result = analyze(&[0.0], &tols, None, &|x: &[f64]| x[0], 4, 1.0, 1);
        assert_eq!(result, Err(ToleranceError::InvalidTolerance { index: 0 }));
    }
}

/// Sensitivity table: the merit spread caused by each tolerance on its own.
///
/// Each row is the statistics of a run with only that tolerance applied, so a
/// row says how much that error alone moves the merit. For a linear merit
/// with independent errors, the row sigmas add in quadrature to the spread of
/// the full analysis.
///
/// # Errors
///
/// The same as [`analyze_seeded`].
pub fn sensitivity_table<M>(
    base: &[f64],
    tolerances: &[Tolerance],
    merit: &M,
    passes: usize,
    pass_limit: f64,
    seed: u64,
) -> Result<Vec<(usize, ToleranceStats)>, ToleranceError>
where
    M: Fn(&[f64]) -> f64,
{
    tolerances
        .iter()
        .map(|t| {
            let stats = analyze_seeded(
                base,
                core::slice::from_ref(t),
                None,
                merit,
                passes,
                pass_limit,
                seed,
            )?;
            Ok((t.index(), stats))
        })
        .collect()
}

#[cfg(test)]
mod sensitivity_tests {
    use super::*;

    #[test]
    fn rows_add_in_quadrature_for_a_linear_merit() {
        // Merit x0 + x1 with independent errors of sigma 0.10 and 0.05. Each row
        // reports its own sigma, and the full spread is their quadrature sum.
        let base = [0.0, 0.0];
        let tols = [
            Tolerance::Normal {
                index: 0,
                sigma: 0.10,
            },
            Tolerance::Normal {
                index: 1,
                sigma: 0.05,
            },
        ];
        let merit = |x: &[f64]| x[0] + x[1];
        let rows = sensitivity_table(&base, &tols, &merit, 20_000, 1.0, 9).unwrap();
        assert_eq!(rows[0].0, 0);
        assert_eq!(rows[1].0, 1);
        assert!(
            (rows[0].1.sigma - 0.10).abs() < 0.005,
            "{}",
            rows[0].1.sigma
        );
        assert!(
            (rows[1].1.sigma - 0.05).abs() < 0.003,
            "{}",
            rows[1].1.sigma
        );
        let full = analyze_seeded(&base, &tols, None, &merit, 20_000, 1.0, 9).unwrap();
        let quadrature = (0.10_f64 * 0.10 + 0.05 * 0.05).sqrt();
        assert!((full.sigma - quadrature).abs() < 0.005, "{}", full.sigma);
    }

    #[test]
    fn invalid_tolerance_is_refused_before_any_row() {
        let tols = [Tolerance::TruncatedNormal {
            index: 0,
            sigma: 0.1,
            limit: -1.0,
        }];
        let result = sensitivity_table(&[0.0], &tols, &|x: &[f64]| x[0], 4, 1.0, 1);
        assert!(matches!(
            result,
            Err(ToleranceError::InvalidTolerance { index: 0 })
        ));
    }
}

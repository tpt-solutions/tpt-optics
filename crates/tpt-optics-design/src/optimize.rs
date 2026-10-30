//! Generic optimizers: damped least squares (Levenberg–Marquardt),
//! simulated annealing, and a real-coded genetic algorithm.
//!
//! All operate on an `f64` parameter vector with optional box bounds and a
//! scalar objective; no globals, deterministic given the seed.

use alloc::vec;
use alloc::vec::Vec;

use tpt_optics_core::Rng;

/// Per-parameter bound pair.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ParamBound {
    /// Lower limit.
    pub lo: f64,
    /// Upper limit.
    pub hi: f64,
}

impl ParamBound {
    /// Clamp a value into the bound.
    pub fn clamp(&self, v: f64) -> f64 {
        v.clamp(self.lo, self.hi)
    }
}

/// Why an optimization run stopped.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum StopReason {
    /// The last accepted step lowered the cost by less than a relative
    /// `1e-12`, or the cost reached exactly zero.
    Converged,
    /// The iteration or generation budget ran out.
    MaxIterations,
    /// No trial step lowered the cost, even at the largest damping.
    Stalled,
    /// The observer asked the run to stop.
    Cancelled,
}

/// What an observer tells the optimiser to do after an iteration.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Control {
    /// Keep going.
    Continue,
    /// Stop now and return the best point so far.
    Cancel,
}

/// Snapshot passed to an observer at the start of each iteration.
#[derive(Copy, Clone, Debug)]
pub struct Progress<'a> {
    /// Zero-based iteration number.
    pub iteration: usize,
    /// Current sum of squared residuals.
    pub cost: f64,
    /// Current parameter vector.
    pub x: &'a [f64],
}

/// Result of an optimization run.
#[derive(Clone, Debug)]
pub struct OptimizeResult {
    /// Best parameter vector found.
    pub x: Vec<f64>,
    /// Objective at `x`.
    pub f: f64,
    /// Objective evaluations spent.
    pub evaluations: usize,
    /// Best objective after each iteration, step, or generation.
    pub history: Vec<f64>,
    /// Why the run stopped.
    pub stop: StopReason,
}

/// Objective trait: `f(x) -> cost` (minimization).
pub trait Objective {
    /// Evaluate the objective.
    fn eval(&self, x: &[f64]) -> f64;
}

impl<F> Objective for F
where
    F: Fn(&[f64]) -> f64,
{
    fn eval(&self, x: &[f64]) -> f64 {
        self(x)
    }
}

/// Central-difference Jacobian of a vector-valued residual `r(x)` with
/// step `h` per component (scaled by `|x_i|` when large).
pub fn numeric_jacobian<R: Fn(&[f64], &mut Vec<f64>)>(
    r: &R,
    x: &[f64],
    n_res: usize,
    h: f64,
) -> Vec<f64> {
    let n = x.len();
    let mut jac = alloc::vec![0.0; n_res * n];
    let mut rp = alloc::vec![0.0; n_res];
    let mut rm = alloc::vec![0.0; n_res];
    let mut xs = x.to_vec();
    for j in 0..n {
        let hj = h.max(1e-9 * libm::fabs(x[j]));
        xs[j] = x[j] + hj;
        r(&xs, &mut rp);
        xs[j] = x[j] - hj;
        r(&xs, &mut rm);
        xs[j] = x[j];
        for i in 0..n_res {
            jac[i * n + j] = (rp[i] - rm[i]) / (2.0 * hj);
        }
    }
    jac
}

/// Solve `(JᵀJ + λ·diag(JᵀJ))·δ = −Jᵀr` via Gaussian elimination.
fn damped_normal_equations(jac: &[f64], r: &[f64], n: usize, lambda: f64) -> Option<Vec<f64>> {
    // Normal equations A = JᵀJ + λ diag, b = −Jᵀ r.
    let mut a = alloc::vec![0.0; n * n];
    let mut b = alloc::vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0;
            for k in 0..r.len() {
                acc += jac[k * n + i] * jac[k * n + j];
            }
            a[i * n + j] = acc + if i == j { lambda * acc.max(1e-30) } else { 0.0 };
        }
    }
    for i in 0..n {
        let mut acc = 0.0;
        for k in 0..r.len() {
            acc += jac[k * n + i] * r[k];
        }
        b[i] = -acc;
    }
    // Gaussian elimination with partial pivoting.
    for col in 0..n {
        let (mut piv, mut prow) = (libm::fabs(a[col * n + col]), col);
        for row in (col + 1)..n {
            let v = libm::fabs(a[row * n + col]);
            if v > piv {
                piv = v;
                prow = row;
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
        let pv = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pv;
            a[row * n + col] = 0.0;
            for k in (col + 1)..n {
                a[row * n + k] -= factor * a[col * n + k];
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = alloc::vec![0.0; n];
    for row in (0..n).rev() {
        let mut sum = b[row];
        for k in (row + 1)..n {
            sum -= a[row * n + k] * x[k];
        }
        x[row] = sum / a[row * n + row];
    }
    Some(x)
}

/// Levenberg–Marquardt minimization of a vector residual's squared norm.
///
/// * `r(x, out)` — residual evaluator writing `out.len()` components,
/// * `x0` — start, `bounds` — box constraints (clamped steps),
/// * `iterations` — max LM iterations.
pub fn damped_least_squares<R>(
    r: &R,
    x0: &[f64],
    n_res: usize,
    bounds: Option<&[ParamBound]>,
    iterations: usize,
) -> OptimizeResult
where
    R: Fn(&[f64], &mut Vec<f64>),
{
    damped_least_squares_observed(r, x0, n_res, bounds, iterations, &mut |_| Control::Continue)
}

/// Levenberg–Marquardt with an observer. The observer sees the state at the
/// start of each iteration and can cancel the run by returning
/// [`Control::Cancel`]. The run then stops with [`StopReason::Cancelled`].
pub fn damped_least_squares_observed<R, C>(
    r: &R,
    x0: &[f64],
    n_res: usize,
    bounds: Option<&[ParamBound]>,
    iterations: usize,
    observer: &mut C,
) -> OptimizeResult
where
    R: Fn(&[f64], &mut Vec<f64>),
    C: FnMut(&Progress<'_>) -> Control,
{
    let jacobian =
        |x: &[f64]| numeric_jacobian(&|xv: &[f64], out: &mut Vec<f64>| r(xv, out), x, n_res, 1e-7);
    lm_observed_with(r, &jacobian, x0, n_res, bounds, iterations, observer)
}

/// Levenberg–Marquardt with a caller-supplied Jacobian.
///
/// `jacobian(x)` returns the `n_res × x.len()` matrix in row-major order. Use
/// it when the derivatives are known (for example from dual-number tracing)
/// and finite differences would be slow or inexact.
pub fn damped_least_squares_with_jacobian<R, J, C>(
    r: &R,
    jacobian: &J,
    x0: &[f64],
    n_res: usize,
    bounds: Option<&[ParamBound]>,
    iterations: usize,
    observer: &mut C,
) -> OptimizeResult
where
    R: Fn(&[f64], &mut Vec<f64>),
    J: Fn(&[f64]) -> Vec<f64>,
    C: FnMut(&Progress<'_>) -> Control,
{
    lm_observed_with(r, jacobian, x0, n_res, bounds, iterations, observer)
}

fn lm_observed_with<R, J, C>(
    r: &R,
    jacobian: &J,
    x0: &[f64],
    n_res: usize,
    bounds: Option<&[ParamBound]>,
    iterations: usize,
    observer: &mut C,
) -> OptimizeResult
where
    R: Fn(&[f64], &mut Vec<f64>),
    J: Fn(&[f64]) -> Vec<f64>,
    C: FnMut(&Progress<'_>) -> Control,
{
    let n = x0.len();
    let mut x = x0.to_vec();
    clamp_all(&mut x, bounds);
    let mut resid = alloc::vec![0.0; n_res];
    r(&x, &mut resid);
    let mut cost: f64 = resid.iter().map(|v| v * v).sum();
    let mut lambda = 1e-3_f64;
    let mut evals = 1usize;
    let mut history = Vec::with_capacity(iterations);
    let mut stop = StopReason::MaxIterations;
    for iteration in 0..iterations {
        if cost == 0.0 {
            stop = StopReason::Converged;
            break;
        }
        let progress = Progress {
            iteration,
            cost,
            x: &x,
        };
        if observer(&progress) == Control::Cancel {
            stop = StopReason::Cancelled;
            break;
        }
        let jac = jacobian(&x);
        let Some(step) = damped_normal_equations(&jac, &resid, n, lambda) else {
            lambda *= 10.0;
            history.push(cost);
            if lambda >= 1e8 {
                stop = StopReason::Stalled;
                break;
            }
            continue;
        };
        let mut trial = x.clone();
        for (t, s) in trial.iter_mut().zip(step.iter()) {
            *t += s;
        }
        clamp_all(&mut trial, bounds);
        let mut trial_resid = alloc::vec![0.0; n_res];
        r(&trial, &mut trial_resid);
        let trial_cost: f64 = trial_resid.iter().map(|v| v * v).sum();
        evals += 2;
        if trial_cost < cost {
            let relative = (cost - trial_cost) / cost;
            x = trial;
            resid = trial_resid;
            cost = trial_cost;
            lambda = (lambda * 0.3).max(1e-10);
            history.push(cost);
            if relative < 1e-12 {
                stop = StopReason::Converged;
                break;
            }
        } else {
            lambda = (lambda * 5.0).min(1e8);
            history.push(cost);
            if lambda >= 1e8 {
                stop = StopReason::Stalled;
                break;
            }
        }
    }
    OptimizeResult {
        x,
        f: cost,
        evaluations: evals,
        history,
        stop,
    }
}

fn clamp_all(x: &mut [f64], bounds: Option<&[ParamBound]>) {
    if let Some(bounds) = bounds {
        for (v, b) in x.iter_mut().zip(bounds.iter()) {
            *v = b.clamp(*v);
        }
    }
}

/// Simulated annealing with exponential cooling and per-parameter step
/// scales. Deterministic for a given seed.
pub fn simulated_annealing<O: Objective>(
    f: &O,
    x0: &[f64],
    bounds: Option<&[ParamBound]>,
    steps: usize,
    t0: f64,
    step_scale: &[f64],
    seed: u64,
) -> OptimizeResult {
    let mut rng = Rng::new(seed);
    let mut x = x0.to_vec();
    clamp_all(&mut x, bounds);
    let mut fx = f.eval(&x);
    let mut best = x.clone();
    let mut f_best = fx;
    let mut evals = 1usize;
    let mut history = Vec::with_capacity(steps);
    for k in 0..steps {
        let t = t0 * libm::exp(-(k as f64) / (steps as f64 / 6.0));
        let mut trial = x.clone();
        for (v, s) in trial.iter_mut().zip(step_scale.iter()) {
            *v += rng.normal(0.0, *s);
        }
        clamp_all(&mut trial, bounds);
        let ft = f.eval(&trial);
        evals += 1;
        let delta = ft - fx;
        if delta < 0.0 || (t > 1e-300 && rng.next_f64() < libm::exp(-delta / t)) {
            x = trial;
            fx = ft;
            if fx < f_best {
                f_best = fx;
                best.clone_from(&x);
            }
        }
        history.push(f_best);
    }
    OptimizeResult {
        x: best,
        f: f_best,
        evaluations: evals,
        history,
        stop: StopReason::MaxIterations,
    }
}

/// Real-coded genetic algorithm with tournament selection, blend crossover
/// (BLX-α) and Gaussian mutation. Minimizes `f`.
pub fn genetic_minimize<O: Objective>(
    f: &O,
    bounds: &[ParamBound],
    population: usize,
    generations: usize,
    seed: u64,
) -> OptimizeResult {
    let mut rng = Rng::new(seed);
    let n = bounds.len();
    let mut evals = 0usize;
    // Initial population (uniform in bounds).
    let mut pop: Vec<Vec<f64>> = Vec::with_capacity(population);
    for _ in 0..population {
        let ind: Vec<f64> = bounds.iter().map(|b| rng.range(b.lo, b.hi)).collect();
        pop.push(ind);
    }
    let mut scores: Vec<f64> = pop
        .iter()
        .map(|ind| {
            evals += 1;
            f.eval(ind)
        })
        .collect();
    let mut history = Vec::with_capacity(generations);
    for _ in 0..generations {
        // Elitism: carry the best individual unchanged.
        let mut best_i = 0;
        for (i, s) in scores.iter().enumerate() {
            if merit_key(*s) < merit_key(scores[best_i]) {
                best_i = i;
            }
        }
        let mut next = vec![pop[best_i].clone()];
        let mut next_scores = vec![scores[best_i]];
        while next.len() < population {
            // Tournament selection (size 3).
            let mut pick = || {
                let mut best = (rng.next_f64() * population as f64) as usize % population;
                for _ in 0..2 {
                    let cand = (rng.next_f64() * population as f64) as usize % population;
                    if merit_key(scores[cand]) < merit_key(scores[best]) {
                        best = cand;
                    }
                }
                best
            };
            let p1 = pick();
            let p2 = pick();
            // BLX-α crossover (α = 0.5).
            let alpha = 0.5;
            let mut child = alloc::vec![0.0; n];
            for j in 0..n {
                let lo = pop[p1][j].min(pop[p2][j]);
                let hi = pop[p1][j].max(pop[p2][j]);
                let span = hi - lo;
                child[j] = rng.range(lo - alpha * span, hi + alpha * span);
            }
            // Gaussian mutation (10 % chance per gene, 10 % of span).
            for j in 0..n {
                if rng.next_f64() < 0.1 {
                    let span = bounds[j].hi - bounds[j].lo;
                    child[j] += rng.normal(0.0, 0.1 * span);
                }
                child[j] = bounds[j].clamp(child[j]);
            }
            evals += 1;
            let s = f.eval(&child);
            next.push(child);
            next_scores.push(s);
        }
        pop = next;
        scores = next_scores;
        history.push(scores.iter().copied().fold(f64::INFINITY, f64::min));
    }
    let mut best_i = 0;
    for (i, s) in scores.iter().enumerate() {
        if merit_key(*s) < merit_key(scores[best_i]) {
            best_i = i;
        }
    }
    OptimizeResult {
        x: pop[best_i].clone(),
        f: scores[best_i],
        evaluations: evals,
        history,
        stop: StopReason::MaxIterations,
    }
}

/// Ordering key for merit values: NaN (an undefined or blocked trace) sorts
/// after every finite value, so it is never selected as best.
fn merit_key(f: f64) -> f64 {
    if f.is_nan() {
        f64::INFINITY
    } else {
        f
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn lm_solves_rosenbrock_ish() {
        // Least-squares form of the Rosenbrock banana:
        // r = [10(x2 − x1²), (1 − x1)] → minimum at (1, 1).
        let r = |x: &[f64], out: &mut Vec<f64>| {
            out[0] = 10.0 * (x[1] - x[0] * x[0]);
            out[1] = 1.0 - x[0];
        };
        let res = damped_least_squares(&r, &[-1.2, 1.0], 2, None, 200);
        assert!(close(res.x[0], 1.0, 1e-4), "x = {:?}", res.x);
        assert!(close(res.x[1], 1.0, 1e-3), "x = {:?}", res.x);
        assert!(res.f < 1e-8);
    }

    #[test]
    fn lm_with_bounds_respects_box() {
        let r = |x: &[f64], out: &mut Vec<f64>| {
            out[0] = x[0] - 5.0; // wants x0 = 5
            out[1] = x[1] * x[1] * 0.01; // wants x1 = 0
        };
        let bounds = [
            ParamBound { lo: -2.0, hi: 2.0 },
            ParamBound { lo: -1.0, hi: 1.0 },
        ];
        let res = damped_least_squares(&r, &[0.0, 0.5], 2, Some(&bounds), 60);
        assert!(res.x[0] <= 2.0 + 1e-12);
        assert!(close(res.x[0], 2.0, 1e-9), "x = {:?}", res.x);
    }

    #[test]
    fn annealing_minimizes_sphere() {
        // Sphere in a box: minimum at origin.
        let f = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>();
        let bounds = [
            ParamBound { lo: -2.0, hi: 2.0 },
            ParamBound { lo: -2.0, hi: 2.0 },
        ];
        let steps = [0.5, 0.5];
        let res = simulated_annealing(&f, &[1.5, -1.5], Some(&bounds), 20_000, 2.0, &steps, 42);
        assert!(res.f < 1e-3, "f = {}", res.f);
    }

    #[test]
    fn genetic_minimizes_rastrigin_ish() {
        // Multimodal: sin-ripple over a bowl — GA should beat the local
        // basin and get near the origin.
        let f = |x: &[f64]| {
            x.iter()
                .map(|v| v * v + 2.0 * (1.0 - libm::cos(3.0 * v)))
                .sum::<f64>()
        };
        let bounds = [
            ParamBound { lo: -3.0, hi: 3.0 },
            ParamBound { lo: -3.0, hi: 3.0 },
        ];
        let res = genetic_minimize(&f, &bounds, 60, 80, 7);
        assert!(res.f < 0.5, "f = {} at {:?}", res.f, res.x);
        // Determinism.
        let res2 = genetic_minimize(&f, &bounds, 60, 80, 7);
        assert!(close(res.f, res2.f, 1e-12));
    }
}

#[cfg(test)]
mod nan_tests {
    use super::*;

    #[test]
    fn nan_merits_are_never_reported_as_best() {
        // Points with x < 0 are "blocked" (NaN); the GA must still return the
        // best finite point, never the NaN.
        let f = |x: &[f64]| if x[0] < 0.0 { f64::NAN } else { x[0] * x[0] };
        let bounds = [ParamBound { lo: -1.0, hi: 1.0 }];
        let r = genetic_minimize(&f, &bounds, 20, 30, 3);
        assert!(r.f.is_finite(), "best f = {}", r.f);
        assert!(r.x[0] >= 0.0);
    }
}

#[cfg(test)]
mod observer_tests {
    use super::*;

    fn linear(x: &[f64], out: &mut Vec<f64>) {
        out.clear();
        out.push(x[0] - 1.0);
        out.push(x[1] + 2.0);
    }

    #[test]
    fn observer_can_cancel_the_run() {
        let mut seen = Vec::new();
        let res = damped_least_squares_observed(&linear, &[0.0, 0.0], 2, None, 50, &mut |p| {
            seen.push(p.iteration);
            if p.iteration >= 2 {
                Control::Cancel
            } else {
                Control::Continue
            }
        });
        assert_eq!(res.stop, StopReason::Cancelled);
        assert_eq!(seen, vec![0, 1, 2]);
    }

    #[test]
    fn observer_sees_the_cost_falling() {
        let mut costs = Vec::new();
        let _ = damped_least_squares_observed(&linear, &[0.0, 0.0], 2, None, 50, &mut |p| {
            costs.push(p.cost);
            Control::Continue
        });
        assert!(
            costs.len() >= 2 && costs.windows(2).all(|w| w[1] <= w[0]),
            "{costs:?}"
        );
    }
}

#[cfg(test)]
mod stop_reason_tests {
    use super::*;

    /// Linear residual with its minimum at `(1, −2)`.
    fn linear(x: &[f64], out: &mut Vec<f64>) {
        out.clear();
        out.push(x[0] - 1.0);
        out.push(x[1] + 2.0);
    }

    #[test]
    fn damped_least_squares_reports_convergence_and_history() {
        let r = linear;
        let res = damped_least_squares(&r, &[0.0, 0.0], 2, None, 50);
        assert_eq!(res.stop, StopReason::Converged);
        assert!(res.f < 1e-20, "{}", res.f);
        assert_ne!(res.history.len(), 0);
        assert!(
            res.history.windows(2).all(|w| w[1] <= w[0]),
            "{:?}",
            res.history
        );
    }

    #[test]
    fn running_out_of_iterations_is_reported() {
        let r = linear;
        let res = damped_least_squares(&r, &[0.0, 0.0], 2, None, 1);
        assert_eq!(res.stop, StopReason::MaxIterations);
    }

    #[test]
    fn a_flat_residual_is_reported_as_stalled() {
        // Zero Jacobian: no step can lower the cost, so damping grows until it
        // is refused.
        let r = |_x: &[f64], out: &mut Vec<f64>| {
            out.clear();
            out.push(1.0);
        };
        let res = damped_least_squares(&r, &[0.0], 1, None, 500);
        assert_eq!(res.stop, StopReason::Stalled);
    }
}

#[cfg(test)]
mod supplied_jacobian_tests {
    use super::*;

    fn residual(x: &[f64], out: &mut Vec<f64>) {
        out.clear();
        out.push(x[0] * x[0] - 2.0);
        out.push(x[0] + 3.0 * x[1] - 5.0);
    }

    fn analytic(x: &[f64]) -> Vec<f64> {
        // Row-major 2 × 2: d/dx0, d/dx1 of the two residuals.
        alloc::vec![2.0 * x[0], 0.0, 1.0, 3.0]
    }

    #[test]
    fn analytic_jacobian_gives_the_same_solution_as_finite_differences() {
        let numeric = damped_least_squares(&residual, &[1.0, 1.0], 2, None, 60);
        let exact = damped_least_squares_with_jacobian(
            &residual,
            &analytic,
            &[1.0, 1.0],
            2,
            None,
            60,
            &mut |_| Control::Continue,
        );
        assert!((numeric.f - exact.f).abs() < 1e-18 + 1e-12 * numeric.f.max(1.0));
        for (a, b) in numeric.x.iter().zip(&exact.x) {
            assert!((a - b).abs() < 1e-8, "{a} vs {b}");
        }
    }
}

//! Parameter sweeps: full-factorial grids and Latin-hypercube samples over a
//! box of parameter bounds.
//!
//! Both return parameter vectors that can be fed to any evaluator, such as a
//! merit function or a tolerance run. Latin-hypercube samples are seeded, so
//! the same seed gives the same sample on every platform.

use alloc::vec::Vec;

use tpt_optics_core::rng::Rng;

/// Why a sweep could not be built.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SweepError {
    /// No parameters, or the sample count is zero.
    Empty,
    /// A bound is non-finite or `lo > hi`, or a step count does not match the
    /// number of parameters.
    InvalidBounds,
    /// The Sobol generator supports up to [`SOBOL_MAX_DIM`] parameters.
    UnsupportedDimension,
}

impl core::fmt::Display for SweepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Empty => write!(f, "sweep needs at least one parameter and one point"),
            Self::InvalidBounds => write!(f, "sweep bounds must be finite with lo <= hi"),
            Self::UnsupportedDimension => write!(
                f,
                "Sobol sequences support at most {SOBOL_MAX_DIM} parameters"
            ),
        }
    }
}

impl core::error::Error for SweepError {}

fn check_bounds(bounds: &[(f64, f64)]) -> Result<(), SweepError> {
    if bounds.is_empty() {
        return Err(SweepError::Empty);
    }
    if bounds
        .iter()
        .any(|&(lo, hi)| !(lo.is_finite() && hi.is_finite() && lo <= hi))
    {
        return Err(SweepError::InvalidBounds);
    }
    Ok(())
}

/// Full-factorial grid: `steps[i]` evenly spaced values for parameter `i`,
/// spanning `bounds[i] = (lo, hi)` inclusively. The first parameter varies
/// fastest. A single step gives the midpoint of the interval.
///
/// # Errors
/// [`SweepError::Empty`] when there are no parameters or a step count is
/// zero, and [`SweepError::InvalidBounds`] for a bad bound or a mismatched
/// `steps` length.
pub fn grid(bounds: &[(f64, f64)], steps: &[usize]) -> Result<Vec<Vec<f64>>, SweepError> {
    check_bounds(bounds)?;
    if steps.len() != bounds.len() {
        return Err(SweepError::InvalidBounds);
    }
    if steps.contains(&0) {
        return Err(SweepError::Empty);
    }
    let axes: Vec<Vec<f64>> = bounds
        .iter()
        .zip(steps)
        .map(|(&(lo, hi), &n)| {
            if n == 1 {
                alloc::vec![0.5 * (lo + hi)]
            } else {
                (0..n)
                    .map(|k| lo + (hi - lo) * k as f64 / (n - 1) as f64)
                    .collect()
            }
        })
        .collect();
    let total: usize = steps.iter().product();
    let mut out = Vec::with_capacity(total);
    for index in 0..total {
        let mut rest = index;
        let point = axes
            .iter()
            .map(|axis| {
                let v = axis[rest % axis.len()];
                rest /= axis.len();
                v
            })
            .collect();
        out.push(point);
    }
    Ok(out)
}

/// Largest number of parameters the Sobol generator supports.
pub const SOBOL_MAX_DIM: usize = 10;

/// Joe–Kuo direction-number data for dimensions 2 to 10: the degree `s`, the
/// packed polynomial coefficients `a`, and the initial odd integers `m₁…m_s`.
const SOBOL_TABLE: [(usize, u64, &[u64]); 9] = [
    (1, 0, &[1]),
    (2, 1, &[1, 3]),
    (3, 1, &[1, 3, 1]),
    (3, 2, &[1, 1, 1]),
    (4, 1, &[1, 1, 3, 3]),
    (4, 4, &[1, 3, 5, 13]),
    (5, 2, &[1, 1, 5, 5, 17]),
    (5, 4, &[1, 1, 5, 5, 5]),
    (5, 7, &[1, 1, 7, 11, 19]),
];

/// The 32 direction numbers `v[0..32]` of dimension `dim` (0-based).
fn direction_numbers(dim: usize) -> [u32; 32] {
    let mut v = [0u32; 32];
    if dim == 0 {
        for (i, slot) in v.iter_mut().enumerate() {
            *slot = 1u32 << (31 - i);
        }
        return v;
    }
    let (s, a, init) = SOBOL_TABLE[dim - 1];
    // m[i] is m_{i+1}; the first s values come from the table.
    let mut m = [0u64; 32];
    for (i, &value) in init.iter().enumerate() {
        m[i] = value;
    }
    for i in s..32 {
        // m_i = 2^s·m_{i−s} ⊕ m_{i−s} ⊕ Σ_k 2^k·a_k·m_{i−k}, with a_k the k-th
        // bit of the packed polynomial, and 1-based indices in the formula.
        let mut value = (m[i - s] << s) ^ m[i - s];
        for k in 1..s {
            if (a >> (s - 1 - k)) & 1 == 1 {
                value ^= m[i - k] << k;
            }
        }
        m[i] = value;
    }
    for (i, slot) in v.iter_mut().enumerate() {
        *slot = (m[i] << (31 - i)) as u32;
    }
    v
}

/// Sobol sequence of `n` points in the box `bounds`, in Gray-code order.
///
/// The first `2^k` points form a stratified set: every elementary interval
/// of volume `2^-k` in a pair of coordinates holds exactly one point. The
/// sequence is deterministic, so no seed is needed.
///
/// # Errors
/// [`SweepError::Empty`] for `n == 0` or no parameters,
/// [`SweepError::InvalidBounds`] for a bad bound, and
/// [`SweepError::UnsupportedDimension`] above [`SOBOL_MAX_DIM`] parameters.
pub fn sobol(bounds: &[(f64, f64)], n: usize) -> Result<Vec<Vec<f64>>, SweepError> {
    check_bounds(bounds)?;
    if n == 0 {
        return Err(SweepError::Empty);
    }
    let dims = bounds.len();
    if dims > SOBOL_MAX_DIM {
        return Err(SweepError::UnsupportedDimension);
    }
    let directions: Vec<[u32; 32]> = (0..dims).map(direction_numbers).collect();
    let mut x = alloc::vec![0u32; dims];
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let point = bounds
            .iter()
            .zip(&x)
            .map(|(&(lo, hi), &xi)| lo + (hi - lo) * (xi as f64 / 4_294_967_296.0))
            .collect();
        out.push(point);
        // Next point: XOR the direction number of the rightmost zero bit of k.
        let c = (!(k as u64)).trailing_zeros() as usize;
        for (xi, v) in x.iter_mut().zip(&directions) {
            *xi ^= v[c.min(31)];
        }
    }
    Ok(out)
}

/// Latin-hypercube sample of `n` points in the box `bounds`, seeded.
///
/// Each parameter's range is cut into `n` equal strata, and every stratum is
/// used exactly once. The position inside a stratum is random.
///
/// # Errors
/// [`SweepError::Empty`] when `n` is zero or there are no parameters, and
/// [`SweepError::InvalidBounds`] for a bad bound.
pub fn latin_hypercube(
    bounds: &[(f64, f64)],
    n: usize,
    seed: u64,
) -> Result<Vec<Vec<f64>>, SweepError> {
    check_bounds(bounds)?;
    if n == 0 {
        return Err(SweepError::Empty);
    }
    let mut rng = Rng::new(seed);
    let mut columns: Vec<Vec<f64>> = Vec::with_capacity(bounds.len());
    for &(lo, hi) in bounds {
        // Fisher–Yates shuffle of the stratum indices.
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = (rng.next_u64() % (i as u64 + 1)) as usize;
            order.swap(i, j);
        }
        let column = order
            .into_iter()
            .map(|stratum| {
                let u = (stratum as f64 + rng.next_f64()) / n as f64;
                lo + (hi - lo) * u
            })
            .collect();
        columns.push(column);
    }
    Ok((0..n)
        .map(|row| columns.iter().map(|c| c[row]).collect())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_varies_the_first_parameter_fastest() {
        let pts = grid(&[(0.0, 1.0), (10.0, 30.0)], &[2, 3]).unwrap();
        assert_eq!(pts.len(), 6);
        assert_eq!(pts[0], alloc::vec![0.0, 10.0]);
        assert_eq!(pts[1], alloc::vec![1.0, 10.0]);
        assert_eq!(pts[2], alloc::vec![0.0, 20.0]);
        assert_eq!(pts[5], alloc::vec![1.0, 30.0]);
    }

    #[test]
    fn grid_with_one_step_uses_the_midpoint() {
        let pts = grid(&[(2.0, 4.0)], &[1]).unwrap();
        assert_eq!(pts, alloc::vec![alloc::vec![3.0]]);
    }

    #[test]
    fn grid_rejects_bad_input() {
        assert_eq!(grid(&[], &[]), Err(SweepError::Empty));
        assert_eq!(grid(&[(0.0, 1.0)], &[0]), Err(SweepError::Empty));
        assert_eq!(grid(&[(1.0, 0.0)], &[2]), Err(SweepError::InvalidBounds));
        assert_eq!(grid(&[(0.0, 1.0)], &[2, 2]), Err(SweepError::InvalidBounds));
    }

    #[test]
    fn latin_hypercube_hits_every_stratum_once() {
        let bounds = [(0.0, 1.0), (-5.0, 5.0)];
        let n = 16;
        let pts = latin_hypercube(&bounds, n, 7).unwrap();
        assert_eq!(pts.len(), n);
        for (d, &(lo, hi)) in bounds.iter().enumerate() {
            let mut strata: Vec<usize> = pts
                .iter()
                .map(|p| {
                    assert!(p[d] >= lo && p[d] <= hi, "{p:?}");
                    (((p[d] - lo) / (hi - lo)) * n as f64)
                        .floor()
                        .min((n - 1) as f64) as usize
                })
                .collect();
            strata.sort_unstable();
            assert_eq!(strata, (0..n).collect::<Vec<_>>(), "dimension {d}");
        }
    }

    #[test]
    fn latin_hypercube_is_reproducible_for_a_seed() {
        let bounds = [(0.0, 1.0)];
        assert_eq!(
            latin_hypercube(&bounds, 8, 42).unwrap(),
            latin_hypercube(&bounds, 8, 42).unwrap()
        );
        assert_ne!(
            latin_hypercube(&bounds, 8, 42).unwrap(),
            latin_hypercube(&bounds, 8, 43).unwrap()
        );
    }

    /// Whether the first `2^m` points put exactly one point in every
    /// elementary box `[i/2^a, (i+1)/2^a) × [j/2^b, (j+1)/2^b)` with `a + b = m`.
    fn is_net(points: &[Vec<f64>], dims: (usize, usize), m: u32) -> bool {
        for a in 0..=m {
            let b = m - a;
            let mut counts = alloc::vec![0usize; 1usize << m];
            for p in points {
                let i = (p[dims.0] * (1u64 << a) as f64) as usize;
                let j = (p[dims.1] * (1u64 << b) as f64) as usize;
                counts[(i << b) + j] += 1;
            }
            if counts.iter().any(|&c| c != 1) {
                return false;
            }
        }
        true
    }

    #[test]
    fn sobol_first_block_is_stratified_in_every_dimension() {
        // Each 1-D projection of the first 2^4 points has one point per bin.
        let pts = sobol(&[(0.0, 1.0); 6], 16).unwrap();
        for d in 0..6 {
            let mut bins: Vec<usize> = pts.iter().map(|p| (p[d] * 16.0) as usize).collect();
            bins.sort_unstable();
            assert_eq!(bins, (0..16).collect::<Vec<_>>(), "dimension {d}");
        }
    }

    #[test]
    fn sobol_first_two_dimensions_form_a_net() {
        // The first two Sobol dimensions are a (0, 2)-sequence in base 2.
        let pts = sobol(&[(0.0, 1.0); 2], 16).unwrap();
        assert!(is_net(&pts, (0, 1), 4));
    }

    #[test]
    fn sobol_points_lie_in_the_box() {
        let pts = sobol(&[(-1.0, 3.0), (10.0, 20.0)], 64).unwrap();
        assert_eq!(pts.len(), 64);
        assert!(pts
            .iter()
            .all(|p| p[0] >= -1.0 && p[0] < 3.0 && p[1] >= 10.0 && p[1] < 20.0));
    }

    #[test]
    fn sobol_rejects_too_many_dimensions() {
        assert_eq!(
            sobol(&[(0.0, 1.0); 11], 4),
            Err(SweepError::UnsupportedDimension)
        );
    }

    #[test]
    fn latin_hypercube_rejects_zero_points() {
        assert_eq!(latin_hypercube(&[(0.0, 1.0)], 0, 1), Err(SweepError::Empty));
    }
}

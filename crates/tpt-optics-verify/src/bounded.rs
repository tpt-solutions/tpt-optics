//! A small bounded checker: exhaustive enumeration over a finite grid.
//!
//! This is the light alternative to a model checker. Each input is an
//! interval with a fixed number of evenly spaced samples, endpoints always
//! included. The property runs on every point of the product grid, in a fixed
//! order, so a failure is reproducible and reports the exact input.
//!
//! It is a finite check, not a proof. It says nothing between grid points, so
//! choose the grid to put samples where the risk is (near the edges of the
//! domain and near any singularity).

use alloc::string::String;

/// One input axis: a closed interval sampled at `steps` evenly spaced points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Axis {
    /// Lower end of the interval (sampled exactly).
    pub lo: f64,
    /// Upper end of the interval (sampled exactly).
    pub hi: f64,
    /// Number of samples, at least one. With one sample the axis is `lo`.
    pub steps: usize,
}

impl Axis {
    /// An axis over `[lo, hi]` with `steps` samples.
    #[must_use]
    pub const fn new(lo: f64, hi: f64, steps: usize) -> Self {
        Self { lo, hi, steps }
    }

    /// The `i`-th sample. The last one is exactly `hi`, not a rounded value.
    #[must_use]
    pub fn sample(&self, i: usize) -> f64 {
        if self.steps <= 1 {
            self.lo
        } else if i + 1 >= self.steps {
            self.hi
        } else {
            self.lo + (self.hi - self.lo) * (i as f64) / ((self.steps - 1) as f64)
        }
    }
}

/// The first grid point where a property fails.
#[derive(Clone, Debug, PartialEq)]
pub struct Counterexample<const N: usize> {
    /// The inputs at the failing point, in axis order.
    pub point: [f64; N],
    /// The message the property returned.
    pub message: String,
}

/// Run `prop` on every point of the grid spanned by `axes`.
///
/// The property returns `Ok(())` when the invariant holds at that point, or
/// `Err(message)` when it does not. A property that does not apply at some
/// point (a domain condition, say) should return `Ok(())` there.
///
/// # Returns
/// The number of points checked, or the first [`Counterexample`].
///
/// # Panics
/// If any axis has zero steps.
///
/// # Errors
/// [`Counterexample`] for the first failing point in the enumeration order,
/// which varies the first axis fastest.
pub fn check<const N: usize, F>(axes: [Axis; N], mut prop: F) -> Result<usize, Counterexample<N>>
where
    F: FnMut([f64; N]) -> Result<(), String>,
{
    assert!(
        axes.iter().all(|a| a.steps > 0),
        "each axis needs at least one sample"
    );
    let mut index = [0_usize; N];
    let mut cases = 0;
    loop {
        let point = core::array::from_fn(|k| axes[k].sample(index[k]));
        cases += 1;
        if let Err(message) = prop(point) {
            return Err(Counterexample { point, message });
        }
        // Advance the odometer: the first axis changes fastest.
        let mut k = 0;
        loop {
            if k == N {
                return Ok(cases);
            }
            index[k] += 1;
            if index[k] < axes[k].steps {
                break;
            }
            index[k] = 0;
            k += 1;
        }
    }
}

/// [`check`], panicking with the failing point when the invariant breaks.
///
/// For use in tests.
///
/// # Panics
/// On the first failing grid point, or if an axis has zero steps.
#[track_caller]
pub fn assert_holds<const N: usize, F>(axes: [Axis; N], prop: F)
where
    F: FnMut([f64; N]) -> Result<(), String>,
{
    if let Err(c) = check(axes, prop) {
        panic!("bounded check failed at {:?}: {}", c.point, c.message);
    }
}

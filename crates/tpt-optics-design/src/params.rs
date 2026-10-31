//! Named parameters for optimisation: a set of `name = value` pairs that
//! maps to the plain `f64` vector the optimisers use.
//!
//! Naming parameters makes a design file or an optimisation log readable.
//! Names must be unique, and lookups are exact.

use alloc::string::{String, ToString};
use alloc::vec::Vec;
use core::fmt;

use crate::optimize::{damped_least_squares_observed, Control, OptimizeResult, Progress};

/// Why a [`ParamSet`] operation failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ParamError {
    /// Two parameters share a name.
    DuplicateName(String),
    /// No parameter has this name.
    UnknownName(String),
    /// The value list does not match the number of parameters.
    LengthMismatch {
        /// Number of parameters in the set.
        expected: usize,
        /// Number of values given.
        found: usize,
    },
}

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateName(n) => write!(f, "parameter '{n}' is defined twice"),
            Self::UnknownName(n) => write!(f, "no parameter named '{n}'"),
            Self::LengthMismatch { expected, found } => {
                write!(f, "expected {expected} values, got {found}")
            }
        }
    }
}

impl core::error::Error for ParamError {}

/// Ordered, named parameter values.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ParamSet {
    names: Vec<String>,
    values: Vec<f64>,
}

impl ParamSet {
    /// Build a set from `(name, value)` pairs.
    ///
    /// # Errors
    /// [`ParamError::DuplicateName`] if a name appears twice.
    pub fn new(pairs: &[(&str, f64)]) -> Result<Self, ParamError> {
        let mut names: Vec<String> = Vec::with_capacity(pairs.len());
        for (name, _) in pairs {
            if names.iter().any(|n| n == name) {
                return Err(ParamError::DuplicateName((*name).to_string()));
            }
            names.push((*name).to_string());
        }
        Ok(ParamSet {
            names,
            values: pairs.iter().map(|(_, v)| *v).collect(),
        })
    }

    /// Value of a parameter, or `None` for an unknown name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<f64> {
        self.index_of(name).map(|i| self.values[i])
    }

    /// Set a parameter's value.
    ///
    /// # Errors
    /// [`ParamError::UnknownName`] if no parameter has this name.
    pub fn set(&mut self, name: &str, value: f64) -> Result<(), ParamError> {
        let i = self
            .index_of(name)
            .ok_or_else(|| ParamError::UnknownName(name.to_string()))?;
        self.values[i] = value;
        Ok(())
    }

    /// Parameter names, in order.
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Parameter values, in the same order as [`ParamSet::names`].
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// A copy of this set with new values, in the same order.
    ///
    /// # Errors
    /// [`ParamError::LengthMismatch`] if `values` has the wrong length.
    pub fn with_values(&self, values: &[f64]) -> Result<Self, ParamError> {
        if values.len() != self.values.len() {
            return Err(ParamError::LengthMismatch {
                expected: self.values.len(),
                found: values.len(),
            });
        }
        Ok(ParamSet {
            names: self.names.clone(),
            values: values.to_vec(),
        })
    }

    /// Position of a named parameter.
    #[must_use]
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }
}

/// Levenberg–Marquardt on a named parameter set.
///
/// `r(params, out)` writes the residuals for the named values. The returned
/// set holds the best values found, and the observer works as in
/// [`damped_least_squares_observed`].
pub fn optimize_named<R, C>(
    params: &ParamSet,
    r: &R,
    n_res: usize,
    iterations: usize,
    observer: &mut C,
) -> (ParamSet, OptimizeResult)
where
    R: Fn(&ParamSet, &mut Vec<f64>),
    C: FnMut(&Progress<'_>) -> Control,
{
    let template = params.clone();
    let wrapped = |x: &[f64], out: &mut Vec<f64>| {
        // Lengths match by construction, so this cannot fail.
        if let Ok(p) = template.with_values(x) {
            r(&p, out);
        }
    };
    let result =
        damped_least_squares_observed(&wrapped, params.values(), n_res, None, iterations, observer);
    let best = params
        .with_values(&result.x)
        .unwrap_or_else(|_| params.clone());
    (best, result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optimize::StopReason;

    #[test]
    fn names_resolve_to_values() {
        let mut p = ParamSet::new(&[("r1", 0.1), ("r2", -0.1)]).unwrap();
        assert_eq!(p.get("r2"), Some(-0.1));
        p.set("r1", 0.2).unwrap();
        assert_eq!(p.values(), &[0.2, -0.1]);
        assert_eq!(p.get("missing"), None);
    }

    #[test]
    fn duplicate_and_unknown_names_are_errors() {
        assert_eq!(
            ParamSet::new(&[("a", 1.0), ("a", 2.0)]),
            Err(ParamError::DuplicateName("a".into()))
        );
        let mut p = ParamSet::new(&[("a", 1.0)]).unwrap();
        assert_eq!(p.set("b", 0.0), Err(ParamError::UnknownName("b".into())));
    }

    #[test]
    fn with_values_checks_the_length() {
        let p = ParamSet::new(&[("a", 1.0), ("b", 2.0)]).unwrap();
        assert_eq!(
            p.with_values(&[1.0]),
            Err(ParamError::LengthMismatch {
                expected: 2,
                found: 1
            })
        );
    }

    #[test]
    fn optimize_named_recovers_named_values() {
        let p = ParamSet::new(&[("a", 0.0), ("b", 0.0)]).unwrap();
        let r = |q: &ParamSet, out: &mut Vec<f64>| {
            out.clear();
            out.push(q.get("a").unwrap() - 1.0);
            out.push(q.get("b").unwrap() + 2.0);
        };
        let (best, res) = optimize_named(&p, &r, 2, 50, &mut |_| Control::Continue);
        assert_eq!(res.stop, StopReason::Converged);
        assert!((best.get("a").unwrap() - 1.0).abs() < 1e-9);
        assert!((best.get("b").unwrap() + 2.0).abs() < 1e-9);
    }
}

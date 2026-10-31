# tpt-optics-design

Optical system optimization and tolerance analysis for
[tpt-optics](../../README.md). It provides lens merit functions built on
sequential ray tracing, operand-based optimization with exact gradients,
gradient-based and stochastic optimizers, Zernike wavefront-error and MTF
analysis, parameter sweeps, glass selection, and Monte-Carlo tolerancing with
compensators.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes (`alloc` required; see the `parallel` feature below)
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Merit functions.** RMS spot radius over a set of fractional ring
  weights, chromatic focal shift across a set of trace contexts, lateral
  colour, and distortion at a field angle.
- **Operands.** Each operand is one weighted residual. `optimize_system`
  drives a set of operands by damped least squares. Operands that read a
  traced ray have exact derivatives through the geometric crate's dual-number
  tracing (`gradient` module).
- **Damped least squares.** Levenberg–Marquardt-style optimizer over residual
  vectors, with optional box bounds, a numeric Jacobian when no analytic
  gradient is available, and an observer hook for progress reporting.
- **Stochastic optimizers.** Simulated annealing and a real-coded genetic
  algorithm for merits that are not smooth or not differentiable.
- **Named parameters.** `ParamSet` maps `name = value` pairs to the optimizer's
  parameter vector, so merits are written in terms of surfaces and gaps.
- **Wavefront merit.** Hexapolar pupil sampling, a least-squares fit of
  OSA/ANSI Zernike coefficients, and the residual RMS wavefront error. This
  replaces the earlier spot or OPD proxy as the wavefront merit.
- **MTF and Seidel.** System MTF computed from the wavefront, and the
  spherical-aberration Seidel sum from the exact on-axis wavefront.
- **Analysis and self-check.** `analyze` gives paraxial properties and, per
  field, the spot, ray fan, OPD, distortion and wavefront error, with CSV,
  JSON and SVG output. `self_check` lists invariants for a system.
- **Glass selection.** `rank_achromats` ranks crown–flint pairs for an
  achromatic doublet, and `genetic_achromat` searches the catalog with a
  genetic algorithm.
- **Sweeps.** Full-factorial grids, Sobol sequences and Latin-hypercube
  samples over a box of parameters.
- **Tolerance analysis.** Monte-Carlo sensitivity with normal, uniform and
  truncated-normal perturbations and an optional compensator. Reports the
  mean merit, its standard deviation, the 90th percentile and the failure
  fraction. `analyze_seeded` makes a run reproducible, and `sensitivity_table`
  reports the effect of each tolerance.

## Installation

```toml
[dependencies]
tpt-optics-design = "0.1"
```

Or enable the `design` feature of the umbrella crate
[`tpt-optics`](../tpt-optics). That feature also enables `geometric`.

## Quick start

Minimize the Rosenbrock residuals with damped least squares. The residual
closure writes one entry per residual into the vector:

```rust
use tpt_optics_design::optimize::damped_least_squares;

let residuals = |x: &[f64], r: &mut Vec<f64>| {
    r.clear();
    r.push(1.0 - x[0]);
    r.push(10.0 * (x[1] - x[0] * x[0]));
};

let result = damped_least_squares(&residuals, &[-1.2, 1.0], 2, None, 200);
assert!((result.x[0] - 1.0).abs() < 1e-3);
assert!((result.x[1] - 1.0).abs() < 1e-3);
```

A plain closure also works as an `Objective` for the stochastic optimizers,
because `Objective` is implemented for any `Fn(&[f64]) -> f64`.

## Feature flags

| Feature | Enables |
|---|---|
| `parallel` | `tolerance::analyze_parallel`, which runs Monte-Carlo passes on the rayon thread pool. Needs `std`. |
| `serde` | `Serialize` and `Deserialize` for the analysis report and self-check types. |

Both are off by default.

## Conventions

- Optimizers minimize their objective. Parameters are plain `f64` slices.
- `ParamBound { lo, hi }` gives box limits. `clamp` maps a value into them.
- Lengths are metres and angles are degrees for field angles, as in the
  `geometric` crate's merit inputs.

## Module map

| Module | Contents |
|---|---|
| `merit` | `rms_spot`, `spot_residuals`, `chromatic_focal_shift`, `lateral_colour`, `distortion` |
| `operands` | `Operand`, `RayAxis`, `optimize_system`, `residuals` |
| `gradient` | `optimize_exact`, `operand_gradient`, `system_jacobian`, `parameter_space`, `tolerance_budget` |
| `optimize` | `damped_least_squares` (and its `_observed` and `_with_jacobian` variants), `simulated_annealing`, `genetic_minimize`, `numeric_jacobian`, `Objective`, `ParamBound`, `OptimizeResult`, `StopReason` |
| `params` | `ParamSet`, `ParamError`, `optimize_named` |
| `wavefront` | `hexapolar_grid`, `wavefront_samples`, `fit_zernike`, `wavefront_rms`, `PupilSample` |
| `mtf` | `mtf_curves`, `system_mtf`, `MtfCurves` |
| `seidel` | `spherical_coefficient` |
| `analysis` | `analyze`, `AnalysisReport`, `FieldReport` |
| `check` | `self_check`, `Check` |
| `glassmap` | `rank_achromats`, `MapSpec`, `GlassPair`, `genetic_achromat` |
| `sweep` | `grid`, `sobol`, `latin_hypercube`, `SweepError` |
| `tolerance` | `Tolerance` (`Normal`, `Uniform`, `TruncatedNormal`), `ToleranceStats`, `analyze`, `analyze_seeded`, `analyze_parallel`, `sensitivity_table` |

## Accuracy and limits

- Gradient-based optimizers converge to local minima. Use the stochastic
  optimizers or several starting points when the merit is multimodal.
- Monte-Carlo tolerancing reports statistics from the number of passes you
  request. Increase the pass count when you need tight percentile estimates.
- The Zernike fit uses a fixed hexapolar sampling of the pupil. Its accuracy
  depends on the ring count and the Zernike order you choose.

## `no_std`

The crate is `#![no_std]` and builds on `tpt-optics-core`,
`tpt-optics-geometric` and `tpt-optics-wave`. Random sampling uses the
deterministic RNG from `tpt-optics-core`, so runs are reproducible from a seed.
The `parallel` feature brings in `std` through rayon.

## Testing

```sh
cargo test -p tpt-optics-design
cargo test -p tpt-optics-design --features parallel,serde
cargo build -p tpt-optics-design --no-default-features
```

Reference checks against a Cooke triplet live in `tests/cooke_triplet_reference.rs`.

## Related crates

- [`tpt-optics-geometric`](../tpt-optics-geometric): the sequential traces these merits evaluate
- [`tpt-optics-wave`](../tpt-optics-wave): the FFT behind the MTF
- [`tpt-optics-core`](../tpt-optics-core): the RNG and numeric types
- [`tpt-optics-io`](../tpt-optics-io): import existing designs from Zemax or CODE V as starting points

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

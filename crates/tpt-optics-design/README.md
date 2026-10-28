# tpt-optics-design

Optical system optimization and tolerance analysis for
[tpt-optics](../../README.md). It provides lens merit functions built on
sequential ray tracing, gradient-based and stochastic optimizers, a Zernike
wavefront-error fit, and Monte-Carlo tolerancing with compensators.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes (`alloc` required)
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Merit functions.** RMS spot radius over a set of fractional ring
  weights, chromatic focal shift across a set of trace contexts, and
  distortion at a field angle.
- **Damped least squares.** Levenberg–Marquardt-style optimizer over residual
  vectors, with optional box bounds and a numeric Jacobian when no analytic
  gradient is available.
- **Stochastic optimizers.** Simulated annealing and a real-coded genetic
  algorithm for merits that are not smooth or not differentiable.
- **Wavefront merit.** Hexapolar pupil sampling, a least-squares fit of
  OSA/ANSI Zernike coefficients, and the residual RMS wavefront error. This
  replaces the earlier spot or OPD proxy as the wavefront merit.
- **Tolerance analysis.** Monte-Carlo sensitivity with normal and uniform
  perturbations and an optional compensator. Reports the mean merit, its
  standard deviation, the 90th percentile and the failure fraction.

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

## Conventions

- Optimizers minimize their objective. Parameters are plain `f64` slices.
- `ParamBound { lo, hi }` gives box limits. `clamp` maps a value into them.
- Lengths are metres and angles are degrees for field angles, as in the
  `geometric` crate's merit inputs.

## Module map

| Module | Contents |
|---|---|
| `merit` | `rms_spot`, `spot_residuals`, `chromatic_focal_shift`, `distortion` |
| `optimize` | `damped_least_squares`, `simulated_annealing`, `genetic_minimize`, `numeric_jacobian`, `Objective`, `ParamBound`, `OptimizeResult` |
| `wavefront` | `hexapolar_grid`, `wavefront_samples`, `fit_zernike`, `wavefront_rms`, `PupilSample` |
| `tolerance` | `Tolerance`, `ToleranceStats`, `analyze` |

## Accuracy and limits

- Gradient-based optimizers converge to local minima. Use the stochastic
  optimizers or several starting points when the merit is multimodal.
- Monte-Carlo tolerancing reports statistics from the number of passes you
  request. Increase the pass count when you need tight percentile estimates.
- The Zernike fit uses a fixed hexapolar sampling of the pupil. Its accuracy
  depends on the ring count and the Zernike order you choose.

## `no_std`

The crate is `#![no_std]` and builds on `tpt-optics-core` and
`tpt-optics-geometric`. Random sampling uses the deterministic RNG from
`tpt-optics-core`, so runs are reproducible from a seed.

## Testing

```sh
cargo test -p tpt-optics-design
cargo build -p tpt-optics-design --no-default-features
```

## Related crates

- [`tpt-optics-geometric`](../tpt-optics-geometric): the sequential traces these merits evaluate
- [`tpt-optics-core`](../tpt-optics-core): the RNG and numeric types
- [`tpt-optics-io`](../tpt-optics-io): import existing designs from Zemax or CODE V as starting points

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

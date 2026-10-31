# tpt-optics-core

Foundation layer of the [tpt-optics](../../README.md) computational optics suite.

`tpt-optics-core` provides the unit-safe physical quantities, refractive index
representations, dispersion models and material database that every other
tpt-optics crate builds on. It also ships the numerical scaffolding those crates
share: complex numbers, dual numbers for forward-mode automatic
differentiation, small dense linear algebra, Bessel functions and a
deterministic random number generator.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes (uses `alloc` for tabulated data and material names)
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Unit-safe quantities.** `Wavelength`, `Frequency`, `OpticalPathLength`,
  `Angle` and `Kelvin` are newtypes over SI base storage, with constructors and
  accessors in the units optical engineers use (nm, µm, THz, degrees).
- **Dispersion models.** Sellmeier, Cauchy, Schott, Drude–Lorentz and
  tabulated (nearest and linear interpolation) models, all behind the
  `Dispersion` trait. Presets include BK7 and fused silica.
- **Material database.** Named `Material` entries with temperature
  dependence (`dn/dT`), a preset database, and Edlén air and vacuum models.
- **Refractive indices.** Real and complex (absorbing) indices, with
  absorption coefficients and normal-incidence reflectance.
- **Fresnel coefficients.** Amplitude and energy coefficients using the
  complex-angle formulation, plus Brewster and critical angles.
- **Numerics.** Complex numbers and dual numbers generic over the `Scalar`
  trait, 2×2 and 3-vector types (including the Gaussian-beam *q* law),
  Bessel functions of the first and second kind, and a xoshiro256\*\* RNG.

## Installation

```toml
[dependencies]
tpt-optics-core = "0.1"
```

Most users should depend on the umbrella crate [`tpt-optics`](../tpt-optics)
instead and enable the `core` feature.

## Quick start

```rust
use tpt_optics_core::prelude::*;

// Wavelength of the sodium d line.
let wl = Wavelength::from_nm(587.6);

// Refractive index of BK7 at that wavelength via the Sellmeier model.
let bk7 = Sellmeier::bk7();
let n = bk7.index(wl).re();
assert!((n - 1.5168).abs() < 1e-4);
```

Brewster and critical angles from the Fresnel module:

```rust
use tpt_optics_core::fresnel::{brewster_angle, critical_angle};

// Air to glass (n = 1.5): Brewster angle is atan(1.5).
let theta_b = brewster_angle(1.0, 1.5);
assert!((theta_b - 1.5_f64.atan()).abs() < 1e-12);

// Glass to air has a critical angle; air to glass does not.
assert!(critical_angle(1.5, 1.0).is_some());
assert!(critical_angle(1.0, 1.5).is_none());
```

## Module map

| Module | Contents |
|---|---|
| `units` | `Wavelength`, `Frequency`, `OpticalPathLength`, `Angle`, `Kelvin`, `Length`, `Curvature`, `Power`, physical constants |
| `dispersion` | `Dispersion` trait and the Sellmeier, Cauchy, Schott, Drude–Lorentz and tabulated models |
| `index` | `RefractiveIndex` and `ComplexRefractiveIndex` |
| `material` | `Material`, `MaterialDatabase`, `EdlenAir`, `Vacuum`, `ConstIndex` |
| `fresnel` | Fresnel amplitudes and energies, Brewster and critical angles |
| `complex` | `Complex<T>` |
| `dual` | `Dual` numbers for forward-mode automatic differentiation |
| `linalg` | `Vec2`, `Vec3`, `Matrix2` |
| `bessel` | Bessel functions J₀, J₁, Jₙ and K₀, K₁, Kₙ |
| `rng` | Deterministic xoshiro256\*\* generator |
| `num` | The `Scalar` trait shared by generic kernels |

Import the common items with `use tpt_optics_core::prelude::*;`.

## Accuracy and scope

The Bessel functions are from-scratch series and asymptotic implementations
with accuracy of roughly 1e-8 over the argument ranges the optics kernels use.
They are not general-purpose special-function libraries. The dispersion
formulas are the classic glass equations, defined in micrometres, and the
crate converts to and from `Wavelength` internally.

## `no_std`

The crate is `#![no_std]` and uses `alloc`. Transcendental functions come from
the pure-Rust [`libm`](https://crates.io/crates/libm) crate. Dependencies are
MIT or MIT/Apache-2.0 dual licensed.

## Testing

```sh
cargo test -p tpt-optics-core
cargo build -p tpt-optics-core --no-default-features
```

Property-based tests for the optical laws built on this crate live in
[`tpt-optics-verify`](../tpt-optics-verify).

## Related crates

- [`tpt-optics-geometric`](../tpt-optics-geometric): ray tracing built on these quantities
- [`tpt-optics-wave`](../tpt-optics-wave): diffraction and Gaussian beams
- [`tpt-optics-polarization`](../tpt-optics-polarization): Jones, Mueller and Stokes calculus
- [`tpt-optics-waveguide`](../tpt-optics-waveguide): mode solvers and photonic devices
- [`tpt-optics-design`](../tpt-optics-design): merit functions and optimizers
- [`tpt-optics-io`](../tpt-optics-io): Zemax and CODE V file formats and glass catalogs
- [`tpt-optics-verify`](../tpt-optics-verify): proptest strategies and invariants

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

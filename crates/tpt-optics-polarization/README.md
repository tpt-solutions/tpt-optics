# tpt-optics-polarization

Polarization optics for [tpt-optics](../../README.md): Jones and Mueller
calculus, Stokes vectors for partially polarized light, polarization ray
tracing through Fresnel interfaces, and birefringent, dichroic and optically
active media.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Jones calculus.** 2-component complex field vectors and 2×2 matrices for
  fully polarized light: polarizers, waveplates, rotators and diattenuators.
- **Stokes vectors.** Four-component intensity description that covers
  partially polarized light. Provides degree of polarization, azimuth and
  ellipticity.
- **Mueller calculus.** 4×4 real matrices that include depolarization. A
  coherency-matrix bridge converts a Jones matrix into a Mueller matrix.
- **Lu–Chipman decomposition.** Splits a Mueller matrix into a depolarizer, a
  retarder and a diattenuator, `M = M_Δ · M_R · M_D`.
- **Polarization ray tracing.** Tracks s and p reference frames through
  Fresnel interfaces, with per-element local operators and window transmission.
- **Crystal and medium models.** Uniaxial retardance, dichroic polarizers,
  optical activity (Biot's law rotation), walk-off and birefringence.

## Installation

```toml
[dependencies]
tpt-optics-polarization = "0.1"
```

Or enable the `polarization` feature of the umbrella crate
[`tpt-optics`](../tpt-optics).

## Quick start

Pass linearly polarized light through a quarter-wave plate at 45°. The output
is circularly polarized, so its Stokes degree of polarization is 1:

```rust
use std::f64::consts::FRAC_PI_4;
use tpt_optics_polarization::prelude::*;

let input = JonesVector::linear_x();
let qwp = JonesMatrix::quarter_wave(FRAC_PI_4);
let output = input.apply(&qwp);

assert!((output.intensity() - 1.0).abs() < 1e-12);

let stokes = StokesVector::from_jones(&output);
assert!((stokes.degree_of_polarization() - 1.0).abs() < 1e-12);
```

## Conventions

- Jones vectors are complex amplitudes in a fixed x–y basis, and angles are
  radians measured from x.
- Stokes vectors follow the `(S0, S1, S2, S3)` convention. `S0` is intensity.
- Retardances and rotations are in radians. Thicknesses and index differences
  are in metres and dimensionless units. `uniaxial_retardance` takes the
  wavelength as a `Wavelength` from `tpt-optics-core`.

## Module map

| Module | Contents |
|---|---|
| `jones` | `JonesVector`, `JonesMatrix`, standard elements (polarizer, waveplate, rotator, diattenuator) |
| `stokes` | `StokesVector` with degree of polarization, azimuth and ellipticity |
| `mueller` | `MuellerMatrix`, including the Jones-to-Mueller bridge and `apply_to_jones` |
| `decomposition` | `LuChipman`, `lu_chipman`: depolarizer, retarder and diattenuator of a Mueller matrix |
| `ray_trace` | `PolarizedElement`, `PolarizedSystem`, Fresnel-window reflection and transmission |
| `crystals` | `uniaxial_retardance`, `retarder_thickness`, `extraordinary_walkoff`, `dichroic_polarizer`, `optical_activity` |

Import the common items with `use tpt_optics_polarization::prelude::*;`.

## Scope

Jones calculus describes fully polarized light only. Use the Stokes and
Mueller forms for partially polarized or depolarized beams. Polarization ray
tracing applies Fresnel coefficients at each interface, so check the window
model against your own stack before relying on it for thin-film work.

## Testing

```sh
cargo test -p tpt-optics-polarization
cargo build -p tpt-optics-polarization --no-default-features
```

The Fresnel energy-conservation invariant is checked by property tests in
[`tpt-optics-verify`](../tpt-optics-verify).

## Related crates

- [`tpt-optics-core`](../tpt-optics-core): Fresnel coefficients and the `Wavelength` type
- [`tpt-optics-geometric`](../tpt-optics-geometric): the ray paths that polarization ray tracing follows
- [`tpt-optics-wave`](../tpt-optics-wave): scalar and coherent field propagation

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

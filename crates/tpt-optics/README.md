# tpt-optics

Umbrella crate for [tpt-optics](../../README.md), a pure-Rust computational
optics and photonics suite. It re-exports the eight subcrates behind Cargo
features, so you compile only the parts you use.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes. The umbrella crate is `#![no_std]`, and subcrates that
  need the standard library, such as `tpt-optics-io`, bring it in only when
  their feature is enabled.
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## What is included

- Geometric ray tracing, sequential and non-sequential
- Wave optics: FFT, diffraction, Gaussian beams, beam propagation, PSF and MTF
- Polarization: Jones, Mueller and Stokes calculus
- Photonics: slab and fiber modes, effective index method, coupled-mode theory,
  photonic crystals, integrated devices, a scalar FEM mode solver
- Optical design: merit functions, damped least squares, simulated annealing,
  genetic optimization and Monte-Carlo tolerancing
- File formats: Zemax `.zmx` text, CODE V `.seq`, AGF glass catalogs
- Verification: proptest strategies and optical invariants

Everything is written from scratch in Rust, with no C or C++ FFI and no
external optical solver.

## Installation

```toml
[dependencies]
tpt-optics = { version = "0.1", features = ["full"] }
```

Or select only what you need:

```toml
[dependencies]
tpt-optics = { version = "0.1", default-features = false, features = ["core", "polarization"] }
```

## Features

| Feature | Subcrate | Enables |
|---|---|---|
| `core` (default) | [`tpt-optics-core`](../tpt-optics-core) | Units, dispersion, materials, Fresnel, complex and dual numbers, Bessel functions |
| `geometric` (default) | [`tpt-optics-geometric`](../tpt-optics-geometric) | Sequential and non-sequential tracing, surfaces, ABCD matrices |
| `wave` | [`tpt-optics-wave`](../tpt-optics-wave) | FFT, diffraction, Gaussian beams, BPM, PSF/OTF/MTF, coherence. Also enables `geometric`. |
| `polarization` | [`tpt-optics-polarization`](../tpt-optics-polarization) | Jones, Mueller and Stokes calculus, polarization ray tracing, crystals |
| `waveguide` | [`tpt-optics-waveguide`](../tpt-optics-waveguide) | Slab and fiber modes, EIM, coupled-mode theory, photonic crystals, FEM |
| `design` | [`tpt-optics-design`](../tpt-optics-design) | Merit functions, optimizers, tolerancing, wavefront fit. Also enables `geometric`. |
| `io` | [`tpt-optics-io`](../tpt-optics-io) | Zemax `.zmx`, CODE V `.seq`, AGF glass catalogs. Also enables `geometric`. |
| `verify` | [`tpt-optics-verify`](../tpt-optics-verify) | proptest strategies and optical invariants. Also enables `geometric` and `polarization`. |
| `full` | all of the above | Every feature |

Each subcrate is re-exported under its own name, for example
`tpt_optics::geometric` or `tpt_optics::design`.

## Quick start

With the default features, the prelude brings in the core and geometric types:

```rust
use tpt_optics::prelude::*;

// Wavelength of the helium-neon laser line.
let wl = Wavelength::from_nm(632.8);

// Refractive index of BK7 at that wavelength via the Sellmeier model.
let bk7 = Sellmeier::bk7();
let n = bk7.index(wl).re();
assert!((n - 1.5151).abs() < 1e-4);
```

The prelude re-exports the common items of the `core`, `geometric`,
`polarization`, `wave` and `waveguide` subcrates. The `design`, `io` and
`verify` crates have no prelude. Use their module paths, for example
`tpt_optics::design::optimize::damped_least_squares`.

## `no_std`

The umbrella crate is `#![no_std]` and `#![forbid(unsafe_code)]`. The
`core`, `geometric`, `polarization`, `wave`, `waveguide` and `design` features
build without the standard library. The `io` and `verify` features require
`std`.

## Development

From the repository root, the `justfile` recipes run the same gate as CI:

```sh
just check   # fmt, clippy, tests, cargo-deny and the no_std gate
just test    # tests only
```

## Related documentation

- [Repository README](../../README.md): the full suite overview and crate table
- [CHANGELOG](CHANGELOG.md): release history for this crate
- [Specification](../../spec.txt): design goals and the crate layout

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

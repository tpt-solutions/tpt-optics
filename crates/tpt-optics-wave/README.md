# tpt-optics-wave

Wave optics for [tpt-optics](../../README.md): scalar diffraction, Fourier
optics, Gaussian beams, beam propagation and partial coherence. Everything
rests on a from-scratch radix-2 FFT.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes (`alloc` required)
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Radix-2 FFT.** Complex 1-D and 2-D transforms in both directions.
- **Scalar diffraction.** Angular-spectrum propagation, single-FFT Fresnel
  propagation and Fraunhofer (far-field) propagation.
- **Gaussian beams.** q-parameter beam representation, the ABCD law for
  beam transformation, Rayleigh range and thin-lens mode matching.
- **Fourier optics.** Lens Fourier transforms, and PSF, OTF and MTF computed
  from a pupil. The circular-pupil Hopkins reference is included.
- **Beam propagation method.** Split-step BPM in 2-D and 3-D, with
  configurable transverse and longitudinal grids.
- **Coherence theory.** Gaussian Schell-model mutual coherence and
  van Cittert–Zernike far-field coherence.

## Installation

```toml
[dependencies]
tpt-optics-wave = "0.1"
```

Or enable the `wave` feature of the umbrella crate
[`tpt-optics`](../tpt-optics). That feature also enables `geometric`.

## Quick start

Propagate a Gaussian beam by a Rayleigh range and check its radius grows by
√2. A tiny example that uses only the Gaussian-beam helpers:

```rust
use tpt_optics_core::prelude::Wavelength;
use tpt_optics_wave::gaussian::rayleigh_range;

let wl = Wavelength::from_nm(632.8);
let w0 = 1.0e-3; // 1 mm waist radius

// Rayleigh range z_R = π n w0² / λ in air.
let z_r = rayleigh_range(w0, wl, 1.0);
assert!(z_r > 4.0 && z_r < 5.0);
```

Two-dimensional FFT round trip on a power-of-two grid. This example also
needs `tpt-optics-core` in your `[dependencies]` for `Complex`:

```rust
use tpt_optics_core::prelude::Complex;
use tpt_optics_wave::fft::{fft2, ifft2};

let n = 8;
let mut field: Vec<Complex<f64>> = (0..n * n)
    .map(|i| Complex::real(i as f64))
    .collect();
let original = field.clone();

fft2(&mut field, n, n, false);
ifft2(&mut field, n, n);

for (a, b) in field.iter().zip(original.iter()) {
    assert!(a.approx_eq(*b, 1e-9));
}
```

## Grid conventions

Field grids are square, power-of-two samples with uniform spacing. Lengths are
metres. Functions that take a grid length `n` or dimensions `n × m` assert
that the dimensions are powers of two.

## Module map

| Module | Contents |
|---|---|
| `fft` | Radix-2 complex FFT, 1-D and 2-D, forward and inverse |
| `diffraction` | `angular_spectrum`, `fresnel_sfft`, `fraunhofer` |
| `gaussian` | `GaussianBeam`, `rayleigh_range`, `mode_match_thin_lens` |
| `fourier` | `lens_fourier_transform`, `psf_from_pupil`, `mtf_from_psf` |
| `bpm` | `bpm_2d`, `bpm_3d`, `Bpm2dConfig`, `Bpm3dConfig` |
| `coherence` | `gaussian_schell_mu`, `van_cittert_zernike` |

Import the common items with `use tpt_optics_wave::prelude::*;`.

## Accuracy and limits

- The FFT and propagators are sampled methods. Accuracy depends on grid
  resolution, aliasing and the domain size, so check convergence for your
  problem.
- BPM is a paraxial, slowly-varying-envelope method. It does not model
  strongly backscattered light.
- The Hopkins PSF/MTF reference is for a circular, unobscured pupil.

## `no_std`

The crate is `#![no_std]` and depends on `tpt-optics-core` and
`tpt-optics-geometric`, both of which are `no_std`. Transcendental functions
come from [`libm`](https://crates.io/crates/libm).

## Testing

```sh
cargo test -p tpt-optics-wave
cargo build -p tpt-optics-wave --no-default-features
```

Integration tests live in `tests/integration.rs`.

## Related crates

- [`tpt-optics-core`](../tpt-optics-core): units, complex numbers and Bessel functions
- [`tpt-optics-geometric`](../tpt-optics-geometric): ray-level optics and ABCD matrices
- [`tpt-optics-polarization`](../tpt-optics-polarization): polarization state for vector fields
- [`tpt-optics-waveguide`](../tpt-optics-waveguide): guided-wave mode solving

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

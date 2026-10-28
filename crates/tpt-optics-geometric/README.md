# tpt-optics-geometric

Geometric optics for [tpt-optics](../../README.md): ray tracing through
lenses, mirrors and freeform surfaces, ABCD paraxial analysis, and
non-sequential Monte-Carlo stray-light tracing.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes (`alloc` required)
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Sequential ray tracing**, both exact and paraxial, through refracting,
  reflecting and diffracting surfaces. Aperture checks and optical path length
  accumulation are included.
- **Surface library.** Plane, sphere, conic and even asphere, Zernike freeform
  and polynomial-phase diffractive surfaces. Generic surfaces are intersected
  with Newton–Raphson refinement.
- **ABCD matrix formalism** for paraxial systems: free space, thin lenses,
  refracting surfaces, mirrors and telescopes. Covers Gaussian beam
  propagation and resonator stability.
- **Non-sequential Monte-Carlo tracing** with specular, Lambertian and
  Henyey–Greenstein scattering, for stray-light and illumination analysis.
- **Differentiable tracing.** Every trace kernel is generic over `Scalar`.
  Seeding a `Dual` parameter gives exact first derivatives, which the design
  crate uses for gradient-based optimization.
- **Global placement helpers** in the `frame` module.

## Installation

```toml
[dependencies]
tpt-optics-geometric = "0.1"
```

Or enable the `geometric` feature (on by default) of the umbrella crate
[`tpt-optics`](../tpt-optics).

## Quick start

Paraxial imaging with ABCD matrices. Light travels 150 mm through free space
and then meets a thin lens of focal length 100 mm. Compose the elements with
`mul_matrix`, where the element on the right acts first:

```rust
use tpt_optics_geometric::abcd::{free_space, thin_lens};

let system = thin_lens::<f64>(0.100).mul_matrix(free_space::<f64>(0.150));
```

Refraction through a single spherical surface uses the same `Matrix2` type and
is available as `refracting_surface(n1, n2, curvature)`.

## Units and conventions

Lengths are metres and angles are radians. A `SequentialSystem` is described
in its own local frame. Surface 0 sits at z = 0, and each surface's vertex lies
`thickness` along +z from the previous vertex. Use the `frame` module to place
systems in a global frame.

## Module map

| Module | Contents |
|---|---|
| `sequential` | `SequentialSystem`, `SeqSurface`, `SurfaceAction`, `Medium`, `TraceContext`, `ParaxialSystem` |
| `surface` | `Plane`, `Sphere`, `Asphere`, `Zernike` and the `OpticalSurface` trait |
| `ray` | `Ray<T>` with normalisation, propagation and frame transforms |
| `snell` | `refract`, `reflect` and `diffract` vector laws |
| `abcd` | ABCD matrices for free space, thin lenses, interfaces and mirrors |
| `nseq` | `NsScene`, `NsObject` and `ScatterModel` for non-sequential tracing |
| `differentiable` | Dual-number helpers for gradient computation through a trace |
| `frame` | Local to global coordinate frames |

Import the common items with `use tpt_optics_geometric::prelude::*;`.

## Testing

```sh
cargo test -p tpt-optics-geometric
cargo build -p tpt-optics-geometric --no-default-features
```

Property tests (Snell reversibility, trace reversibility, ABCD symplecticity)
live in `tests/properties.rs`. The shared invariants are defined in
[`tpt-optics-verify`](../tpt-optics-verify).

## Related crates

- [`tpt-optics-core`](../tpt-optics-core): units, materials and dispersion used by every surface
- [`tpt-optics-design`](../tpt-optics-design): merit functions that trace through these systems
- [`tpt-optics-wave`](../tpt-optics-wave): wave-optics counterparts to ray-level Gaussian beams
- [`tpt-optics-io`](../tpt-optics-io): import and export of Zemax and CODE V lens prescriptions

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

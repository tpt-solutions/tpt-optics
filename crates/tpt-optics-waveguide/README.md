# tpt-optics-waveguide

Photonic waveguide modeling for [tpt-optics](../../README.md): slab and fiber
mode solving, the effective index method, coupled-mode theory, photonic
crystal band structure, integrated devices and a scalar finite-element mode
solver.

- **License:** MIT OR Apache-2.0
- **`no_std`:** yes (`alloc` required)
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Slab waveguides.** Symmetric slab TE and TM modes from an exact
  transcendental solve, plus planar-stack transfer matrices for multilayer
  stacks.
- **Fiber modes.** Weak-guidance LP modes of step-index fibers, and of
  power-law graded-index fibers.
- **Effective index method.** Reduces a 2-D channel waveguide to two slab
  problems to estimate its effective index.
- **Coupled-mode theory.** Directional couplers, Bragg gratings and ring
  resonators, with coupler power, grating reflectivity and ring resonance and
  free-spectral-range helpers.
- **Photonic crystals.** 1-D Bloch band structure of a Bragg stack, and 2-D
  plane-wave-expansion TE bands.
- **Integrated devices.** Mach–Zehnder transmission and arrayed-waveguide
  grating (AWG) response.
- **Scalar FEM mode solver.** A 2-D finite-element eigenmode solver on
  rectangular domains with P1 triangles, for the fundamental mode of an
  arbitrary index profile.

## Installation

```toml
[dependencies]
tpt-optics-waveguide = "0.1"
```

Or enable the `waveguide` feature of the umbrella crate
[`tpt-optics`](../tpt-optics).

## Quick start

Find the TE modes of a 4 µm slab with core index 1.5 and cladding index 1.45
at 1.55 µm:

```rust
use tpt_optics_core::prelude::Wavelength;
use tpt_optics_waveguide::prelude::*;
use tpt_optics_waveguide::slab::SlabPolarization;

let slab = SlabLayer {
    n_core: 1.5,
    n_clad: 1.45,
    thickness: 4.0e-6,
};
let modes = slab_modes(&slab, Wavelength::from_um(1.55), SlabPolarization::Te);

// The fundamental mode is guided, so its index lies between the cladding
// and core indices.
let n_eff = modes[0].n_eff;
assert!(n_eff > 1.45 && n_eff < 1.5);
```

This example also needs `tpt-optics-core` in your `[dependencies]` for
`Wavelength`.

## Conventions

- `n_eff` is the effective index β/k₀. A mode is guided when its `n_eff`
  lies between the cladding and core indices.
- Lengths are metres and wavelengths are `Wavelength` values from
  `tpt-optics-core`.
- The fiber and FEM solvers are weak-guidance or scalar approximations, so
  they are not full vector solutions. See the accuracy notes below.

## Module map

| Module | Contents |
|---|---|
| `slab` | `SlabLayer`, `SlabMode`, `SlabPolarization`, `slab_modes`, planar-stack transfer matrices |
| `fiber` | `LpMode`, `fiber_lp_modes` (step-index) |
| `graded` | `graded_lp_modes` (power-law graded-index) |
| `eim` | `channel_wguide_neff`, the effective index method |
| `coupled` | `coupler_power`, `grating_reflectivity`, `ring_resonances`, `ring_fsr` |
| `phc` | `bragg_stack_bands` (1-D Bloch), `pwem_te_bands` (2-D plane-wave expansion) |
| `devices` | `mzi_transmission`, `awg_response` |
| `fem` | `RectDomain`, `fem_fundamental_mode` (scalar P1 finite elements) |

Import the common items with `use tpt_optics_waveguide::prelude::*;`.

## Accuracy and limits

- **Weak guidance.** The LP-mode solvers assume small index contrast.
  High-contrast guides, such as silicon wires, need a vector solver.
- **Scalar FEM.** The FEM solver is a scalar weak-guidance formulation. A
  vectorial upgrade is reserved for the `tpt-fem` dependency named in
  `spec.txt`.
- **Plane-wave cutoffs.** Photonic crystal bands depend on the number of
  plane waves kept, so check convergence.
- **Device models.** The integrated-device helpers are response models, not
  full propagation simulations. Check the function documentation for the
  assumptions each one makes.

## `no_std`

The crate is `#![no_std]`. The `std` feature exists for future std-only helpers
but currently enables nothing.

## Testing

```sh
cargo test -p tpt-optics-waveguide
cargo build -p tpt-optics-waveguide --no-default-features
```

## Related crates

- [`tpt-optics-core`](../tpt-optics-core): wavelengths, refractive indices and dispersion
- [`tpt-optics-wave`](../tpt-optics-wave): free-space and beam propagation
- [`tpt-optics-polarization`](../tpt-optics-polarization): polarization of guided modes

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

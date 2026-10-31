# tpt-optics-waveguide

Photonic waveguide modeling for [tpt-optics](../../README.md): slab and fiber
mode solving, the effective index method, coupled-mode theory, photonic
crystal band structure and defects, thin-film stacks, integrated devices, and
scalar and full-vector finite-element mode solvers.

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
  problems to estimate its effective index, for TE and TM polarizations.
- **Coupled-mode theory.** Directional couplers, Bragg gratings and ring
  resonators, with coupler power, grating reflectivity, ring resonances,
  free spectral range, Q and finesse helpers.
- **Photonic crystals.** 1-D Bloch band structure of a Bragg stack, 2-D
  plane-wave-expansion TE bands, and point and line defects in 2-D lattices
  through supercell expansion.
- **Thin-film stacks.** Characteristic-matrix response at oblique incidence,
  for s and p polarization, with complex indices so absorbing layers are
  handled. For a lossless stack, reflectance plus transmittance is one.
- **Integrated devices.** Mach–Zehnder transmission and arrayed-waveguide
  grating (AWG) response.
- **Scalar FEM mode solver.** A 2-D finite-element eigenmode solver on
  rectangular domains with P1 triangles, for the fundamental mode of an
  arbitrary index profile.
- **Full-vector FEM mode solver.** A hybrid edge–nodal solver (first-order
  Nédélec edges for the transverse field, P1 nodes for the longitudinal field)
  that returns the effective index of a vector mode, with convergence and
  residual reporting.
- **Presets.** A single-mode step fiber (SMF-28 values), a silicon ring
  resonator, and a photonic-crystal line defect.

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
- Solver families differ in approximation. The `fiber`, `graded` and scalar
  `fem` solvers are weak-guidance or scalar. `vector_fem` is the full-vector
  solver. Pick the one that matches your index contrast. See the accuracy
  notes below.

## Module map

| Module | Contents |
|---|---|
| `slab` | `SlabLayer`, `SlabMode`, `SlabPolarization`, `slab_modes`, planar-stack transfer matrices (`StackLayer`, `stack_r_t`) |
| `fiber` | `LpMode`, `fiber_lp_modes` (step-index) |
| `graded` | `graded_lp_modes` (power-law graded-index) |
| `eim` | `ChannelGuide`, `channel_wguide_neff`, `channel_wguide_neff_pol`, the effective index method |
| `coupled` | `coupler_power`, `coupling_length`, `grating_reflectivity`, `bragg_wavelength`, `ring_resonances`, `ring_fsr`, `ring_q`, `ring_finesse`, `ring_add_drop` |
| `phc` | `bragg_stack_bands` (1-D Bloch), `pwem_te_bands` (2-D plane-wave expansion) |
| `supercell` | `Rod`, `point_defect_rods`, `supercell_te_freqs`, `supercell_tm_freqs`, `bulk_gap_from_supercell` |
| `thin_film` | `FilmLayer`, `FilmPower`, `film_response`, `quarter_wave_layer` |
| `devices` | `mzi_transmission`, `mzi_phase`, `mzi_fsr`, `awg_response`, `awg_channels` |
| `fem` | `RectDomain`, `fem_fundamental_mode` (scalar P1 finite elements) |
| `vector_fem` | `VectorMode`, `vector_mode_neff` (hybrid edge–nodal full-vector FEM) |
| `eigen` | `symmetric_eigenvalues`, a dense symmetric eigenvalue solver (Householder tridiagonalisation, then the implicit QL algorithm) |
| `presets` | `StepFiber`, `SMF28`, `fiber_modes`, `SiliconRing`, `pc_line_defect_rods` |

Import the common items with `use tpt_optics_waveguide::prelude::*;`. The
prelude covers the coupled-mode helpers, the device responses, the EIM, the
fiber and graded solvers, the photonic-crystal bands, the slab solver and
the FEM entry points. Other modules are reached by path.

## Accuracy and limits

- **Weak guidance.** The LP-mode solvers (`fiber`, `graded`) assume small
  index contrast. High-contrast guides, such as silicon wires, need
  `vector_fem`.
- **Scalar FEM.** `fem` is a scalar weak-guidance formulation. Use
  `vector_fem` when polarization or high contrast matters. Check its
  `converged` and `residual` fields before relying on a result. An evanescent
  solution is returned as an error, not as an index.
- **Plane-wave cutoffs.** Photonic crystal bands and supercell defect modes
  depend on the number of plane waves kept, so check convergence.
- **Thin-film model.** The stack is planar and uniform in each layer. It does
  not model interface roughness.
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

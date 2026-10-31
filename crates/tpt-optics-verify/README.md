# tpt-optics-verify

Verification harnesses for [tpt-optics](../../README.md). The crate collects
the optical laws the rest of the suite must obey as checkable predicates,
reusable proptest strategies for generating physical inputs, a bounded grid
checker, and a system-level check that runs the invariants on a whole
sequential lens system. The critical algorithms (ray–sphere intersection,
lossless energy split) are checked by property tests over random inputs.

- **License:** MIT OR Apache-2.0
- **`no_std`:** the library is `#![no_std]`. It depends on proptest with its
  default features, which pulls in `std`.
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Strategies.** proptest generators for bounded floats, unit incident
  directions, rays, refractive indices, index pairs in descending order,
  surface curvatures, conic constants, unit-modulus complex numbers and
  Sellmeier coefficient sets.
- **Invariants.** Optical laws as predicates that return `Result<(), String>`,
  so a failing case explains itself:
  - Snell's law is reversible, and conserves tangential momentum
  - Fresnel energy is conserved (R + T = 1 for lossless interfaces)
  - ABCD matrices are symplectic, and the determinant of a refraction matrix
    equals the index ratio
  - Ray traces through a system are reversible
  - Optical path length is invariant under frame rotation
- **Bounded checks.** `bounded::check` enumerates a finite grid, always
  including the endpoints, and reports the first counterexample. It is a
  light stand-in for a model checker. It checks the grid you give it, not
  every input.
- **System checks.** `system_check::verify_system` runs the trace invariants
  (reversibility, OPL frame invariance) over a whole `SequentialSystem` for a
  set of rays, and returns one `InvariantCheck` per check.

## Installation

Add the crate as a development dependency in the crate you want to verify:

```toml
[dev-dependencies]
tpt-optics-verify = "0.1"
```

## Quick start

Call an invariant directly. It returns `Ok(())` when the law holds within the
tolerance you give:

```rust
use tpt_optics_verify::invariants::fresnel_energy_conserved;

// Air to glass at 30 degrees incidence: reflected plus transmitted energy is 1.
assert!(fresnel_energy_conserved(1.0, 1.5, 0.5236, 1e-9).is_ok());
```

Use a strategy in a proptest to check the same law over many random inputs:

```rust
use proptest::prelude::*;
use tpt_optics_verify::{invariants::fresnel_energy_conserved, strategies};

proptest! {
    #[test]
    fn fresnel_energy_holds(
        (n1, n2) in strategies::index_pair_descending(),
        theta in strategies::f64_in(0.0, 1.2),
    ) {
        prop_assert!(fresnel_energy_conserved(n1, n2, theta, 1e-9).is_ok());
    }
}
```

## Module map

| Module | Contents |
|---|---|
| `strategies` | proptest generators: `f64_in`, `incident_direction`, `incident_ray`, `refractive_index`, `index_pair_descending`, `curvature`, `conic_constant`, `unit_complex`, `sellmeier_coeffs` |
| `invariants` | `snell_reversible`, `snell_tangential_momentum`, `abcd_symplectic`, `abcd_refraction_det`, `trace_reversibility`, `fresnel_energy_conserved`, `opl_rotation_invariance`, and `d_line` (the 587.6 nm reference wavelength) |
| `bounded` | `check`, `assert_holds`, `Axis`, `Counterexample`: exhaustive grid checks with endpoints always sampled |
| `system_check` | `verify_system`, `InvariantCheck`: the trace invariants run on a whole sequential system |

## Running the checks

```sh
cargo test -p tpt-optics-verify
```

- `tests/properties.rs` exercises every invariant with randomized inputs,
  including ray–sphere intersection and the lossless energy split. Proptest
  regression seeds are checked in alongside it.
- `tests/bounded.rs` runs the grid checker on its own behaviour.
- `tests/system_cases.rs` runs `verify_system` on sequential systems beyond the
  singlet, such as a mirror and a glass plate.

## Choosing tolerances

Each invariant takes a tolerance. Use a value near the floating-point
precision of the operation being checked, for example `1e-9` for a single
Fresnel evaluation. Loosen it only when a longer chain of operations
accumulates rounding error, and say so in the test.

## Related crates

- [`tpt-optics-core`](../tpt-optics-core): the physical quantities and models being checked
- [`tpt-optics-geometric`](../tpt-optics-geometric): ray tracing and ABCD matrices under test
- [`tpt-optics-polarization`](../tpt-optics-polarization): Fresnel and polarization calculus under test

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

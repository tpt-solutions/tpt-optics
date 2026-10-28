# Changelog

All notable changes to `tpt-optics-waveguide` are documented here. Format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the photonics crate.

### Added

- Symmetric slab TE and TM modes, and planar-stack transfer matrices.
- Weak-guidance LP modes of step-index fibers, and of power-law graded-index
  fibers (`graded` module).
- Effective index method for channel waveguides.
- Coupled-mode theory for directional couplers, gratings and ring resonators.
- 1-D Bloch band structure and 2-D plane-wave-expansion photonic-crystal bands.
- Mach–Zehnder interferometer and arrayed-waveguide grating device models.
- Scalar 2-D finite-element mode solver on P1 triangles.
- `no_std` + `alloc` support, with an empty `std` feature reserved for future
  std-only helpers.

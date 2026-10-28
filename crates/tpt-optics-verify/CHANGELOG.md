# Changelog

All notable changes to `tpt-optics-verify` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the verification crate.

### Added

- Shared proptest strategies for bounded floats, incident directions and rays,
  refractive indices, index pairs, surface curvatures, conic constants,
  unit-modulus complex numbers and Sellmeier coefficients.
- Optical invariants as checkable predicates: Snell reversibility and momentum,
  Fresnel energy conservation, ABCD symplecticity, refraction determinant
  equal to the index ratio, trace reversibility, and OPL rotation invariance.
- Kani bounded-model-checking harnesses for nearest-positive-root intersection
  and lossless beam-splitter energy conservation, compiled under `cfg(kani)`.
- Property tests in `tests/properties.rs`, with checked-in regression seeds.

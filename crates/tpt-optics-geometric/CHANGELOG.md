# Changelog

All notable changes to `tpt-optics-geometric` are documented here. Format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the geometric optics crate.

### Added

- Sequential ray tracing, exact and paraxial, through refracting, reflecting
  and diffracting surfaces, with aperture checks and optical path length.
- Surface library: plane, sphere, conic and even asphere, Zernike freeform,
  and polynomial-phase diffractive surfaces, with Newton–Raphson
  intersection for generic surfaces.
- Entrance pupil and paraxial system analysis.
- ABCD matrix formalism for free space, thin lenses, refracting surfaces,
  mirrors and telescopes, including Gaussian beam propagation.
- Non-sequential Monte-Carlo tracing with specular, Lambertian and
  Henyey–Greenstein scattering, and stray-light reporting.
- Differentiable tracing: trace kernels are generic over `Scalar`, so seeding a
  `Dual` parameter gives exact first derivatives.
- Global placement helpers in the `frame` module.
- `no_std` + `alloc` support.

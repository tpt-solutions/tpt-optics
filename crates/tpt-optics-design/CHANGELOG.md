# Changelog

All notable changes to `tpt-optics-design` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Zernike wavefront-error merit (`wavefront` module). `wavefront_samples`
  traces a hexapolar pupil grid and references each ray's OPD to a plane
  through the chief image point. `fit_zernike` least-squares fits OSA/ANSI
  Zernike coefficients, and `wavefront_rms` gives the residual RMS error. It
  replaces the spot/OPD proxy as the wavefront merit.
- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the optical design crate.

### Added

- Merit functions: RMS spot radius, chromatic focal shift, and distortion.
- Damped least squares (Levenberg–Marquardt style), with a numeric Jacobian
  fallback and optional box bounds.
- Simulated annealing and a real-coded genetic optimizer.
- Monte-Carlo tolerancing with normal and uniform perturbations and
  compensators.
- `no_std` + `alloc` support.

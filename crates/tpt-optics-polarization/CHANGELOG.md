# Changelog

All notable changes to `tpt-optics-polarization` are documented here. Format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the polarization crate.

### Added

- Jones calculus: vectors, matrices, polarizers, waveplates, rotators and
  diattenuators.
- Stokes vectors with degree of polarization, azimuth and ellipticity.
- Mueller calculus, with a coherency-matrix bridge from Jones to Mueller.
- Polarization ray tracing with s and p frames through Fresnel interfaces,
  retarders, and walk-off.
- Dichroism and optical activity.
- `no_std` + `alloc` support.

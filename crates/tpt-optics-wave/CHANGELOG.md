# Changelog

All notable changes to `tpt-optics-wave` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the wave optics crate.

### Added

- Radix-2 FFT, 1-D and 2-D, forward and inverse.
- Angular-spectrum, Fresnel and Fraunhofer propagation.
- Gaussian q-parameter beams, mode matching through a thin lens, and Rayleigh
  range.
- Lens Fourier transform, and PSF, OTF and MTF from a pupil, with the
  Hopkins circular-pupil reference.
- Split-step beam propagation method in 2-D and 3-D.
- Gaussian Schell-model mutual coherence and van Cittert–Zernike coherence.
- `no_std` + `alloc` support.

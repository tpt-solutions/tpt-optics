# Changelog

All notable changes to `tpt-optics-core` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the foundation layer.

### Added

- Unit-safe quantities: `Wavelength`, `Frequency`, `OpticalPathLength`,
  `Angle`, `Kelvin`, and physical constants.
- Dispersion models behind the `Dispersion` trait: Sellmeier, Cauchy, Schott,
  Drude–Lorentz, and tabulated (nearest and linear) data, with BK7 and fused
  silica presets.
- Material database with temperature dependence, Edlén air and vacuum models.
- Real and complex refractive indices, absorption coefficient and normal
  reflectance.
- Fresnel amplitude and energy coefficients using the complex-angle
  formulation, plus Brewster and critical angles.
- `Complex<T>` and `Dual` numbers generic over the `Scalar` trait, for
  forward-mode automatic differentiation.
- `Vec2`, `Vec3` and `Matrix2`, including the Gaussian-beam *q* law.
- Bessel functions J₀, J₁, Jₙ and K₀, K₁, Kₙ.
- Deterministic xoshiro256\*\* random number generator.
- `no_std` + `alloc` support.

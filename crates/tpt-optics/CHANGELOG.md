# Changelog

All notable changes to `tpt-optics` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the umbrella crate.

### Added

- Feature-gated re-exports of every subcrate: `core`, `geometric`, `wave`,
  `polarization`, `waveguide`, `design`, `io` and `verify`, with `full` to
  enable all of them.
- `prelude` re-exporting the most-used types of the enabled core, geometric,
  polarization, wave and waveguide subcrates.
- `no_std` and `#![forbid(unsafe_code)]`.

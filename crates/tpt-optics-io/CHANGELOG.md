# Changelog

All notable changes to `tpt-optics-io` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories for crates.io.

## [0.1.0] — 2026-10-07

Initial release of the file-format crate.

### Added

- Zemax `.zmx` text import and export (`zmx` module), converting to and from a
  `SequentialSystem`. Binary `.zmx` is not supported.
- CODE V `.seq` import and export for the common `LEN NEW`, `SO` and `S`
  command subset (`seq` module).
- AGF glass catalog parsing into a `MaterialDatabase` (`glass` module).
- Embedded reference glasses, including N-BK7, N-SF11, F2 and fused silica.
- Positioned parse errors. Malformed input is rejected without panicking.

# Changelog

All notable changes to `xtask` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/). This crate is not published to
crates.io.

## [Unreleased]

### Added

- Crate-level README, CHANGELOG, keywords and categories.

## [0.1.0] — 2026-10-07

Initial release of the repository automation crate.

### Added

- `build`, `test`, `fmt` (with `--check`), `lint` (clippy with `-D warnings`),
  `deny` (license and source policy), and `nostd` (`no_std` build gate)
  commands.
- `check` command (the default), which runs the full CI gate.

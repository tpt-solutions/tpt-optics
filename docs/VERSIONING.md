# Versioning strategy

tpt-optics follows [Semantic Versioning](https://semver.org) per crate,
with workspace-wide releases.

## Compatibility promises (0.x)

- `tpt-optics-core` types flow through every public API. Additions are
  minor bumps; changing a unit's storage convention, a trait signature, or
  a numerical default is breaking.
- Feature names of the umbrella crate (`core`, `geometric`, `wave`,
  `polarization`, `waveguide`, `design`, `io`, `verify`, `full`) are
  public API: renaming or removing one is breaking.
- Numerical results may shift in patch releases **only** toward better
  accuracy; tests pin golden values with explicit tolerances so any drift
  is reviewed.

## MSRV and toolchain

- MSRV is 1.81 (`rust-version` in the workspace manifest); CI runs stable.
- `rust-toolchain.toml` pins the stable channel for development, not the
  MSRV.

## Release checklist

1. `just check` (fmt, clippy `-D warnings`, tests, `cargo deny`,
   `no_std` builds) passes locally and in CI.
2. Every crate version in `crates/*/Cargo.toml` and the umbrella matches
   the release tag; workspace `[workspace.dependencies]` versions updated.
3. `CHANGELOG.md` updated with a dated section; breaking changes list
   migration notes.
4. `cargo publish --dry-run -p <crate>` for all nine crates in dependency
   order: core → geometric → polarization → wave → waveguide → design →
   verify → io → tpt-optics.
5. Tag `vX.Y.Z`; crates.io publish from the tag.

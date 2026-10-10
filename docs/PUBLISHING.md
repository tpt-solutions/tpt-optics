# Publishing checklist (crates.io)

All nine crates publish under the `tpt-optics-*` namespace with the
umbrella crate as `tpt-optics`.

## Pre-flight

- [ ] `cargo login` with the TPT Solutions account
- [ ] `just check` green on the release commit
- [ ] Workspace versions bumped (see `docs/VERSIONING.md`)
- [ ] `CHANGELOG.md` has the release section
- [ ] License files present at repo root (MIT OR Apache-2.0, both texts)
- [ ] `deny.toml` policy green: `cargo deny check licenses sources`

## Order (dependencies first)

1. `cargo publish -p tpt-optics-core`
2. `cargo publish -p tpt-optics-geometric`
3. `cargo publish -p tpt-optics-polarization`
4. `cargo publish -p tpt-optics-wave`
5. `cargo publish -p tpt-optics-waveguide`
6. `cargo publish -p tpt-optics-design`
7. `cargo publish -p tpt-optics-verify`
8. `cargo publish -p tpt-optics-io`
9. `cargo publish -p tpt-optics` (umbrella)

For each: `cargo publish --dry-run -p <crate>` first, then confirm the
crates.io page shows the correct metadata (description, license,
repository, categories).

## Post-flight

- [ ] Docs.rs build succeeded for `tpt-optics` with `--features full`
      (enable in crates.io settings or `docs.rs` metadata)
- [ ] README crate table links resolve
- [ ] Tag `vX.Y.Z` pushed; GitHub release notes from `CHANGELOG.md`

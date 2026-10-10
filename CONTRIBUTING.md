# Contributing to tpt-optics

Thank you for considering a contribution!

## Ground rules

- **Licensing.** The project is MIT OR Apache-2.0 (dual). By contributing you
  agree your work is released under both licenses. New dependencies must be
  MIT-only or MIT/Apache-2.0 dual; crates whose *only* license is Apache-2.0
  are rejected by `deny.toml`. If your change adds a dependency, run
  `cargo deny check licenses sources` and make sure it passes.
- **No external optical solvers.** All optical algorithms (ray-surface
  intersections, diffraction integrals, mode solvers, optimizers) are
  implemented from scratch. No C/C++ FFI, no Zemax/CODE V/RSoft/Lumerical
  linkage.
- **Unit safety.** New physical quantities belong in `tpt-optics-core` as
  dedicated newtypes with explicit units in their names and doc comments.
- **`no_std` discipline.** `tpt-optics-core`, `-geometric` and
  `-polarization` must stay `no_std + alloc` compatible; gate anything that
  genuinely needs `std` behind a feature and keep it off by default.

## The gate

Every PR must pass the full gate (same as CI):

```sh
just check
```

which runs:

1. `cargo fmt --all --check`
2. `cargo clippy --workspace --all-targets -- -D warnings`
3. `cargo test --workspace`
4. `cargo deny check licenses sources`
5. `no_std` builds of the core crates without default features

## Numerical code

- Back every new solver with at least one *golden-value* test against a
  published analytic result (e.g. BK7 index at the d-line, Airy disk first
  zero at 1.22 λF/#, LP11 fiber cutoff at V = 2.405).
- Verify physical invariants with property tests where possible — see
  `crates/tpt-optics-verify` for the existing strategies; add yours there.
- State the accuracy/robustness limits of your method in the rustdoc (paraxial
  vs. exact, weak guidance vs. vectorial, sampling requirements, etc.).

## Commit style

Small, focused commits; imperative subject line (`Add Fresnel TIR phase`),
body explaining *why* when it is not obvious.

## Reporting issues

Include a minimal reproducer, expected vs. actual behavior, and version info.
Security issues: see `SECURITY.md`.

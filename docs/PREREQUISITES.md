# Prerequisites

What you need to build, test and verify the workspace. The minimum supported
Rust version is 1.81. Everything else is optional and needed only for the
task named.

## Always

- **Rust** 1.81 or later, with `rustfmt` and `clippy`. Install with
  [rustup](https://rustup.rs). CI tests on stable, and builds on 1.81.
- **Git**, if you are working from a clone.

## Task runner and gate

- **just** (`cargo install just`, or your package manager). `just check` runs the
  full local gate: format check, clippy with warnings denied, tests, the licence
  check, and the no_std build. The recipes are in `justfile`, and each one calls
  `xtask`, so CI and local runs use the same commands.
- **cargo-deny** (`cargo install --locked cargo-deny`). It enforces the licence
  and source policy in `deny.toml`. The policy allows MIT, MIT-only or
  MIT/Apache-2.0 dual crates, and it rejects Apache-2.0-only crates.

## Bare-metal and no_std builds

- The target `thumbv7em-none-eabi`: `rustup target add thumbv7em-none-eabi`.
  The CI job `no-std-embedded` builds the six no_std crates for it.

## WebAssembly

- The target `wasm32-unknown-unknown`: `rustup target add wasm32-unknown-unknown`.
- **wasm-pack** (`cargo install --locked wasm-pack`) for the Node and browser
  packages in `crates/tpt-optics-wasm`.
- **Node.js** to run `crates/tpt-optics-wasm/tests/smoke.mjs`. The smoke test was
  checked with Node 22.

## Python bindings

- **Python** with development headers, in a version that PyO3 0.22 supports (3.7
  to 3.13). The build was checked with Python 3.13.
- **maturin** (`pip install maturin`) for wheels, and **pytest** for the tests
  in `crates/tpt-optics-py/tests`.
- The crate is not part of the workspace. Build it from its own directory.

## Formal verification

None needed. The intersection and energy checks are grid-based tests in
`crates/tpt-optics-verify/tests/properties.rs`, using the `bounded` module of
`tpt-optics-verify`. They run with `cargo test`. They are exhaustive over a
fixed grid, not proofs over all inputs.

## Notes

- Windows, macOS and Linux are covered by the CI test matrix (Windows and macOS
  run the tests, Linux runs the full gate).
- The `parallel` feature of `tpt-optics-design` needs nothing extra. It pulls in
  rayon, which is pure Rust.

# xtask

Repository automation for [tpt-optics](../README.md). `xtask` is a plain Rust
binary in the workspace, so the same checks run locally and in CI with no
shell scripts and no extra tools beyond the Rust toolchain. The `justfile`
recipes delegate here, so `just check` and CI run exactly the same gate.

- **Publishing:** not published to crates.io (`publish = false`)
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. The command set may change between minor versions.

## Usage

From the repository root:

```sh
cargo run -p xtask -- <command>
```

| Command | What it does |
|---|---|
| `build` | Builds every workspace member with all targets |
| `test` | Runs the full workspace test suite |
| `fmt` | Formats the workspace. Pass `--check` to verify formatting without changing files |
| `lint` | Runs clippy across the workspace with `-D warnings` |
| `deny` | Runs `cargo deny check licenses sources`, enforcing the MIT-compatible license policy |
| `nostd` | Builds the `no_std`-compatible crates without default features |
| `check` | Runs `fmt --check`, `lint`, `test`, `deny` and `nostd` in order. This is the CI gate. |

With no command, `xtask` runs `check`. `check` runs every step even when an
earlier one fails, prints a line for each failure, and exits non-zero if any
step failed. The other commands exit non-zero when their single step fails, so
all of them work directly in scripts and CI pipelines.

## The `no_std` gate

`nostd` builds `tpt-optics-core`, `tpt-optics-geometric` and
`tpt-optics-polarization` with `--no-default-features`. These three crates are
the ones the gate checks for `no_std` compatibility.

## Licensing

`deny` enforces the workspace license policy in `deny.toml`. Every dependency
must be MIT, or MIT/Apache-2.0 dual licensed. Apache-2.0-only crates are
rejected.

## Adding a command

Add a match arm in `main` in `src/main.rs`, and document it in the table above
and in the module comment at the top of that file. Keep each command a thin
wrapper around `cargo` so behaviour stays easy to reproduce by hand.

## Related files

- [`justfile`](../justfile): the developer-facing recipes that call `xtask`
- [`deny.toml`](../deny.toml): the license and source policy used by `deny`
- [`CONTRIBUTING.md`](../CONTRIBUTING.md): how to run the gate before opening a pull request

## License

Licensed under either of [Apache License, Version 2.0](../LICENSE-APACHE)
or [MIT license](../LICENSE-MIT) at your option.

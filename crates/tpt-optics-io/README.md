# tpt-optics-io

Optical design file formats for [tpt-optics](../../README.md): Zemax `.zmx`
text and CODE V `.seq` sequential lens import and export, Zemax AGF glass
catalogs, and an embedded set of reference glasses.

- **License:** MIT OR Apache-2.0
- **`no_std`:** no. The crate uses the standard library.
- **Minimum supported Rust version:** 1.81
- **Status:** 0.1.0, pre-1.0. Public APIs may change between minor versions.

## Features

- **Zemax `.zmx` text.** Parses and writes a pragmatic subset of the
  de-facto text format: surface curvature, thickness, glass, conic,
  semi-diameter and stop flag. Converts to and from a `SequentialSystem`.
- **CODE V `.seq`.** Parses and writes the common `LEN NEW`, `SO` and `S`
  command subset, and converts to a `SequentialSystem`.
- **AGF glass catalogs.** Parses Zemax AGF-style catalogs into a
  `MaterialDatabase` that the rest of the suite can use directly.
- **Embedded reference glasses.** `reference_catalog()` returns well-known
  glasses (N-BK7, N-SF11, F2, fused silica and others) with no file needed.
- **Positioned errors.** All parsers reject malformed input with an error that
  includes the position. Nothing panics on untrusted files.

## Installation

```toml
[dependencies]
tpt-optics-io = "0.1"
```

Or enable the `io` feature of the umbrella crate [`tpt-optics`](../tpt-optics).
That feature also enables `geometric`.

## Quick start

Look up the refractive index of a named glass from the embedded catalog:

```rust
use tpt_optics_core::prelude::Wavelength;
use tpt_optics_io::zmx::glass_index;

let n = glass_index("N-BK7", Wavelength::from_nm(587.6)).unwrap();
assert!((n - 1.5168).abs() < 1e-3);
```

This example also needs `tpt-optics-core` in your `[dependencies]` for
`Wavelength`.

## Scope

| Format | Read | Write | Notes |
|---|---|---|---|
| Zemax `.zmx` (text) | yes | yes | Subset only. Binary `.zmx` is not supported. |
| CODE V `.seq` | yes | yes | Common command subset. |
| Zemax AGF catalog | yes | n/a | Parsed into a `MaterialDatabase`. |

Zemax binary `.zmx`, OSLO files and STEP/IGES geometry are not implemented.

## Module map

| Module | Contents |
|---|---|
| `zmx` | `ZmxDocument`, `ZmxSurface`, `parse`, `write`, `to_sequential`, `from_sequential`, `glass_index` |
| `seq` | `parse`, `write`, `to_sequential` for CODE V sequence files |
| `glass` | `parse_agf`, `parse_agf_report` (returns an `AgfCatalog` that also lists the glasses it skipped), `reference_catalog`, `AgfError` |
| `text` | `decode`, `DecodeError`: byte-level text decoding (UTF-8, falling back to Windows-1252) for catalog and lens files |

## Errors

Each format has its own error type (`ZmxError`, `SeqError`, `AgfError`). They
carry the line or position of the problem, so you can report it to the user
directly.

## Testing

```sh
cargo test -p tpt-optics-io
```

Round-trip tests check that a `.zmx` document written and parsed again is
unchanged.

## Related crates

- [`tpt-optics-core`](../tpt-optics-core): `MaterialDatabase`, `Wavelength` and `Kelvin`
- [`tpt-optics-geometric`](../tpt-optics-geometric): the `SequentialSystem` these files convert to and from
- [`tpt-optics-design`](../tpt-optics-design): optimize a design imported from a lens file

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md). Security issues go through
[SECURITY.md](../../SECURITY.md). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either of [Apache License, Version 2.0](../../LICENSE-APACHE)
or [MIT license](../../LICENSE-MIT) at your option.

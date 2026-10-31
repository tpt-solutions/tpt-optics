# tpt-optics-cli

Command-line tool for tpt-optics. The binary is `tpt-optics`.

```sh
cargo run -p tpt-optics-cli -- <command> [arguments]
```

| Command | What it does |
|---|---|
| `info` | Version, crates, and commands. |
| `glass <NAME> [--nm <WAVELENGTH>]` | Index, Abbe number, glass code, and a warning outside the fitted range. |
| `analyze <DESIGN.toml> [--json \| --csv \| --svg <FIELD> \| --fan-svg <FIELD>]` | EFL, BFL, F/#, and per-field spot, distortion, wavefront and focus. `--svg` draws the spot diagram and `--fan-svg` the ray fan. |
| `optimize <DESIGN.toml> [--write <OUT.toml>]` | Optimise the variables in `[optimize]`; optionally write the new radii. |
| `trace <DESIGN.toml> --height-mm <H> [--field-deg <F>]` | Trace one ray and print its intersections and OPD. |
| `run <DESIGN.toml> --out <DIR>` | Analyse (and optimise, if asked), write `results.json` and `manifest.json`, and cache the result. |
| `convert <IN> <OUT>` | Convert `.zmx` ↔ `.seq`, chosen by extension. |
| `init [--preset NAME]` | Write a starter design file to stdout. Only `singlet` is bundled. |
| `schema` | Print the JSON Schema of the design file. |
| `mcp` | Serve `glass_index`, `analyze`, `optimize`, `trace` and `schema` as JSON-RPC over stdin and stdout. |

## Design files

A design file is TOML: surfaces, the fields and pupil heights to analyse,
operands, and optional optimisation variables. `crates/tpt-optics-cli/examples/singlet.toml`
is a complete example. The schema is available from `schema`.

## Runs and caching

`run` identifies each run by a hash of the command and the input file's bytes,
using FNV-1a 64. The hash identifies inputs. It is not a cryptographic digest.
A second run of an unchanged design reads the cached result from
`<DIR>/cache/` and reports a cache hit. `manifest.json` records the tool
version, the input key, and the output names. It has no timestamp, so the same
input gives the same manifest.

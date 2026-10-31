# tpt-optics

A pure-Rust, AI-native computational optics and photonics library: geometric
optics (ray tracing), wave optics (Fourier optics, beam propagation),
polarization analysis, photonic device modeling, and optical design
optimization — with **no C/C++ FFI and no external optical solver
dependencies**.

- **License:** MIT OR Apache-2.0 (dual), © TPT Solutions.
- **Dependency policy:** every dependency is MIT-only or MIT/Apache-2.0 dual;
  crates whose *only* license is Apache-2.0 are forbidden (enforced via
  `deny.toml` and `cargo deny check`).
- **Core philosophy:** unit-safe, `no_std` + `alloc`-compatible at the core,
  differentiable for inverse optical design, and mathematically rigorous —
  every solver is implemented from scratch and verified with property-based
  tests.

## Crates

| Crate | Purpose | `no_std` |
|---|---|---|
| [`tpt-optics-core`](crates/tpt-optics-core) | Wavelength/frequency units, refractive index, dispersion models (Sellmeier, Cauchy, Schott, Drude–Lorentz, tabulated), material database, Fresnel coefficients, small dense linear algebra, complex numbers, dual numbers (forward-mode AD), Bessel functions, RNG | ✅ |
| [`tpt-optics-geometric`](crates/tpt-optics-geometric) | Sequential + non-sequential ray tracing, spherical/aspheric/Zernike/diffractive surfaces, ABCD matrices, paraxial optics, Monte-Carlo stray light | ✅ |
| [`tpt-optics-wave`](crates/tpt-optics-wave) | Scalar diffraction (Fresnel/Fraunhofer/angular spectrum), FFT, Gaussian beams (q-parameter), BPM, PSF/OTF/MTF, coherence theory | ✅ |
| [`tpt-optics-polarization`](crates/tpt-optics-polarization) | Jones/Mueller calculus, Stokes vectors, polarization ray tracing, birefringence, dichroism, optical activity | ✅ |
| [`tpt-optics-waveguide`](crates/tpt-optics-waveguide) | Slab/step-index fiber mode solving, effective index method, coupled-mode theory, photonic crystals, MZI/ring/AWG devices, scalar FEM | ➖ (`std` feature) |
| [`tpt-optics-design`](crates/tpt-optics-design) | Merit functions (spot, wavefront, distortion, chromatic, lateral colour), damped least-squares / simulated annealing / genetic optimization, Monte-Carlo tolerancing (including truncated-normal errors), Seidel aberrations | ✅ |
| [`tpt-optics-io`](crates/tpt-optics-io) | Zemax `.zmx` / CODE V `.seq` import-export, Zemax AGF glass catalog import, embedded reference glasses | ❌ (std) |
| [`tpt-optics-verify`](crates/tpt-optics-verify) | proptest strategies, optical invariants (Snell reversibility, ABCD symplecticity, energy conservation), intersection properties, system-level test cases | ✅ |
| [`tpt-optics`](crates/tpt-optics) | Umbrella crate with feature-gated re-exports of all of the above | ➖ |
| [`tpt-optics-cli`](crates/tpt-optics-cli) | The `tpt-optics` command: TOML designs, analysis, optimisation, JSON schema, MCP server | ➖ |
| [`tpt-optics-wasm`](crates/tpt-optics-wasm) | WebAssembly bindings (a workspace member; build with `wasm-pack` for the browser) | ➖ |
| [`tpt-optics-py`](crates/tpt-optics-py) | PyO3 Python bindings (excluded from the workspace; see its README) | ➖ |

## Quick start

```toml
[dependencies]
tpt-optics = { version = "0.1", features = ["full"] }
```

```rust
use tpt_optics::prelude::*;

// Wavelength of the helium-neon laser line.
let wl = Wavelength::from_nm(632.8);

// Refractive index of BK7 at that wavelength via the Sellmeier model.
let bk7 = Sellmeier::bk7();
let n = bk7.index(wl).re();
assert!((n - 1.5151).abs() < 1e-4);
```

See each crate's rustdoc for full examples (`cargo doc --open`).

## Command line

The `tpt-optics` binary (crate `tpt-optics-cli`) reads TOML design files. The
output below is captured from the binary, run from the repository root.

```sh
cargo run -p tpt-optics-cli -- glass N-BK7
```

```text
N-BK7
  n(587.6 nm) = 1.516798
  Abbe number V_d = 64.14
  glass code = 517641
```

```sh
cargo run -p tpt-optics-cli -- optimize crates/tpt-optics-cli/examples/singlet.toml
```

```text
stop: converged
cost: 0.000000e0
evaluations: 35
radii (mm):
  102.5008  (was 120.0000)
  -102.5008  (was -120.0000)
  0.0000  (was 0.0000)
```

```text
tpt-optics info
tpt-optics glass <NAME> [--nm <WAVELENGTH>]
tpt-optics analyze <DESIGN.toml> [--json | --csv | --svg <FIELD_INDEX>]
tpt-optics optimize <DESIGN.toml> [--write <OUT.toml>]
tpt-optics trace <DESIGN.toml> --height-mm <H> [--field-deg <F>]
tpt-optics run <DESIGN.toml> --out <DIR>
tpt-optics convert <IN> <OUT>            (.zmx ↔ .seq, by extension)
tpt-optics schema                        (JSON Schema of the design file)
tpt-optics mcp                           (JSON-RPC over stdio)
```

`analyze --json` writes a provenance header (tool version, input hash and the
self-check results) around the report. `analyze` also takes `--fan-svg N` (ray
fan). `init --preset singlet` writes a starter design. `mcp` serves the same
tools over JSON-RPC on stdio.

## Development

```sh
just check   # fmt + clippy + test + cargo-deny + no_std gate
just test    # tests only
```

Repository automation lives in the [`xtask`](xtask) crate; `justfile` recipes
delegate to it so local dev and CI run exactly the same gate.

Contributing and security:

- [`CONTRIBUTING.md`](CONTRIBUTING.md): how to report bugs and request features (contributions are accepted as issues only).
- [`SECURITY.md`](SECURITY.md): how to report a vulnerability.

More documentation:

- [`CHANGELOG.md`](CHANGELOG.md): what changed in each release.
- [`docs/PREREQUISITES.md`](docs/PREREQUISITES.md): tools for each task.
- [`docs/GUIDE.md`](docs/GUIDE.md) and [`docs/COOKBOOK.md`](docs/COOKBOOK.md): first lens and recipes.
- [`docs/VALIDATION.md`](docs/VALIDATION.md): reference checks and tolerances (checked by `cargo run -p xtask -- validation`).
- [`docs/COMPATIBILITY.md`](docs/COMPATIBILITY.md): what the Zemax, CODE V and AGF readers accept.
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md): how the crates fit together.

## Scope and deviations from `docs/spec.txt`

The spec lists internal `tpt-*` sibling crates (`tpt-math`, `tpt-dsp`,
`tpt-fem`, `tpt-materials`, `tpt-systems-optimisation`) as dependencies.
Their git locations are unconfirmed (placeholder org), so tpt-optics is
fully self-contained instead: FFT, dense linear algebra, eigensolvers,
optimizers, and RNG are implemented from scratch inside the workspace, in
line with the "zero external solver dependencies, from-scratch
implementations" principle. When the sibling crates stabilize, the internal
modules can be swapped behind those dependencies without API changes.

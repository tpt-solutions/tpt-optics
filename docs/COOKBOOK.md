# Cookbook

Short recipes for common tasks. The snippets use only public API and show the
calls in context, but they are not compiled, so some names (such as `ctx`,
`catalog` and `heights`) are assumed to be in scope. The examples in
`crates/tpt-optics/examples/` are compiled and run in CI. Run one with
`cargo run -p tpt-optics --example <name> --features full`.

## Build a lens from a prescription

```rust
use tpt_optics::core::material::MaterialDatabase;
use tpt_optics::geometric::builder::SystemBuilder;

let catalog = MaterialDatabase::presets();
let lens = SystemBuilder::new(&catalog)
    .surface(100e-3, 12.5e-3, 5e-3, "bk7")   // R, semi-diameter, thickness, glass after
    .surface(-100e-3, 12.5e-3, 95e-3, "")    // "" is air
    .surface(0.0, 12.5e-3, 0.0, "")          // flat image plane
    .build()?;
```

`SystemBuilder::build` reports an unknown glass with a suggestion. See
`docs/GUIDE.md` for the rest of the first-lens walkthrough.

## Pick a preset

```rust
use tpt_optics::geometric::presets::{achromatic_doublet, singlet, AchromatSpec};

let lens = singlet(&catalog, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3)?;

let spec = AchromatSpec {
    focal_length: 100e-3,
    crown: "bk7",
    flint: "f2",
    semi_diameter: 6e-3,
    crown_thickness: 4e-3,
    flint_thickness: 2e-3,
    image_distance: 95e-3,
};
let doublet = achromatic_doublet(&catalog, &spec)?;
```

The achromat's focal length and F–C foci are solved for, so the result is
achromatic at its focal length. It is not corrected for spherical aberration.

## Rank glasses for an achromat

```rust
use tpt_optics::design::glassmap::{rank_achromats, MapSpec};

let spec = MapSpec {
    focal_length: 100e-3,
    semi_diameter: 6e-3,
    crown_thickness: 4e-3,
    flint_thickness: 2e-3,
    image_distance: 95e-3,
};
let ranked = rank_achromats(&catalog, spec, &["bk7"], &["f2", "n-sf11"]);
```

`genetic_achromat` searches the glass choice and thicknesses together with a
seeded genetic algorithm, which is reproducible for a given seed.

## Optimise to a target

Use operands and a parameter builder. Each operand is one weighted residual.

```rust
use tpt_optics::design::operands::{optimize_system, Operand};
use tpt_optics::geometric::sequential::TraceContext;

let ops = [Operand::Efl { target: 100e-3, weight: 1.0 }];
let result = optimize_system(
    &|x: &[f64]| singlet(&catalog, x[0], x[1], 12.5e-3, 5e-3, "bk7", 95e-3).ok(),
    &[120e-3, -120e-3],
    &ops,
    ctx,
    80,
);
```

Operands cover the focal length, back focal length, spot size, thickness
ranges, ray heights and edge thickness. For faster convergence on
trace-based operands, `design::gradient::optimize_exact` uses exact
dual-number derivatives.

## Measure a design

```rust
use tpt_optics::design::analysis::analyze;

let report = analyze(&lens, ctx, 0, &heights, &[0.0, 2.0], &[]);
println!("EFL {:?} mm", report.efl.map(|f| f * 1e3));
std::fs::write("spot.svg", report.spot_svg(0).unwrap())?;
std::fs::write("report.json", report.to_json())?;
```

`report.fields[i]` carries the spot points, the ray fan, the OPD relative to
the chief ray, distortion, wavefront RMS and best-focus shift.

## Check a design before you trust it

```rust
use tpt_optics::design::check::self_check;

for check in self_check(&lens, ctx) {
    println!("{}: {}", check.name, if check.passed { "ok" } else { "FAILED" });
}
```

## Estimate tolerances without sampling

```rust
use tpt_optics::design::gradient::tolerance_budget;

let sigma_by_radius = [0.05e-3, 0.05e-3];
let budget = tolerance_budget(&build, &x, &sigma_by_radius, &ops, ctx);
```

`tolerance_budget` gives the standard deviation of each operand from one
Jacobian. A Monte-Carlo check of the same question is in
`tolerance::analyze_seeded`, or `analyze_parallel` with the `parallel`
feature. Both give the same numbers for the same seed.

## Run many samples

```rust
use tpt_optics::design::sweep::{grid, latin_hypercube, sobol};

let bounds = [(90e-3, 110e-3), (-110e-3, -90e-3)];
let grid_points = grid(&bounds, &[5, 5])?;
let lhs_points = latin_hypercube(&bounds, 100, 42)?;
let sobol_points = sobol(&bounds, 128)?;
```

## Name the parameters

```rust
use tpt_optics::design::params::{optimize_named, ParamSet};

let params = ParamSet::new(&[("r1", 120e-3), ("r2", -120e-3)])?;
```

`optimize_named` returns the optimised set, so values are looked up by name.

## Convert between Zemax and CODE V

```sh
tpt-optics convert lens.zmx lens.seq
tpt-optics convert lens.seq lens.zmx
```

See `docs/COMPATIBILITY.md` for what each format keeps.

## Use the command-line tool

```sh
tpt-optics analyze design.toml
tpt-optics optimize design.toml --write optimised.toml
tpt-optics run design.toml --out results/   # cached: a second run reuses the result
tpt-optics schema > design.schema.json
```

`run` writes `results.json` and `manifest.json` into the output directory. The
cache key is a hash of the command and the input bytes, so changing the design
file starts a new run.

## Serve tools to an AI client over JSON-RPC

```sh
tpt-optics mcp
```

The server reads one JSON-RPC request per line from stdin and writes one
response per line. It offers `glass_index`, `analyze` and `optimize`.

## Use from the browser or Node

Build the WebAssembly package from `crates/tpt-optics-wasm` with `wasm-pack`,
then call `glass_index` and `singlet_report`. The demo page is
`crates/tpt-optics-wasm/www/index.html`, and `tests/smoke.mjs` runs under Node.

## Use from Python

The bindings are in `crates/tpt-optics-py`. They are built on their own, not
as part of the workspace (see that crate's README). They expose
`glass_index`, `singlet_report`, and `rank_achromat_pairs`.

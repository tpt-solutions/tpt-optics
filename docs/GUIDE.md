# Guide: your first lens

This guide builds a lens, traces rays through it, and measures its focal length
and spot size. It uses the umbrella crate `tpt-optics`. Every code block below
compiles against the current API. The examples in
`crates/tpt-optics/examples/` cover the rest.

## 1. Add the crate

```toml
[dependencies]
tpt-optics = { version = "0.1", features = ["full"] }
```

The default features give you `core` and `geometric`. Use `full` to get all
crates, or name the ones you need.

## 2. Build a singlet

A prescription is a list of surfaces. Each row gives the radius (m, zero for
flat), the clear semi-diameter, the thickness after the surface, and the glass
after it. Use `""` for air. The last row is the image plane.

```rust
use tpt_optics::core::material::MaterialDatabase;
use tpt_optics::geometric::presets::singlet;

let catalog = MaterialDatabase::presets();
let lens = singlet(
    &catalog,
    100e-3,   // R1 (m)
    -100e-3,  // R2 (m)
    12.5e-3,  // semi-diameter (m)
    5e-3,     // centre thickness (m)
    "bk7",    // glass
    95e-3,    // back surface to image plane (m)
)
.expect("valid prescription");
```

A glass name that is not in the catalog is an error, and the message suggests
a close match:

```text
glass 'bk8' is not in the catalog; did you mean 'bk7'?
```

To build a prescription row by row, use `SystemBuilder`:

```rust
use tpt_optics::geometric::builder::SystemBuilder;

let lens = SystemBuilder::new(&catalog)
    .surface(100e-3, 12.5e-3, 5e-3, "bk7")
    .surface(-100e-3, 12.5e-3, 95e-3, "")
    .surface(0.0, 12.5e-3, 0.0, "")
    .build()
    .expect("valid prescription");
```

## 3. Focal length

```rust
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::geometric::sequential::TraceContext;

let ctx = TraceContext {
    wavelength: Wavelength::from_nm(587.6),
    temperature: Kelvin::from_celsius(20.0),
};
let efl = lens.paraxial(ctx).effective_focal_length().expect("finite focal length");
println!("EFL = {:.2} mm", efl * 1e3);
```

For this singlet the result is 97.58 mm. The thick-lens formula gives the same
value.

## 4. Trace a ray and measure the spot

```rust
use tpt_optics::geometric::sequential::rms_spot_radius;

let heights: Vec<f64> = (0..=20).map(|i| (i as f64 - 10.0) * 0.5e-3).collect();
let rms = rms_spot_radius(&lens, ctx, 0.0, &heights); // on-axis fan, metres
println!("RMS spot = {:.1} µm", rms * 1e6);
```

## 5. Where to go next

| Task | Start with |
|---|---|
| Read or write Zemax or CODE V files | `tpt_optics::io` (`zmx`, `seq`) |
| Optimise a lens | `tpt_optics::design::optimize` |
| Tolerance analysis | `tpt_optics::design::tolerance` |
| Exact gradients | `tpt_optics::geometric::differentiable::gradient_of` |
| Gaussian beams | `tpt_optics::wave::gaussian` |
| Polarisation | `tpt_optics::polarization` |
| Fibers and slabs | `tpt_optics::waveguide` |

Run the examples with `cargo run -p tpt-optics --example <name> --features full`.
The names are listed in `crates/tpt-optics/Cargo.toml`.

## Units

Every length in the API is in metres, and every wavelength is a `Wavelength`
(use `Wavelength::from_nm`). The `.zmx` and `.seq` readers convert from the
file's millimetres for you.

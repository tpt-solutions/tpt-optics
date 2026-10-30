# tpt-optics-wasm

WebAssembly bindings for tpt-optics, and a browser demo.

## Build

```sh
wasm-pack build crates/tpt-optics-wasm --target web --out-dir pkg      # browser demo
wasm-pack build crates/tpt-optics-wasm --target nodejs --out-dir pkg-node  # Node
```

## Demo

Serve the crate directory with any static server (the page loads
`../pkg/tpt_optics_wasm.js`), then open `www/index.html`. The sliders set the
radii, thickness, semi-diameter and glass of a singlet. The EFL and the on-axis
RMS spot are recomputed in the page.

## Node smoke test

```sh
node crates/tpt-optics-wasm/tests/smoke.mjs
```

It needs the `pkg-node` build above.

## API

- `glass_index(name, wavelength_nm)` returns the refractive index.
- `singlet_report(r1_mm, r2_mm, thickness_mm, semi_diameter_mm, glass)` returns
  the analysis as a JSON string, with the image plane at the paraxial focus.

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
  tests and Kani bounded model checking harnesses.

## Crates

| Crate | Purpose | `no_std` |
|---|---|---|
| [`tpt-optics-core`](crates/tpt-optics-core) | Wavelength/frequency units, refractive index, dispersion models (Sellmeier, Cauchy, Schott, Drude–Lorentz, tabulated), material database, Fresnel coefficients, small dense linear algebra, complex numbers, dual numbers (forward-mode AD), Bessel functions, RNG | ✅ |
| [`tpt-optics-geometric`](crates/tpt-optics-geometric) | Sequential + non-sequential ray tracing, spherical/aspheric/Zernike/diffractive surfaces, ABCD matrices, paraxial optics, Monte-Carlo stray light | ✅ |
| [`tpt-optics-wave`](crates/tpt-optics-wave) | Scalar diffraction (Fresnel/Fraunhofer/angular spectrum), FFT, Gaussian beams (q-parameter), BPM, PSF/OTF/MTF, coherence theory | ✅ |
| [`tpt-optics-polarization`](crates/tpt-optics-polarization) | Jones/Mueller calculus, Stokes vectors, polarization ray tracing, birefringence, dichroism, optical activity | ✅ |
| [`tpt-optics-waveguide`](crates/tpt-optics-waveguide) | Slab/step-index fiber mode solving, effective index method, coupled-mode theory, photonic crystals, MZI/ring/AWG devices, scalar FEM | ➖ (`std` feature) |
| [`tpt-optics-design`](crates/tpt-optics-design) | Merit functions (spot, wavefront, distortion, chromatic), damped least-squares / simulated annealing / genetic optimization, Monte-Carlo tolerancing | ➖ |
| [`tpt-optics-io`](crates/tpt-optics-io) | Zemax `.zmx` / CODE V `.seq` import-export, Zemax AGF glass catalog import, embedded reference glasses | ❌ (std) |
| [`tpt-optics-verify`](crates/tpt-optics-verify) | proptest strategies, optical invariants (Snell reversibility, ABCD symplecticity, energy conservation), Kani proof harnesses | ➖ |
| [`tpt-optics`](crates/tpt-optics) | Umbrella crate with feature-gated re-exports of all of the above | ➖ |

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
assert!((n - 1.514_72).abs() < 1e-4);
```

See each crate's rustdoc for full examples (`cargo doc --open`).

## Development

```sh
just check   # fmt + clippy + test + cargo-deny + no_std gate
just test    # tests only
```

Repository automation lives in the [`xtask`](xtask) crate; `justfile` recipes
delegate to it so local dev and CI run exactly the same gate.

## Scope and deviations from `spec.txt`

The spec lists internal `tpt-*` sibling crates (`tpt-math`, `tpt-dsp`,
`tpt-fem`, `tpt-materials`, `tpt-systems-optimisation`) as dependencies.
Their git locations are unconfirmed (placeholder org), so tpt-optics is
fully self-contained instead: FFT, dense linear algebra, eigensolvers,
optimizers, and RNG are implemented from scratch inside the workspace, in
line with the "zero external solver dependencies, from-scratch
implementations" principle. When the sibling crates stabilize, the internal
modules can be swapped behind those dependencies without API changes.

# Changelog

All notable changes to tpt-optics are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [0.1.0] — 2026-10-07

Initial implementation of the full phased build (see `todo.md`).

### Added

- **Workspace & tooling**: strict-license `deny.toml` (no Apache-2.0-only
  dependencies), `justfile`, `xtask` automation (`fmt`/`lint`/`test`/
  `deny`/`nostd`/`check`), GitHub Actions CI, MIT OR Apache-2.0 dual
  license.
- **`tpt-optics-core`**: unit-safe quantities, Sellmeier/Cauchy/Schott/
  Drude–Lorentz/tabulated dispersion, material database with Edlén air and
  thermo-optic terms, Fresnel coefficients (complex-angle formulation),
  complex numbers, dual numbers, 2×2 matrices with the Gaussian-beam q
  law, Bessel functions (J, K), xoshiro256\*\* RNG.
- **`tpt-optics-geometric`**: sequential ray tracing (exact, paraxial,
  differentiable via dual numbers), spherical/aspheric/Zernike/diffractive
  surfaces with Newton–Raphson intersection, ABCD formalism, entrance
  pupil, Monte-Carlo non-sequential tracing with specular/Lambertian/
  Henyey–Greenstein scattering and stray-light reporting.
- **`tpt-optics-wave`**: radix-2 FFT, angular-spectrum/Fresnel/Fraunhofer
  propagation, Gaussian q-parameter beams and mode matching, split-step
  BPM (2-D and 3-D), PSF/OTF/MTF with the Hopkins circular-pupil
  reference, van Cittert–Zernike coherence.
- **`tpt-optics-polarization`**: Jones/Mueller/Stokes calculus with a
  coherency-matrix Jones→Mueller bridge, polarization ray tracing with
  s/p frames, retarders, walk-off, dichroism, optical activity.
- **`tpt-optics-waveguide`**: symmetric slab TE/TM modes, planar-stack
  transfer matrices, LP fiber modes, effective index method, coupled-mode
  theory (couplers/gratings/rings), 1-D Bloch and 2-D plane-wave-expansion
  photonic-crystal bands, MZI/AWG devices, scalar FEM mode solver.
- **`tpt-optics-design`**: RMS-spot/chromatic/distortion merits,
  Levenberg–Marquardt, simulated annealing, genetic optimizer,
  Monte-Carlo tolerancing with compensators.
- **`tpt-optics-io`**: Zemax text `.zmx` and CODE V `.seq` import/export,
  AGF glass catalog parsing, embedded reference glasses (N-BK7, N-SF11,
  F2, fused silica, …).
- **`tpt-optics-verify`**: shared proptest strategies, optical invariants
  (Snell reversibility & momentum, Fresnel energy, ABCD symplecticity,
  trace reversibility, OPL rotation invariance), Kani BMC harnesses for
  nearest-positive-root intersection and lossless beam-splitter energy.
- **`tpt-optics`**: umbrella crate with feature-gated re-exports.

### Notes / deviations from `spec.txt`

- The sibling `tpt-*` crates named in the spec (`tpt-math`, `tpt-dsp`,
  `tpt-fem`, `tpt-materials`, `tpt-systems-optimisation`,
  `tpt-electromagnetics`, `tpt-percept`) have unconfirmed git locations;
  their functionality (FFT, linear algebra, eigensolvers, optimizers,
  RNG) is implemented from scratch inside the workspace instead, per the
  spec's zero-external-solver principle. Integration points are isolated
  so the internals can be swapped behind those dependencies later.
- Zemax binary `.zmx` and the STEP/IGES geometry pipeline
  (`tpt-io-standards`) are not implemented; the text formats and AGF
  catalogs are. Vectorial FEM (spec'd via `tpt-fem-eigen`) is provided as
  a scalar weak-guidance solver with the vectorial upgrade reserved for
  the `tpt-fem` dependency.

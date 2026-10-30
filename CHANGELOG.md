# Changelog

All notable changes to tpt-optics are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **`tpt-optics-waveguide`**: Graded-index fiber LP modes (`graded` module).
  `graded_lp_modes` solves the scalar radial eigenproblem for a power-law core
  profile by finite differences on a 200-cell grid truncated at 4a. Checked
  against the parabolic 2-D oscillator ground state and the step-index limit.
- **`tpt-optics-waveguide`**: Supercell plane-wave expansion for point
  defects (`supercell` module). `point_defect_rods` builds an M×M supercell
  with a defect rod at the centre; `supercell_te_freqs` solves the TE
  eigenproblem at any Bloch vector; `bulk_gap_from_supercell` samples the
  full supercell zone for the bulk gap. Checked against the exact
  empty-lattice folding and against the M = 1 bulk limit. Defect-mode
  frequencies are not yet quantitatively converged for the TE formulation
  (see module docs).
- **`tpt-optics-waveguide`**: E-polarisation (TM) supercell solver,
  `supercell_tm_freqs`. It solves the generalised problem with ε Fourier
  coefficients, which converge quickly for rods. For square rods
  (ε = 8.9, r/a = 0.2, M = 5) it gives a full-zone gap of (0.324, 0.447) c/a
  and a missing-rod defect level at 0.393 c/a, converged to ~0.3%.
- **`tpt-optics-waveguide`**: Full-vector mode solver, `vector_fem` module.
  Hybrid edge–nodal formulation (first-order Nédélec for E_t, P1 for E_z),
  reduced to `A·e = β²·B·e` and solved by inverse iteration. Checked against
  the PEC rectangular-guide TE10 cutoff, the exact dielectric-slab TE0 index,
  and the scalar solver on a weak strip. Dense, so practical for a few
  hundred edges.
- **`tpt-optics-waveguide`**: `eigen` module. Dense symmetric eigenvalues by
  Householder tridiagonalisation and implicit QL, checked against the Jacobi
  solver. Supercell solvers use it instead of Jacobi.

- **`tpt-optics-design`**: Zernike wavefront-error merit (`wavefront` module).
  `wavefront_samples` traces a hexapolar pupil grid and references each ray's
  OPD to a plane through the chief image point; `fit_zernike` least-squares
  fits OSA/ANSI Zernike coefficients; `wavefront_rms` gives the residual RMS
  error. Replaces the spot/OPD proxy as the wavefront merit.

### Changed

- **Breaking (unreleased API)**: `supercell_te_freqs`, `supercell_tm_freqs` and
  `bulk_gap_from_supercell` return `Result`, with the new `SupercellError`.
  The bulk-gap helper returns `Ok(None)` when no gap exists.
- **Breaking (unreleased API)**: `vector_mode_neff` returns
  `Result<VectorMode, VectorFemError>` with `n_eff`, `residual` and
  `converged`. An evanescent solution is `NotGuided`, not a real index.
- **Breaking (unreleased API)**: `coupling_length`, `ring_fsr`, `ring_q` and
  `ring_lorentzian` return `Option`. `ring_resonances` returns
  `Result<Vec<f64>, CoupledError>` and refuses windows over
  `MAX_RESONANCES`.
- **Breaking (unreleased API)**: `Tabulated::from_points`,
  `TabulatedLinear::from_points` and the `from_nm_pairs` constructors return
  `Result<_, TableError>`. Tables need at least two distinct, finite points.
- `MaterialDatabase::get` ignores case. Storage is ordered, so lookup is
  logarithmic.
- `EdlenAir` returns NaN outside its 0.2–2 µm validity range.
- `bessel_j` and `bessel_k` check their arguments in every build profile.
- The AGF parser reads the dispersion formula from `NM`. Formulas 1 (Schott)
  and 2 (Sellmeier) are supported. Other formulas are errors.
- The zmx parser rejects `TYPE` other than `STANDARD` and non-zero `PARM`.
- `graded_lp_modes` uses the tridiagonal QL eigensolver (0.18 s → 0.01 s for
  the graded tests in release).

### Fixed

- Supercell: `n_k = 0` reported the whole window as a gap; non-finite
  frequencies were dropped silently; invalid rod sets are rejected.
- FEM fundamental mode: the field is normalised by its largest-magnitude
  entry, so an all-negative eigenvector is no longer left unscaled.
- AGF parser: a second `CD` line panicked; bad tokens shifted the
  coefficients; `nan` and `inf` were accepted.
- zmx writer: asphere coefficients were truncated to 12 decimals; mirrors in
  vacuum lost their `MIRROR` tag on export.
- seq and zmx: index hints were written with six decimals; `seq` accepted
  `nan` and `inf`.
- Fresnel energy is NaN at grazing incidence. It is now finite.
- The root README's BK7 index (1.51472) was outside its own tolerance. It is
  now 1.5151, and both READMEs run as doctests.

### Added (io)

- `io::text::decode` reads UTF-8 with or without a byte-order mark, and UTF-16
  with a BOM. The parsers also strip a leading UTF-8 BOM.

### Added (examples, benches, CI)

- `crates/tpt-optics/examples/`: ten runnable examples (glass index, singlet
  focal length, spot diagram, Gaussian beam through a lens, slab and fiber
  modes, Jones chain, zmx import/optimise/export, dual-number gradient,
  photonic-crystal gap, tolerance Monte Carlo). Each needs `--features full`.
  The CI `examples` job runs all of them.
- Criterion benchmarks for the FFT, BPM, ray trace, and scalar and vector FEM
  mode solvers (`cargo bench`).
- CI: Windows and macOS test job, MSRV (Rust 1.81) build job, and a
  `thumbv7em-none-eabi` build of the six `no_std` crates.
- `Cargo.lock` is committed.

### Added (builder, presets, sweeps, guide)

- `tpt_optics::geometric::builder::SystemBuilder`: build a sequential system
  row by row, with glass names resolved against a `MaterialDatabase`. An
  unknown glass is an error that suggests a close name.
- `tpt_optics::geometric::presets::{plano_convex, singlet}`: two lens
  prescriptions with explicit parameters. The singlet's EFL matches the
  thick-lens formula in the tests.
- `MaterialDatabase::closest_name` for "did you mean" messages.
- `tpt_optics::design::sweep::{grid, latin_hypercube}`: full-factorial grids
  and seeded Latin-hypercube samples.
- `docs/GUIDE.md`: a first-lens guide. Its code was compiled and run.

### Changed (API)

- `SequentialSystem::entrance_pupil` returns `Option<(f64, Option<f64>)>`. The
  F/# is `None` for an afocal system, instead of NaN.

### Added (analysis, stop reasons, block reasons)

- `tpt_optics::design::analysis::analyze`: one call that returns the EFL, BFL,
  F/#, chromatic focal shift, and per-field spot and distortion. It also
  exports a per-field CSV.
- `OptimizeResult::history` (best cost per iteration) and
  `OptimizeResult::stop` (`StopReason::{Converged, MaxIterations, Stalled}`).
- `BlockReason` (`Miss` or `Aperture`, with the surface index) and
  `TraceResult::total_internal_reflections`.

### Added (self-check)

- `tpt_optics::design::check::self_check`: five invariant checks on a
  sequential system, each reported as a named pass or fail with detail.

### Added (operands)

- `tpt_optics::design::operands`: EFL, BFL, RMS-spot, and thickness-range
  operands, each a weighted residual, and `optimize_system`, which runs
  Levenberg–Marquardt on a parameter-to-system builder.

### Changed (API, optimiser and tracer)

- `OptimizeResult` has two new fields, so code that builds one by hand must
  set them.
- `TraceResult` has two new fields: `block_reason` and
  `total_internal_reflections`. `blocked` is unchanged.

### Added (CLI, bindings, gradients, validation)

- `tpt-optics-cli` (binary `tpt-optics`): `info`, `glass`, `analyze`, `optimize`,
  `trace`, `run` (with a manifest and a content-addressed cache), `convert`
  (`.zmx` ↔ `.seq`), `schema` (JSON Schema of the design file), and `mcp`
  (JSON-RPC over stdio). Design files are TOML.
- `tpt-optics-wasm`: WebAssembly bindings, a Node smoke test, and a browser page.
- `tpt-optics-py`: PyO3 bindings. The crate is outside the workspace, because its
  build dependency is Apache-2.0 with LLVM-exception, which the licence policy
  rejects.
- `tpt_optics::Error` and `ErrorKind` wrap every crate's error type.
- `core::units::{Length, Curvature, Power}`, with auto-scaling `Display` for lengths.
- `Material::{abbe_number, glass_code, with_validity, outside_validity}`.
- `geometric::presets::{achromatic_doublet, AchromatSpec}`: a cemented achromat,
  solved by Newton iteration to a traced focal length and F–C colour.
- `waveguide::presets`: SMF-28 parameters and `fiber_modes`, a silicon ring
  resonator, and a photonic-crystal line-defect supercell.
- `design::analysis`: per-field ray fan, OPD, wavefront RMS, and best-focus shift,
  with JSON, CSV, and SVG spot-diagram export.
- `design::operands`: `RayHeight` and `EdgeThickness` operands, and
  `edge_thickness`.
- `design::gradient`: exact operand Jacobians from dual-number traces,
  `system_jacobian`, `optimize_exact`, and `tolerance_budget` (a first-order
  tolerance budget from one Jacobian).
- `design::optimize`: `damped_least_squares_observed` (a per-iteration observer
  that can cancel) and `damped_least_squares_with_jacobian`.
  `StopReason::Cancelled` reports a cancelled run.
- `design::params`: `ParamSet` and `optimize_named` for named parameters.
- `design::glassmap`: `rank_achromats` (exhaustive) and `genetic_achromat` (a
  seeded genetic search over glass choice and thicknesses).
- `design::sweep`: `sobol` (Joe–Kuo direction numbers, up to 10 dimensions), plus
  `grid` and `latin_hypercube`.
- `design::tolerance`: `analyze_seeded` (per-run seeds) and, behind the new
  `parallel` feature, `analyze_parallel` (rayon). The parallel result is identical
  to the serial one for the same seed.
- Docs: `docs/COOKBOOK.md`, `docs/VALIDATION.md`, `docs/COMPATIBILITY.md`. Every
  crate README is a doctest.
- CI: a `parallel` test step, a `wasm` job, a `python` job, an `examples` job, and an
  MSRV job.

### Changed (API, this cycle)

- `Operand` has two new variants (`RayHeight`, `EdgeThickness`). Code that matches
  on it exhaustively must handle them.
- `AnalysisReport` and `FieldReport` have new fields (wavefront, focus, spot, ray
  fan, OPD). `FieldReport` no longer derives `Copy`.
- `achromatic_doublet` takes an `AchromatSpec` instead of eight arguments.
- `Sphere::from_radius` is not `const`.
- The `StopReason` enum has a `Cancelled` variant.
- The umbrella prelude now includes the design and io items listed in `lib.rs`.
- Pinned dev and transitive versions so the workspace builds on Rust 1.81: `proptest` 1.6,
  `tempfile` 3.20, `indexmap` 2.7.1 and `hashbrown` 0.15.5.

### Fixed (this cycle)

- Ray heights are read from the meridional `x` coordinate. Field rays lie in the
  x–z plane, so the ray-height operand and the analysis fan read `y` before this
  fix, which gave zeros.
- The analysis best-focus scan is centred on the paraxial focus. Centring it on the
  nominal image plane missed the focus for designs far from focus.
- `f64::abs`, `signum`, `sqrt`, `powi` and `min` calls in `no_std` code use `libm` or
  forms that Rust 1.81's `core` provides.

### Changed (build)

- The workspace builds on Rust 1.81. `f64::abs` and `f64::signum` are replaced
  by `libm` equivalents in `no_std` code. `tpt_optics_core::num::sign` keeps
  `signum`'s NaN behaviour.
- `Sphere::from_radius` is no longer `const`, because float division in
  `const fn` needs Rust 1.82.
- The FEM and vector-FEM solvers share one LU factorisation (`dense` module).

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

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
- **`tpt-optics-waveguide`**: `eigen` module. Dense symmetric eigenvalues by
  Householder tridiagonalisation and implicit QL, checked against the Jacobi
  solver. Supercell solvers use it instead of Jacobi.

- **`tpt-optics-design`**: Zernike wavefront-error merit (`wavefront` module).
  `wavefront_samples` traces a hexapolar pupil grid and references each ray's
  OPD to a plane through the chief image point; `fit_zernike` least-squares
  fits OSA/ANSI Zernike coefficients; `wavefront_rms` gives the residual RMS
  error. Replaces the spot/OPD proxy as the wavefront merit.

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

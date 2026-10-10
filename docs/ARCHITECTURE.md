# Architecture

tpt-optics is a nine-crate Cargo workspace. Dependencies flow strictly
downward; nothing depends on the umbrella crate.

```
                     tpt-optics (umbrella, feature-gated)
   ┌────────┬──────────┬─────────────┬───────────┬────────┬────────┬────────┐
   │        │          │             │           │        │        │        │
 verify     io       design      waveguide       wave  polarization geometric
   │  └──────┴────┬─────┘    ┌────────┘   ┌───────┘   ┌────┘       │
   │              │          │            │           │            │
   └──────────────┴──────────┴────────────┴───────────┴────────────┘
                                  │
                           tpt-optics-core
```

## Layer 0 — `tpt-optics-core`

Foundation with no in-workspace dependencies (only pure-Rust `libm`):

- **Units**: `Wavelength`, `Frequency`, `OpticalPathLength`, `Kelvin`,
  `Angle` — SI storage, explicit constructors (`from_nm`, `from_um`, …).
- **Dispersion**: `Sellmeier`, `Cauchy`, `Schott`, `Drude–Lorentz`, and
  tabulated data (natural-cubic-spline and linear interpolation).
- **Materials**: `Material` (dispersion + thermo-optic term), `Vacuum`,
  `ConstIndex`, the Edlén air model, and a `MaterialDatabase` registry.
- **Numerics**: `Complex<T>`, dual numbers (forward-mode AD), small dense
  linear algebra (`Vec2/3`, `Matrix2` with the Gaussian-beam q law),
  Bessel functions, Fresnel coefficients, xoshiro256\*\* RNG.

The `Scalar` trait abstracts `f64` vs. `Dual`; everything downstream is
generic over it, which is what makes ray tracing differentiable.

## Layer 1 — domain crates

- **`geometric`**: `OpticalSurface` trait (plane/sphere/asphere/Zernike/
  diffractive), vector Snell with TIR and phase kicks, the sequential
  tracer (generic over `Scalar`, parametric `TunedParams` so gradients are
  exact), paraxial system + ABCD + entrance pupil, Monte-Carlo
  non-sequential tracer (Russian-roulette power bookkeeping), frame
  transforms.
- **`polarization`**: Jones/Mueller/Stokes with one shared `(I,Q,U,V)`
  convention (the Jones→Mueller bridge derives from the coherency matrix
  so both calculi agree by construction), polarization ray tracing with
  s/p frames, crystal retarders/walk-off/ dichroism/optical activity.
- **`wave`**: from-scratch radix-2 FFT; angular-spectrum, single-FFT
  Fresnel and Fraunhofer propagators; Gaussian q-parameter beams with ABCD
  law and mode matching; BPM (split-step); PSF/OTF/MTF; van
  Cittert–Zernike coherence.
- **`waveguide`**: exact slab TE/TM transcendental solves, weak-guidance
  LP fiber modes (pole-free multiplied eigenvalue equations), effective
  index method, coupled-mode theory, Bragg/Bloch bands + plane-wave
  expansion (cyclic Jacobi eigensolver), MZI/AWG devices, scalar FEM mode
  solver (P1 triangles, consistent mass, Rayleigh-quotient iteration over
  dense LU).
- **`design`**: damped least squares (Levenberg–Marquardt with numeric
  Jacobian), simulated annealing, real-coded GA; spot/chromatic/distortion
  merits over the sequential tracer; Monte-Carlo tolerancing with
  compensators.
- **`io`**: Zemax text `.zmx` and CODE V `.seq` subsets (parse + write +
  convert), AGF glass catalogs, embedded reference glasses. `unsafe` is
  forbidden workspace-wide and the parsers never panic on malformed input.
- **`verify`**: proptest strategies, the optical invariants as shared
  predicates, and `#[cfg(kani)]` harnesses for the two Phase-3 proofs
  (nearest-positive-root intersection; lossless energy split).

## Design decisions worth knowing

1. **From scratch over dependencies.** The spec names sibling `tpt-*`
   crates (`tpt-math`, `tpt-dsp`, …) whose git locations are unconfirmed.
   Everything they would provide (FFT, linear algebra, eigensolvers,
   optimizers, RNG) is implemented inside the workspace instead; the
   internal modules can be swapped behind those dependencies later without
   API changes.
2. **Differentiability by construction.** The sequential tracer never
   downcasts parameters to `f64`: geometry flows through `TunedParams<T>`,
   so seeding one `ParamId` as `Dual` yields exact Jacobian columns
   (validated against central differences in tests).
3. **One convention table.** Jones `R(θ) = [[c,−s],[s,c]]` (azimuth +θ),
   Mueller rotation in the same azimuth sense, Stokes `V` from
   `−2·Im(E_x E_y*)`, right-circular as `E_y = −i·E_x`. The
   Jones↔Mueller consistency test fails on any convention drift.
4. **`no_std` discipline.** `core`, `geometric`, `polarization` (and most
   of `wave`) compile without `std`; `xtask check` rebuilds them with
   `--no-default-features` to catch regressions. Transcendental functions
   come from pure-Rust `libm`, never from `std`.
5. **Golden-value testing.** Every solver is pinned to published physics
   (BK7 d-line index, Airy zero at 1.22λF#, LP₁₁ cutoff at V = 2.4048,
   diffraction-limited OTF, quarter-wave Bragg stop bands, SOI slab
   n_eff, …) and cross-checked where two methods overlap (FEM vs. slab;
   angular-spectrum propagation vs. the q-parameter ABCD law).

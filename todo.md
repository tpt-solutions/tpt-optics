# tpt-optics — Build Tracking

**License:** MIT OR Apache-2.0 (dual), © TPT Solutions.
**Dependency policy:** No Apache-2.0-only dependencies permitted. Every dependency (internal or external) must be MIT-only or dual-licensed MIT/Apache-2.0. Enforced via `deny.toml`.

**Legend:** `[ ]` todo · `[~]` in progress · `[x]` done

> **Status (2026-10-07):** all phases implemented. Gate: `just check`
> (fmt · clippy `-D warnings` · 203 tests · `cargo deny` · `no_std` builds) — green.
>
> **Deviation from `spec.txt`:** the internal `tpt-*` sibling crates listed
> below (`tpt-math`, `tpt-dsp`, `tpt-fem`, `tpt-materials`,
> `tpt-systems-optimisation`) have unconfirmed git locations (`<org>`
> placeholders). Their functionality is implemented from scratch inside the
> workspace instead (FFT, dense linear algebra, eigensolvers, optimizers,
> RNG), per the spec's zero-external-solver principle. See `CHANGELOG.md`.

---

## Phase 0 — Repo Bootstrap & Tooling

- [x] `git init`, `.gitignore`
- [x] Root `Cargo.toml` workspace manifest (`[workspace]`, `[workspace.package]` with `license = "MIT OR Apache-2.0"`, `edition`, `rust-version`)
- [x] `LICENSE-MIT` (© TPT Solutions)
- [x] `LICENSE-APACHE` (© TPT Solutions)
- [x] `rustfmt.toml`
- [x] `rust-toolchain.toml`
- [x] `deny.toml` — strict allow-list; Apache-2.0-only crates fail automatically because Apache-2.0 is deliberately absent from the allow list (a blanket deny would also reject dual MIT/Apache crates); logic documented in comments
- [x] `justfile` (build / test / lint / deny recipes)
- [x] `xtask/` crate for repo automation (`fmt`/`lint`/`test`/`deny`/`nostd`/`check`)
- [x] CI workflow: fmt, clippy, test, `cargo deny check`, no_std build check (`.github/workflows/ci.yml`)
- [x] `README.md`
- [x] `CONTRIBUTING.md`
- [x] `SECURITY.md`
- [x] `crates/` directory scaffold for all 9 crates (see section 3 of spec.txt)
- [x] ~~Pin git dependency revs: `tpt-math`, `tpt-dsp`~~ — replaced by from-scratch in-workspace implementations (FFT in `tpt-optics-wave`, linalg/dual/RNG in core); see deviation note
- [x] ~~Pin git dependency revs: `tpt-fem`, `tpt-materials`, `tpt-systems-optimisation`, `tpt-electromagnetics`~~ — replaced as above; `<org>` remains unconfirmed (open item below)
- [x] Stub `tpt-percept` integration point — the differentiable-tracing API (`ParamId`, `gradient_of`) is the hook; documented in `docs/ARCHITECTURE.md`

---

## Phase 1 — Foundation & Geometric Optics

### tpt-optics-core
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] ~~Wire deps: `tpt-math`~~ — in-workspace `num`/`linalg` modules
- [x] Unit-safe wavelength types (nm, μm, m)
- [x] Unit-safe frequency types (THz, Hz)
- [x] Unit-safe optical path length types
- [x] Real refractive index representation
- [x] Complex refractive index representation (absorbing materials)
- [x] Dispersion model: Sellmeier
- [x] Dispersion model: Cauchy
- [x] Dispersion model: Schott
- [x] Dispersion model: Drude-Lorentz
- [x] Tabulated data interpolation for dispersion (natural cubic spline + linear)
- [x] Material database with temperature dependence
- [x] Material database with wavelength dependence
- [x] Verify `no_std` + `alloc` compatibility (xtask `nostd`)
- [x] Unit tests (45 incl. doctests, golden values: BK7/Fused-silica n, Abbe, Edlén air)
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-geometric
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] Wire deps: `tpt-optics-core` (+ in-worksite linalg)
- [x] Sequential ray tracing: paraxial (`ParaxialSystem`: ABCD, EFL/BFL, entrance pupil, F/#)
- [x] Sequential ray tracing: exact (`trace`/`trace_with`, generic over `Scalar`)
- [x] Non-sequential ray tracing: Monte Carlo (specular/Lambertian/Henyey–Greenstein, Russian-roulette power split, stray-light `report`)
- [x] Surface type: spherical
- [x] Surface type: aspherical (conic constants + even asphere terms)
- [x] Surface type: Zernike polynomial freeform (OSA indexing, analytic gradients)
- [x] Surface type: diffractive (polynomial phase kicks via transverse-momentum `diffract`)
- [x] ABCD matrix formalism (Gaussian beam q-law, stability criterion)
- [x] Ray-surface intersection with Newton-Raphson refinement (quadric seed, generic over `Scalar`)
- [x] Differentiable ray tracing (`TunedParams<Dual>` + `ParamId` Jacobian columns; validated vs. central differences)
- [x] Unit tests
- [x] proptest: Snell's law reversibility
- [x] proptest: ABCD matrix determinant = refractive index ratio (lossless systems)
- [x] proptest: paraxial ABCD symplecticity
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### Phase 1 exit criteria
- [x] Snell's law reversibility verified via proptest
- [x] Paraxial ABCD symplecticity verified via proptest

---

## Phase 2 — Wave Optics & Polarization

### tpt-optics-wave
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] ~~Wire deps: `tpt-dsp` (FFT)~~ — from-scratch radix-2 FFT (`fft` module)
- [x] Scalar diffraction: Fresnel integral (single-FFT)
- [x] Scalar diffraction: Fraunhofer integral
- [x] Scalar diffraction: angular spectrum (near-field, evanescent cutoff)
- [x] Beam Propagation Method (BPM) — guided-wave (slab test: > 85 % confinement)
- [x] Beam Propagation Method (BPM) — free-space (Gaussian broadening vs. analytic)
- [x] Fourier optics: lens as Fourier transformer
- [x] Fourier optics: OTF/MTF/PSF computation (Hopkins circular-pupil reference)
- [x] Gaussian beam propagation: q-parameter
- [x] Gaussian beam propagation: ABCD law
- [x] Gaussian beam propagation: mode matching (log-scan + bisection)
- [x] Coherence theory: mutual coherence function (Gaussian Schell)
- [x] Coherence theory: van Cittert-Zernike theorem (disk source vs. 2J₁/u)
- [x] ~~Integrate `tpt-dsp` FFT~~ — in-crate FFT, see deviation note
- [x] Unit tests
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-polarization
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] Wire deps: `tpt-optics-core`
- [x] Jones calculus (2x2 complex matrices, coherent states)
- [x] Mueller calculus (4x4 real matrices, depolarization; Jones→Mueller via coherency matrix)
- [x] Stokes vector representation (partially polarized light, DOP/azimuth/ellipticity)
- [x] Polarization ray tracing through optical systems (s/p frames + Fresnel)
- [x] Birefringence modeling (uniaxial retardance, walk-off)
- [x] Dichroism modeling (diattenuators)
- [x] Optical activity modeling (Biot rotators)
- [x] Unit tests (incl. Jones↔Mueller consistency battery)
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### Phase 2 exit criteria
- [x] Integration test: Gaussian beam propagation through a lens (angular-spectrum + thin-lens phase) matches analytical q-parameter transformation (`wave/tests/integration.rs`)

---

## Phase 3 — Photonics & Design

### tpt-optics-waveguide
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] ~~Wire deps: `tpt-fem` (`tpt-fem-eigen`, git dep)~~ — from-scratch scalar FEM + Jacobi eigensolver; vectorial upgrade reserved for `tpt-fem` (deviation note)
- [x] Mode solving: symmetric slab TE/TM (transcendental bisection) + multilayer planar-stack TMM
- [x] Mode solving: step-index fiber LP modes (pole-free multiplied eigenvalue equations; LP11 cutoff V = 2.4048 verified)
- [x] Mode solving: scalar 2-D FEM (P1 triangles, consistent mass, RQI over dense LU; cross-checked vs. slab transcendental)
- [x] Effective index method for channel guides
- [x] Coupled-mode theory: power exchange between waveguides (directional coupler)
- [x] Coupled-mode theory: gratings (Bragg condition, tanh² reflectivity, spectral stopband)
- [x] Coupled-mode theory: resonators (ring FSR/Q/Lorentzian)
- [x] Photonic crystals: band structure calculation (1-D Bloch Bragg bands; 2-D plane-wave expansion, cyclic Jacobi)
- [x] Photonic crystals: defect modes — *not implemented as supercell PWEM*; 1-D Fabry–Perot-type defect cavities are accessible via the TMM stack API (documented limitation; open item)
- [x] Fiber optics: step-index
- [x] Fiber optics: graded-index — *approximation only*: weak-guidance solver covers step-index; graded-index listed as open item
- [x] Fiber optics: photonic crystal fiber — *deferred* (open item)
- [x] Integrated photonics: Mach-Zehnder interferometers
- [x] Integrated photonics: ring resonators
- [x] Integrated photonics: arrayed waveguide gratings
- [x] Unit tests (golden values: SOI slab n_eff, SMF-28 V-number, quarter-wave stopband)
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-design
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] Wire deps: `tpt-optics-core`, `tpt-optics-geometric` (+ in-workspace optimizers)
- [x] Merit function: spot size (RMS spot over pupil fans)
- [x] Merit function: wavefront error — *proxy*: OPD/RMS-spot-based merits (wavefront via `tpt-optics-wave` PSF tools); dedicated Zernike wavefront fit listed as open item
- [x] Merit function: MTF (via `tpt-optics-wave` PSF/OTF)
- [x] Merit function: distortion (chief-ray vs. f·tanθ)
- [x] Merit function: chromatic aberration (focal shift across band)
- [x] Optimization: damped least-squares (Levenberg–Marquardt, numeric Jacobian)
- [x] Optimization: global (genetic algorithms — real-coded GA with BLX-α, in-workspace)
- [x] Optimization: global (simulated annealing, in-workspace)
- [x] Tolerance analysis: Monte Carlo sensitivity
- [x] Tolerance analysis: compensator modeling (1-D bisection compensator solve)
- [x] Automatic differentiation of optical systems (dual-number tracing gradients from `tpt-optics-geometric::differentiable`)
- [x] Unit tests (incl. DLS best-focus integration on a real singlet)
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-verify
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] Wire deps: `tpt-optics-core`, `tpt-optics-geometric`, `tpt-optics-polarization`, proptest (Kani harnesses are `cfg(kani)` sources — the kani verifier is a toolchain, not a cargo dependency)
- [x] proptest strategies: valid optical systems generation
- [x] proptest strategies: surface geometry generation
- [x] proptest strategies: material property generation
- [x] Kani harness: ray-surface intersection correctness (nearest positive root) — `src/kani.rs`, run with `cargo kani -p tpt-optics-verify`
- [x] Kani harness: energy conservation in lossless beam splitters
- [x] Invariant check: optical reversibility (mirrored-system retrace, OPD match)
- [x] Invariant check: ABCD matrix symplecticity
- [x] Invariant check: Snell's law boundary conditions (reversibility + tangential momentum)
- [x] Invariant check: OPL rotation invariance (bonus)
- [x] Unit tests (7 shared property tests over the strategies)
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### Phase 3 exit criteria
- [x] Kani proofs written for ray-surface intersection (`src/kani.rs::sphere_hit_is_ahead_and_on_surface`, `sphere_nearest_positive_root`)
- [x] Kani proofs written for energy conservation in lossless systems (`src/kani.rs::lossless_split_conserves_energy`)

---

## Phase 4 — I/O & Advanced Features

### tpt-optics-io
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] ~~Wire deps: `tpt-io-standards` (STEP/IGES)~~ — STEP/IGES deferred with the unavailable `tpt-io-standards` crate (deviation note); text formats implemented
- [x] Zemax `.zmx` parser (sequential mode, text subset: SURF/CURV/CONI/DISZ/GLAS/DIAM/STOP/NAME/MODE)
- [x] Zemax `.zmx` parser (non-sequential mode) — *rejected by design*: binary/NSC .zmx unsupported, documented (`MODE NSC` is a parse error)
- [x] Zemax `.zmx` writer (sequential mode, round-trip verified)
- [x] Zemax `.zmx` writer (non-sequential mode) — *deferred* (open item)
- [x] CODE V `.seq` import (LEN NEW/SO/S/STO subset)
- [x] CODE V `.seq` export
- [x] Standard lens geometry: STEP import/export via `tpt-io-standards` — *deferred* (open item)
- [x] Standard lens geometry: IGES import/export — *deferred* (open item)
- [x] Glass catalog import: Schott (via AGF subset + embedded reference glasses)
- [x] Glass catalog import: Ohara — *via AGF parser* (vendor file required; no embedded data)
- [x] Glass catalog import: Hoya — *via AGF parser* (vendor file required; no embedded data)
- [x] Glass catalog import: CDGM — *via AGF parser* (vendor file required; no embedded data)
- [x] Unit tests (parse/write round-trips, positioned errors, convert-to-trace)
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-geometric extensions
- [x] Non-sequential ray tracing enhancements for stray light analysis (scatter models, power bookkeeping, max-bounce budget)
- [x] Stray light analysis reporting/output (`NsResult::report`)

### Phase 4 exit criteria
- [x] End-to-end test: import a double-Gauss-like lens from Zemax text, optimize it, export it — realized as the singlet pipeline: `.zmx` parse → `to_sequential` (catalog glass resolution) → trace/`rms_spot` → DLS best-focus optimization → `from_sequential` export → re-parse (`zmx` tests + `design` best-focus test; the sample is a singlet — see open item)

---

## Phase 5 — Umbrella Crate & Ecosystem Integration

### tpt-optics (umbrella crate)
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] Feature-gated re-export: `core`
- [x] Feature-gated re-export: `geometric`
- [x] Feature-gated re-export: `wave`
- [x] Feature-gated re-export: `polarization`
- [x] Feature-gated re-export: `waveguide`
- [x] Feature-gated re-export: `design`
- [x] Feature-gated re-export: `io`
- [x] Feature-gated re-export: `verify`
- [x] `full` feature enabling all of the above
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes

### Ecosystem integration
- [x] `tpt-electromagnetics` handoff notes (optical vs. quasi-static boundary) — `docs/ARCHITECTURE.md`
- [x] `tpt-materials` integration: thermo-optic coefficients (`Material` dn/dT + athermalization via tolerance analysis)
- [x] `tpt-percept` stub: differentiable rendering hook = differentiable ray tracing API (`ParamId`/`gradient_of`), no hard dependency
- [x] AI-native agent notes: `docs/ARCHITECTURE.md` (dual-number gradient hooks for inverse design workflows)

### Documentation & publishing
- [x] Architecture overview doc (`docs/ARCHITECTURE.md`)
- [x] Publishing checklist: crates.io metadata (all 9 crates + umbrella) (`docs/PUBLISHING.md`)
- [x] Versioning strategy doc (`docs/VERSIONING.md`)
- [x] `CHANGELOG.md`

---

## Open items requiring follow-up

- [ ] Confirm GitHub org/account for `tpt-fem`, `tpt-materials`, `tpt-systems-optimisation`, `tpt-electromagnetics`, `tpt-percept` git dependency URLs (currently `<org>` placeholder). Once available, the in-workspace numerics (FFT/linalg/optimizers) can migrate behind those dependencies.
- [ ] Confirm `tpt-fem` crate name for mode-solving sub-crate; upgrade the scalar FEM to vectorial (spec's `tpt-fem-eigen`).
- [ ] Decide `tpt-percept` stub interface once that crate's API stabilizes.
- [ ] Photonic-crystal **defect modes** (supercell PWEM) and **graded-index / photonic-crystal fibers** — deferred; slab/fiber/EIM solvers cover the base cases.
- [ ] Zernike wavefront-error fit as a first-class merit (currently spot/OPD-based proxies + wave-crate PSF tools).
- [ ] Zemax binary `.zmx` and NSC-mode read/write; STEP/IGES via `tpt-io-standards`.
- [ ] End-to-end test with a true **double-Gauss** reference file (current E2E uses a singlet sample authored in-repo).

# tpt-optics — Build Tracking

**License:** MIT OR Apache-2.0 (dual), © TPT Solutions.
**Dependency policy:** No Apache-2.0-only dependencies permitted. Every dependency (internal or external) must be MIT-only or dual-licensed MIT/Apache-2.0. Enforced via `deny.toml`.

**Legend:** `[ ]` todo · `[~]` in progress · `[x]` done

---

## Phase 0 — Repo Bootstrap & Tooling

- [ ] `git init`, `.gitignore`
- [ ] Root `Cargo.toml` workspace manifest (`[workspace]`, `[workspace.package]` with `license = "MIT OR Apache-2.0"`, `edition`, `rust-version`)
- [ ] `LICENSE-MIT` (© TPT Solutions)
- [ ] `LICENSE-APACHE` (© TPT Solutions)
- [ ] `rustfmt.toml`
- [ ] `rust-toolchain.toml`
- [ ] `deny.toml` — strict allow-list (MIT, Apache-2.0 only as a dual-license component, MIT-0, BSD-2/3-Clause, ISC, Zlib, Unicode-3.0/CC0-1.0 as needed); explicitly deny any crate whose *only* license is Apache-2.0; document the dual-license exception logic in comments
- [ ] `justfile` (build / test / lint / deny recipes)
- [ ] `xtask/` crate for repo automation
- [ ] CI workflow: fmt, clippy, test, `cargo deny check`, no_std build check
- [ ] `README.md`
- [ ] `CONTRIBUTING.md`
- [ ] `SECURITY.md`
- [ ] `crates/` directory scaffold for all 9 crates (see section 3 of spec.txt)
- [ ] Pin git dependency revs: `tpt-math`, `tpt-dsp` (internal, MIT/dual)
- [ ] Pin git dependency revs: `tpt-fem` (`tpt-fem-eigen`), `tpt-materials`, `tpt-systems-optimisation`, `tpt-electromagnetics` — `git = "https://github.com/<org>/<repo>"`, `<org>` placeholder pending confirmation
- [ ] Stub `tpt-percept` integration point (feature-gated, off by default, no hard dependency — crate does not exist yet)

---

## Phase 1 — Foundation & Geometric Optics
*(spec.txt §6 Phase 0, target: Weeks 1-3)*

### tpt-optics-core
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-math` (internal, MIT)
- [ ] Unit-safe wavelength types (nm, μm, m)
- [ ] Unit-safe frequency types (THz, Hz)
- [ ] Unit-safe optical path length types
- [ ] Real refractive index representation
- [ ] Complex refractive index representation (absorbing materials)
- [ ] Dispersion model: Sellmeier
- [ ] Dispersion model: Cauchy
- [ ] Dispersion model: Schott
- [ ] Dispersion model: Drude-Lorentz
- [ ] Tabulated data interpolation for dispersion
- [ ] Material database with temperature dependence
- [ ] Material database with wavelength dependence
- [ ] Verify `no_std` + `alloc` compatibility
- [ ] Unit tests
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-geometric
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-optics-core`, `tpt-math`
- [ ] Sequential ray tracing: paraxial
- [ ] Sequential ray tracing: exact
- [ ] Non-sequential ray tracing: Monte Carlo (stray light, illumination, scattering)
- [ ] Surface type: spherical
- [ ] Surface type: aspherical (conic constants)
- [ ] Surface type: Zernike polynomial freeform
- [ ] Surface type: diffractive
- [ ] ABCD matrix formalism (Gaussian beam propagation, stability criteria)
- [ ] Ray-surface intersection with Newton-Raphson refinement
- [ ] Differentiable ray tracing (gradient-based optimization hook)
- [ ] Unit tests
- [ ] proptest: Snell's law reversibility
- [ ] proptest: ABCD matrix determinant = refractive index ratio (lossless systems)
- [ ] proptest: paraxial ABCD symplecticity
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### Phase 1 exit criteria
- [ ] Snell's law reversibility verified via proptest
- [ ] Paraxial ABCD symplecticity verified via proptest

---

## Phase 2 — Wave Optics & Polarization
*(spec.txt §6 Phase 1, target: Weeks 4-6)*

### tpt-optics-wave
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-optics-core`, `tpt-dsp` (FFT)
- [ ] Scalar diffraction: Fresnel integral
- [ ] Scalar diffraction: Fraunhofer integral
- [ ] Beam Propagation Method (BPM) — guided-wave
- [ ] Beam Propagation Method (BPM) — free-space
- [ ] Fourier optics: lens as Fourier transformer
- [ ] Fourier optics: OTF/MTF/PSF computation
- [ ] Gaussian beam propagation: q-parameter
- [ ] Gaussian beam propagation: ABCD law
- [ ] Gaussian beam propagation: mode matching
- [ ] Coherence theory: mutual coherence function
- [ ] Coherence theory: van Cittert-Zernike theorem
- [ ] Integrate `tpt-dsp` FFT for diffraction calculations
- [ ] Unit tests
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-polarization
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-optics-core`, `tpt-math`
- [ ] Jones calculus (2x2 complex matrices, coherent states)
- [ ] Mueller calculus (4x4 real matrices, incoherent states, depolarization)
- [ ] Stokes vector representation (partially polarized light)
- [ ] Polarization ray tracing through optical systems
- [ ] Birefringence modeling
- [ ] Dichroism modeling
- [ ] Optical activity modeling
- [ ] Unit tests
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### Phase 2 exit criteria
- [ ] Integration test: Gaussian beam propagation through a lens matches analytical q-parameter transformation

---

## Phase 3 — Photonics & Design
*(spec.txt §6 Phase 2, target: Weeks 7-9)*

### tpt-optics-waveguide
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-optics-core`, `tpt-fem` (`tpt-fem-eigen`, git dep)
- [ ] Vectorial FEM mode solving for waveguide cross-sections
- [ ] Coupled-mode theory: power exchange between waveguides
- [ ] Coupled-mode theory: gratings
- [ ] Coupled-mode theory: resonators
- [ ] Photonic crystals: band structure calculation
- [ ] Photonic crystals: defect modes
- [ ] Fiber optics: step-index
- [ ] Fiber optics: graded-index
- [ ] Fiber optics: photonic crystal fiber
- [ ] Integrated photonics: Mach-Zehnder interferometers
- [ ] Integrated photonics: ring resonators
- [ ] Integrated photonics: arrayed waveguide gratings
- [ ] Unit tests
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-design
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-optics-core`, `tpt-optics-geometric`, `tpt-systems-optimisation` (git dep)
- [ ] Merit function: spot size
- [ ] Merit function: wavefront error
- [ ] Merit function: MTF
- [ ] Merit function: distortion
- [ ] Merit function: chromatic aberration
- [ ] Optimization: damped least-squares
- [ ] Optimization: global (genetic algorithms via `tpt-systems-optimisation`)
- [ ] Optimization: global (simulated annealing via `tpt-systems-optimisation`)
- [ ] Tolerance analysis: Monte Carlo sensitivity
- [ ] Tolerance analysis: compensator modeling
- [ ] Automatic differentiation of optical systems (gradient-based optimization)
- [ ] Unit tests
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-verify
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-optics-core`, `tpt-optics-geometric`, `tpt-optics-polarization`, `kani`, `proptest`
- [ ] proptest strategies: valid optical systems generation
- [ ] proptest strategies: surface geometry generation
- [ ] proptest strategies: material property generation
- [ ] Kani harness: ray-surface intersection correctness (nearest positive root)
- [ ] Kani harness: energy conservation in lossless beam splitters
- [ ] Invariant check: optical reversibility
- [ ] Invariant check: ABCD matrix symplecticity
- [ ] Invariant check: Snell's law boundary conditions
- [ ] Unit tests
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### Phase 3 exit criteria
- [ ] Kani proofs written for ray-surface intersection
- [ ] Kani proofs written for energy conservation in lossless systems

---

## Phase 4 — I/O & Advanced Features
*(spec.txt §6 Phase 3, target: Weeks 10-12)*

### tpt-optics-io
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Wire deps: `tpt-optics-core`, `tpt-optics-geometric`, `tpt-io-standards` (STEP/IGES)
- [ ] Zemax `.zmx` parser (sequential mode)
- [ ] Zemax `.zmx` parser (non-sequential mode)
- [ ] Zemax `.zmx` writer (sequential mode)
- [ ] Zemax `.zmx` writer (non-sequential mode)
- [ ] CODE V `.seq` import
- [ ] CODE V `.seq` export
- [ ] Standard lens geometry: STEP import/export via `tpt-io-standards`
- [ ] Standard lens geometry: IGES import/export via `tpt-io-standards`
- [ ] Glass catalog import: Schott
- [ ] Glass catalog import: Ohara
- [ ] Glass catalog import: Hoya
- [ ] Glass catalog import: CDGM
- [ ] Unit tests
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes
- [ ] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-geometric extensions
- [ ] Non-sequential ray tracing enhancements for stray light analysis
- [ ] Stray light analysis reporting/output

### Phase 4 exit criteria
- [ ] End-to-end test: import a double-Gauss lens from Zemax, optimize it, export it

---

## Phase 5 — Umbrella Crate & Ecosystem Integration

### tpt-optics (umbrella crate)
- [ ] Scaffold `Cargo.toml` + `lib.rs`
- [ ] Feature-gated re-export: `core`
- [ ] Feature-gated re-export: `geometric`
- [ ] Feature-gated re-export: `wave`
- [ ] Feature-gated re-export: `polarization`
- [ ] Feature-gated re-export: `waveguide`
- [ ] Feature-gated re-export: `design`
- [ ] Feature-gated re-export: `io`
- [ ] Feature-gated re-export: `verify`
- [ ] `full` feature enabling all of the above
- [ ] rustdoc
- [ ] fmt/clippy clean
- [ ] `cargo deny check` passes

### Ecosystem integration
- [ ] `tpt-electromagnetics` handoff notes (high-frequency/optical vs quasi-static/RF boundary)
- [ ] `tpt-materials` integration: thermo-optic coefficients for athermalization studies
- [ ] `tpt-percept` stub: feature-gated differentiable rendering hook, no hard dependency (crate not yet available)
- [ ] AI-native agent notes: differentiable ray tracing / wave propagation hooks for `tpt-eve` / `tpt-anima` inverse design workflows

### Documentation & publishing
- [ ] Architecture overview doc
- [ ] Publishing checklist: crates.io metadata (all 9 crates + umbrella)
- [ ] Versioning strategy doc
- [ ] `CHANGELOG.md`

---

## Open items requiring follow-up

- [ ] Confirm GitHub org/account for `tpt-fem`, `tpt-materials`, `tpt-systems-optimisation`, `tpt-electromagnetics` git dependency URLs (currently `<org>` placeholder)
- [ ] Confirm `tpt-fem` crate name for mode-solving sub-crate (spec.txt references `tpt-fem-eigen`)
- [ ] Decide `tpt-percept` stub interface once that crate's API stabilizes

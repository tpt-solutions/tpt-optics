# tpt-optics — Build Tracking

**License:** MIT OR Apache-2.0 (dual), © TPT Solutions.
**Dependency policy:** No Apache-2.0-only dependencies permitted. Every dependency (internal or external) must be MIT-only or dual-licensed MIT/Apache-2.0. Enforced via `deny.toml`.

**Legend:** `[ ]` todo · `[~]` in progress or partly done · `[x]` done · `[-]` rejected or out of scope

> **Status (2026-10-07):** all phases implemented. Gate: `just check`
> (fmt · clippy `-D warnings` · 208 tests · `cargo deny` · `no_std` builds) — green.
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
- [x] Fiber optics: graded-index — power-law profile, scalar finite-difference radial solver (`graded` module; checked against the parabolic oscillator formula and the step-index limit)
- [~] Fiber optics: photonic crystal fiber — *deferred* (open item)
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
- [x] Merit function: wavefront error — Zernike least-squares fit of pupil OPD (`tpt-optics-design::wavefront`: `fit_zernike`, `wavefront_rms`)
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
- [x] Property test: ray-surface intersection correctness (nearest positive root) — `tests/properties.rs` (Kani removed later, see CHANGELOG)
- [x] Property test: energy conservation in lossless beam splitters (was Kani, now `tests/properties.rs`)
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
- [x] Intersection properties written (`tests/properties.rs::sphere_hit_is_ahead_and_on_surface`, `sphere_returns_the_nearest_positive_root`; originally Kani)
- [x] Energy-conservation property written (`tests/properties.rs::lossless_split_conserves_energy_below_tir`; originally Kani)

---

## Phase 4 — I/O & Advanced Features

### tpt-optics-io
- [x] Scaffold `Cargo.toml` + `lib.rs`
- [x] ~~Wire deps: `tpt-io-standards` (STEP/IGES)~~ — STEP/IGES deferred with the unavailable `tpt-io-standards` crate (deviation note); text formats implemented
- [x] Zemax `.zmx` parser (sequential mode, text subset: SURF/CURV/CONI/DISZ/GLAS/DIAM/STOP/NAME/MODE)
- [-] Zemax `.zmx` parser (non-sequential mode) — *rejected by design*: binary/NSC .zmx unsupported, documented (`MODE NSC` is a parse error)
- [x] Zemax `.zmx` writer (sequential mode, round-trip verified)
- [~] Zemax `.zmx` writer (non-sequential mode) — *deferred* (open item)
- [x] CODE V `.seq` import (LEN NEW/SO/S/STO subset)
- [x] CODE V `.seq` export
- [~] Standard lens geometry: STEP import/export via `tpt-io-standards` — *deferred* (open item)
- [~] Standard lens geometry: IGES import/export — *deferred* (open item)
- [x] Glass catalog import: Schott (via AGF subset + embedded reference glasses)
- [~] Glass catalog import: Ohara — *via AGF parser* (vendor file required; no embedded data)
- [~] Glass catalog import: Hoya — *via AGF parser* (vendor file required; no embedded data)
- [~] Glass catalog import: CDGM — *via AGF parser* (vendor file required; no embedded data)
- [x] Unit tests (parse/write round-trips, positioned errors, convert-to-trace)
- [x] rustdoc
- [x] fmt/clippy clean
- [x] `cargo deny check` passes
- [x] Registered in workspace `Cargo.toml` + umbrella crate

### tpt-optics-geometric extensions
- [x] Non-sequential ray tracing enhancements for stray light analysis (scatter models, power bookkeeping, max-bounce budget)
- [x] Stray light analysis reporting/output (`NsResult::report`)

### Phase 4 exit criteria
- [~] End-to-end test: import a double-Gauss-like lens from Zemax text, optimize it, export it — realized as the singlet pipeline (the double-Gauss reference is still open, see below): `.zmx` parse → `to_sequential` (catalog glass resolution) → trace/`rms_spot` → DLS best-focus optimization → `from_sequential` export → re-parse (`zmx` tests + `design` best-focus test; the sample is a singlet — see open item)

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

- [~] Confirm GitHub org for sibling crates: `github.com/tpt-solutions/` (confirmed; the migration onto those crates is not done) — `tpt-fem`, `tpt-materials`, `tpt-systems-optimisation`, `tpt-electromagnetics`, `tpt-io-standards` and `tpt-perception` exist there. The `tpt-fem` workspace lists `crates/tpt-fem` (umbrella) and `crates/tpt-fem-eigen` as members. The in-workspace numerics (FFT/linalg/optimizers) can migrate behind those dependencies once confirmed.
- [x] Confirm `tpt-fem` crate name: `crates/tpt-fem` (umbrella) and `crates/tpt-fem-eigen` are workspace members. Vectorial upgrade done in-workspace: see the `vector_fem` module (hybrid edge–nodal, Nédélec + P1). It checks against the PEC TE10 cutoff (rel. 2e-4 at 16×8 cells), the exact dielectric-slab TE0 index (rel. 1e-4 at 4×64), and the scalar solver on a weak strip (6.6e-4 apart). It is dense, so practical for a few hundred edges.
- [x] Decide `tpt-percept` stub interface (decision: no dependency; the differentiable-tracing hook stays in this repo). The real crate is `tpt-solutions/tpt-perception` (`tpt-percept-*` sub-crates: core, cloud, features, register, slam, fusion, map). It is a 3-D perception/SLAM library with no gradient or `ParamId` API, so the differentiable-tracing hook stays in this repo with no dependency on it. Integration with `tpt-perception` is deferred (decision: none for now).
- [x] Photonic-crystal **defect modes** (supercell PWEM) — E-polarisation (TM) supercell solver with rod-defect Fourier coefficients and a full-zone gap helper. Checked for square rods (ε = 8.9, r/a = 0.2, M = 5): the TM gap is (0.324, 0.447) c/a, and a missing-rod defect gives a localised level at 0.393 c/a, converged to ~0.3% for n_g = 17–29 in the 5×5 cell. The 7×7 cell shifts it by ~1.3% (0.3996 at n_g = 21), so the isolated-defect level is 0.39–0.40 c/a. H-polarisation (TE) defect levels are not quantitative (slow 1/ε convergence), and high-permittivity defects (ε = 40) are not converged. Eigenproblems use Householder + implicit-QL (`eigen` module) instead of Jacobi.
- [ ] **Photonic-crystal fibers** — deferred (decision: keep deferred until a fibre use case is scoped); slab/fiber/EIM/graded-index solvers cover the base cases.
- [x] Zemax binary `.zmx` and NSC-mode read/write; STEP/IGES — decided, out of scope for this build. `tpt-io-standards` has no STEP or IGES API (its crates cover medical, GIS, EDA and other formats), so STEP/IGES would need a from-scratch ISO 10303 / IGES implementation against specs not in this repo. Binary and NSC `.zmx` stay unsupported (documented). Revisit only if a concrete use case and spec source are provided.
- [ ] End-to-end test with a true **double-Gauss** reference file (current E2E uses a singlet sample authored in-repo). Decision: stays open until an external source file is supplied; an in-repo prescription would not be an independent reference.

---

## Review follow-up (2026-10-31)

Findings from a read-only review. Each item got a regression test before its fix.

### Phase 1 — bugs (high)
- [x] `waveguide/supercell.rs` `bulk_gap_from_supercell`: validates `n_k >= 1` and the window, returns `Result`, and rejects non-finite frequencies instead of dropping them.
- [ ] `waveguide/fem.rs` `fem_fundamental_mode`: fixed `0.97·(k0·n_max)²` shift can converge to the wrong mode for strongly guiding cores; use the unused `fundamental_shift`. **Open, not reproduced:** probes on strong (n_core 3.48) and weak guides returned the fundamental with either shift. `fundamental_shift` sits below the guided spectrum, so it is not a safe drop-in without a case that reproduces the failure.
- [x] `waveguide/fem.rs` `fem_fundamental_mode`: unit-max normalisation uses `max`, not max-abs, so it depends on iterate sign; fixed by `normalise_field`.
- [x] `io/glass.rs` AGF parser: second `CD` line overflows `next_cd_slot` and panics; `filter_map(.ok())` shifts coefficients; `nan`/`inf` accepted. Now: `CD` values accumulate, more than six is an error, and bad tokens are errors.
- [x] `io/glass.rs` AGF parser: dispersion-formula number ignored (non-Sellmeier glasses silently wrong); `TD`/`CC` documented but unhandled; check `ED` index layout. Now: formula read from `NM`; 1 (Schott) and 2 (Sellmeier) supported, others are errors. `TD`/`CC` are skipped and documented. The `ED` index layout is still unverified against a vendor file.
- [x] `io/zmx.rs` `write()`: `{:.12}` truncates higher-order asphere coefficients; use `{:e}`. Also index hints use `{v}` (exact) and records no longer end with a blank line.
- [x] `io/zmx.rs` parser: TYPE/PARM (aspheres), STOP and object DISZ silently dropped; warn or error on unsupported keywords. Now: `TYPE` other than `STANDARD` and non-zero `PARM` are errors. `STOP` is kept but does not change the trace, and object `DISZ` is kept but unused (both documented).

### Phase 1 — bugs (medium)
- [x] `vector_fem`: add input validation, convergence/residual reporting and a `Result` return; detect singular LU pivots; guard dense allocation size; stop `sqrt(|β²|)` turning evanescent results into a real `n_eff`. Now returns `VectorMode` (`n_eff`, `residual`, `converged`) or `VectorFemError`; meshes over `MAX_DENSE_UNKNOWNS` are refused before allocation.
- [x] `supercell`: validate `m`, `n_g` oddness, rod symmetry/overlap/radius and epsilon; report loss of positive-definiteness in the Cholesky step instead of clamping.
- [x] `phc::jacobi_eigen` absolute tolerances never converge on `graded.rs` matrices (~1e14); switch `graded.rs` to tridiagonal `eigen::implicit_ql` and time it in release. Done: release test time 0.18 s → 0.01 s. `jacobi_eigen` itself is unchanged and still used by the 2-D PWEM.
- [x] `coupled.rs`: `ring_resonances` unbounded loop, `grating_spectrum` sinh² overflow gives NaN, divide-by-zero in `coupling_length`/`ring_q`/`ring_fsr`/`ring_lorentzian`. Now: `ring_resonances` returns `Result` and refuses windows over `MAX_RESONANCES`; the divide helpers return `Option`; `grating_spectrum` saturates.
- [x] `io`: `seq.rs` accepts nan/inf; `seq` vs `zmx` last-medium mismatch (Vacuum vs Const 1.0); mirror before Vacuum lost on `zmx` export; `Const(v)` written at `{:.6}`. Now: all fixed, and both readers end with air.
- [x] `io`: handle UTF-8 BOM and UTF-16 input; implement `Error` for `ZmxError`/`SeqError`/`AgfError` without the `std` feature. New `io::text::decode`; parsers strip a UTF-8 BOM; errors implement `core::error::Error` unconditionally (stable since 1.81).
- [x] `core`: `MaterialDatabase::get` is case-sensitive and O(n); `Tabulated::from_points` panics on empty/duplicate points; Edlén pole domain; `fresnel_energy` NaN at 90°; `bessel_*` contracts only `debug_assert`. Now: `BTreeMap` with lower-cased keys; `from_points` returns `TableError`; Edlén is NaN outside 0.2–2 µm; Fresnel energy is finite at 90°; Bessel checks are `assert!` in every profile.
- [x] README: root `README.md:48` asserts BK7 n = 1.51472 ± 1e-4 vs 1.5151 in `crates/tpt-optics/README.md:73`; run the snippet, fix it, and make README snippets doctests. Root value corrected to 1.5151 (Sellmeier BK7 gives n ≈ 1.51509); both READMEs are now doctests.

### Phase 1 — low / infra
- [x] Doc/code mismatches: `fem.rs` says CG but uses LU; garbled `phc.rs` header; `devices.rs` "field" vs power; `coupled.rs` Lorentzian vs dip; `Arc::new(())` hack in `glass.rs`; double newline in `zmx::write`. Fixed the `fem.rs` wording, the `coupled.rs` Lorentzian doc, the `Arc` placeholder and the zmx newlines. The `phc.rs` header and `devices.rs` wording no longer contain the reported text.
- [x] Deduplicate dense LU/eigen helpers across `fem`, `vector_fem`, `phc`, `eigen`, `graded`. Done for LU: one `LuFactor` and one `mat_vec` in the private `dense` module, shared by `fem` and `vector_fem`. `eigen` (Householder + QL) and `phc` (Jacobi) are different algorithms and are unchanged.
- [x] CI: add Windows + macOS matrix and an MSRV 1.81 job; real `thumbv7em-none-eabi` no_std build covering wave/waveguide/design; commit `Cargo.lock`. Done: test matrix (Windows, macOS), `thumbv7em-none-eabi` build of the six no_std crates, an MSRV job, an examples job, and `Cargo.lock` committed. The MSRV build needed `proptest` pinned to 1.6 and `tempfile` to 3.20 (dev-dependencies, to keep their newer rustc requirements out of the resolve), and `f64::abs`/`signum` replaced by `libm` equivalents. The lock file pins those versions.
- [x] Benches (criterion) for trace/FFT/BPM/FEM. Three bench targets: `kernels` (FFT 4096, BPM 256×200), `trace` (two-surface slab), `modes` (scalar and vector FEM on a 16×16 mesh). Run with `cargo bench`. No baseline numbers are recorded yet.
- [x] Missing tests: supercell with a real defect, `eigen` edge cases (n=0/1, repeated, NaN), `vector_fem` failure paths, vacuous `eim::narrow_guide_near_cutoff`, io BOM / duplicate-CD / asphere precision / mirror+Vacuum. Supercell defect coverage was already present; the others were added this round.

### Phase 2 — adoption
- [x] `crates/tpt-optics/examples/` 01–10. All ten run in release, and the CI `examples` job runs each one. Scope changes from the list: 03 is a singlet spot diagram, not a Cooke triplet (no verified triplet prescription in-repo); 06 uses Jones calculus only; 07 imports a singlet .zmx, optimises its front radius to a 100 mm EFL, and exports, not a double-Gauss. The examples use illustrative data, and each one says so.
- [~] `presets` module (plano-convex, singlet, doublet achromat, Cooke triplet, double Gauss, SMF-28, ring resonator, PC waveguide) and a `SystemBuilder` over `SequentialSystem::push`. Done: `geometric::builder::SystemBuilder`; `geometric::presets::{plano_convex, singlet, achromatic_doublet}` (the achromat is solved by Newton iteration on the traced d-line focal length and F–C colour, to 0.02 mm); `waveguide::presets::{SMF28, fiber_modes, SiliconRing, pc_line_defect_rods}`. Done: `cooke_triplet_50mm`. Open: the double Gauss. An attempt to design a triplet by optimisation from a rough start drove air gaps to zero, so the result was degenerate and I did not keep it. A verified starting prescription, with bounded optimisation, is needed first.
- [x] Docs: `docs/GUIDE.md` (first lens in 10 lines), cookbook, `docs/VALIDATION.md` generated from golden tests, zmx/CODE V compatibility table, per-crate quickstart and feature matrix; add `design` and `io` to the prelude; document `just`/`cargo-deny`/Kani prerequisites. Done: `docs/GUIDE.md`, `docs/COOKBOOK.md`, `docs/COMPATIBILITY.md`, `docs/PREREQUISITES.md` (just, cargo-deny, the wasm and Python toolchains, Kani), and `docs/VALIDATION.md`. Every crate README is a doctest, so its quick start is compiled. The umbrella prelude includes the design and io items. `cargo run -p xtask -- validation` fails if a test cited in `docs/VALIDATION.md` no longer exists, and CI runs it. The validation page itself is written by hand, and the checker keeps it honest.
- [ ] Complete `docs/PUBLISHING.md` checklist, release 0.1.0, verify docs.rs with `--features full`; revisit "issues only" in `CONTRIBUTING.md`.

### Phase 3 — usability and automation
- [x] Crate-level typed error enum with kind/context/hints; blocked-ray reasons (miss, TIR, aperture); remove `unwrap_or(f64::NAN)` in `sequential.rs`; "did you mean" for unknown glass names. Done: `tpt_optics::Error` wraps each crate's error type, with `kind()` (`ErrorKind::{Parse, Input, Glass, Numerics}`), `hint()`, and `source()`. `entrance_pupil` returns `Option` for the F/#. Traces report `BlockReason` and total internal reflections. Unknown glasses suggest a close name.
- [x] Units: `Length`, `Curvature`, `Power`; consistent degrees vs `Angle`; auto-scaling `Display`. Done: `core::units::{Length, Curvature, Power}` with unit conversions and an auto-scaling `Display` (m, mm, µm, nm). `Angle` already existed. Field-angle arguments are still `f64` degrees, named `*_deg`, so the public API is consistent but not typed.
- [~] Optional `serde` with TOML/JSON schema for systems, optimizer settings, tolerances and results; `tpt-optics` CLI (`info|glass|trace|convert|analyze|optimize|run design.toml`). Done: the CLI (`tpt-optics-cli`) with all those commands; the design file is serde with a JSON Schema (`tpt-optics schema`); results export as JSON; `tpt-optics-design` has an optional `serde` feature covering analysis reports, operands, optimiser results, tolerance statistics and glass-map results. Open: `SequentialSystem` itself is not serialisable, because its media hold shared material handles. Its serde form needs a design decision (serialise the material name and resolve it on load).
- [x] `system.analyze()` bridge (spot, ray fan, OPD, field curvature, distortion, MTF, wavefront) plus SVG plots and CSV/JSON export. Done: `design::analysis::analyze` gives EFL, BFL, F/#, chromatic shift, and per field the spot, ray fan, OPD relative to the chief ray, distortion, wavefront RMS, best-focus shift (a through-focus search), and the system MTF along the tangential and sagittal axes (`design::mtf`, pupil zero-padded to 2n, matching the Hopkins curve to 0.03). Exports: CSV, JSON, an SVG spot diagram, and an SVG ray-fan plot.
- [x] Sweep helpers (grid/LHS/Sobol) and optional `parallel` (rayon) feature with splittable `Rng` so Monte Carlo is bit-identical at any thread count. Done: `design::sweep::{grid, latin_hypercube, sobol}`; the `parallel` feature adds `tolerance::analyze_parallel`, which matches `analyze_seeded` exactly for the same seed. Each run has its own seed derived from the base seed and its index, so the result does not depend on scheduling.
- [x] Optimizer: history, convergence reason, callback/cancel, named parameters. Done: `OptimizeResult::{history, stop}` (`Converged`, `MaxIterations`, `Stalled`, `Cancelled`); `damped_least_squares_observed` with a per-iteration observer that can cancel; `design::params::{ParamSet, optimize_named}` for named parameters.
- [~] Materials: bundled Schott/Ohara/Hoya/CDGM AGF data (check licences vs `deny.toml`), Abbe/glass-code metadata, validity-range warnings. Done: `Material::{abbe_number, glass_code, with_validity, outside_validity}`; BK7 and fused silica carry their fitted ranges; the CLI warns outside them. Open: the bundled vendor data. It needs the catalogue files and a licence decision, which I do not have here.

### Phase 4 — innovation
- [x] Exact-gradient DLS via `Dual`/`gradient_of`; first-order tolerance sensitivity from one dual pass. Done: exact operand derivatives for EFL, BFL (`paraxial_focus_gradients`, dual-number paraxial trace), RMS spot (per-ray `gradient_of` combined through the centroid), and ray height; `design::gradient::optimize_exact` (Levenberg–Marquardt with the exact Jacobian); `tolerance_budget` (first-order standard deviation from one Jacobian, within 15% of Monte Carlo). Thickness and edge-thickness operands use finite differences, which are exact enough for piecewise-linear and sag terms.
- [x] Operand-based merit functions (EFL, thickness, edge, ray-height constraints, weights). Done: `design::operands` with EFL, BFL, RMS spot, thickness range, ray height at a surface, and edge thickness. Each is a weighted residual, and a missing quantity gets a finite penalty. `optimize_system` runs Levenberg–Marquardt on a builder of the system from parameters.
- [ ] `Scalar`-generic wave/waveguide crates for end-to-end differentiable inverse design. Not started. The wave and waveguide solvers (FFT, BPM, FEM, slab, EIM) are written on `f64`, and making them generic is a large refactor with no test for the derivatives yet.
- [~] `system.self_check()` exposing verify-crate invariants; surface Kani proofs in docs. Done: `design::check::self_check` (five structural checks: finite curvatures, positive apertures, finite gaps, paraxial determinant equals the index ratio, finite on-axis path). `tpt_optics_verify::system_check::verify_system` runs the trace invariants on three rays. Reversibility is checked on axis, at 2 mm off axis and at a 2° tilt. The retrace started on the vertex plane, which ignored the sag of a curved last surface and left about 20 µm of path off axis. It now starts at the traced point on that surface, and a singlet passes all six checks. `design` does not depend on `verify` (which pulls in proptest), so the two checks are separate. Kani proofs are described in `docs/PREREQUISITES.md`, but not surfaced in the self-check.
- [x] Run manifests (version, seed, input hash) and content-addressed caching for sweeps. Done: `tpt-optics run` writes `manifest.json` (tool version, input key, outputs, and the hash algorithm) and caches results under `<dir>/cache/<key>`, keyed by an FNV-1a hash of the command and the input bytes. The run has no random seed, because it is deterministic. Sweeps record their seeds in their own calls, which the manifest does not yet include.
- [~] WASM interactive lens-design demo; PyO3 Python bindings. Done: `crates/tpt-optics-wasm` builds for `wasm32-unknown-unknown` and with `wasm-pack`; the Node smoke test passes; the browser page `www/index.html` is written but I have not opened it in a browser. `crates/tpt-optics-py` is a PyO3 module, built and imported in a scratch directory. Per the recommendation, it stays outside the workspace, because pyo3's build dependency is Apache-2.0 with LLVM-exception, which the licence policy rejects. CI builds both.
- [x] Glass-map auto-design (GA + catalog + operands); MCP/JSON-RPC server over the CLI. Done: `design::glassmap::rank_achromats` (an exhaustive ranking) and `genetic_achromat` (a seeded genetic search over glass choice and element thicknesses, with the traced achromat as the objective); `tpt-optics mcp` serves `glass_index`, `analyze` and `optimize` over JSON-RPC on stdio. The ranking and search are over the glasses given, not over shapes or a full design space.

---

## Review follow-up (2026-10-31, second pass)

Source: static read-only review of all crates. Each fix needs a regression test, and the items below are ticked only once a test covers them.

### Phase 1 — silent wrong answers (before the 0.1.0 publish)
- [x] `io/zmx.rs:129-139`: `CURV` extra args are read as asphere terms (they are solve/pickup fields); `UNIT` ignored (inch files 25.4× wrong); quoted/empty tokens fail; `TYPE EVENASPH` + `PARM` unsupported. Done: `CURV` takes curvature only; `UNIT` (MM/CM/M/IN) scales lengths, curvature and asphere terms to mm on read; `NAME`/`GLAS` accept quotes; `TYPE EVENASPH` reads `PARM 1`–`7` into `asphere[1..8]` (`asphere[0]`, the r² term, has no PARM). The writer emits `TYPE EVENASPH` + `PARM` lines in place of the old asphere-in-`CURV` form. `PARM 8` (r¹⁸) is an error, since the crate stores up to r¹⁶.
- [x] `io/seq.rs`: common CODE V commands (`DIM`, `WL`, `EPD`, `YAN`, `RDY`, `THI`, `GLA`) error though docs say unknown ones are tolerated; `" ref"` arm can never match; inline `!` comments not stripped; `to_sequential` drops surface 0 when `SO` is absent; mirror medium disagrees with the zmx reader. Done: `DIM`/`WL`/`EPD`/`YAN` are tolerated; `RDY`/`THI`/`GLA` are errors because they change geometry; the dead arm is removed; trailing `!` comments are stripped; an `S` before `SO` is an error; mirrors keep the travelling medium, as in zmx.
- [x] `io/glass.rs`: one unsupported dispersion formula (3–13) aborts the whole AGF catalog (skip and collect warnings); `ED` field 5 misread as dn/dT, `TD` line skipped; `LD` validity ignored; wrong line number on the last glass. Done: `parse_agf_report` returns the catalogue plus skipped glasses (each with its `NM` line); `parse_agf` returns the database and drops those with a documented note; `TD` D0 (×10⁻⁶ K⁻¹) gives dn/dT; `ED` no longer does; `LD` sets the validity range. The `ED` field layout is still unverified against a vendor file.
- [x] `io/text.rs`: add cp1252/Latin-1 fallback; UTF-32 LE BOM is detected as UTF-16. Done: non-UTF-8 input decodes as Windows-1252 (never fails); UTF-32 LE/BE BOMs are decoded, and the UTF-32 mark is checked before UTF-16. `DecodeError::InvalidUtf8` is removed.
- [x] `design/operands.rs:~141`: `RayHeight` returns `p.x`, always 0 for the chief ray; add an axis selector or use y. Done: `RayAxis::{X, Y}` on `Operand::RayHeight`, defaulting to `X` so earlier files keep their meaning; the CLI `ray_height` takes `axis = "y"`. The gradient follows the same axis.
- [x] `design/tolerance.rs:~170`: `summarize` panics when `passes == 0`; out-of-range tolerance index panics. Validate and return `Result`. Done: `analyze`, `analyze_seeded` and `analyze_parallel` return `Result<_, ToleranceError>` (`NoRuns`, `IndexOutOfRange`).
- [x] `SequentialSystem::validate()`: mismatched `media`/`surfaces`/`thickness` lengths and a non-final infinite thickness panic or give NaN. Call it from trace entry points, `spec.rs` and the CLI. Also `TunedParams::seeded` with an out-of-range `ParamId`. Done: `validate()` returns `SystemError`; `trace_with` traces an invalid system as blocked with no points; the zmx and seq readers report it as a parse error. The CLI and `SystemBuilder` already reject non-finite gaps. `TunedParams::seeded` returns `Option`, and `gradient_of` gives NaN for an unknown `ParamId`.
- [x] `geometric/nseq.rs:~232`: `launched_power = ray_count` ignores per-ray source power. Done: `launched_power` sums the power of each ray the source launches.

### Phase 2 — correctness of analysis (off-axis, fans)
- [~] Fans are all sagittal (`analysis.rs:288-340`, `sequential.rs:722-763`): add a 2-D hexapolar pupil grid with area weighting and a tangential fan. Done: `FieldReport::ray_fan_tangential` (pupil y, image y), so the report has both fans. Open: the 2-D hexapolar grid with area weights, which needs the tangential rays in the spot and the RMS.
- [~] `design/wavefront.rs`: hard-coded `start_z = -0.1`, no ray aiming, so off-axis wavefront/MTF are wrong; `analysis.rs:~315` uses the on-axis `chief_opd` for every field. Done: the off-axis OPD now uses each field's own chief ray (fixed in analysis.rs). Open: ray aiming and the start plane.
- [ ] Stop assumed at surface 0, object at infinity (`field_ray`, `distortion`, `entrance_pupil`, `analyze`, `system_mtf`): add `Aperture`/`Field`/`Stop` specs, ray aiming, finite objects, vignetting, wavelength weights.
- [ ] `sequential.rs:516-545`: TIR in `Diffract` uncounted; evanescent order silently refracts (add a `BlockReason`); `n_cur` not updated after reflect/TIR.
- [~] `design/optimize.rs:~280-310`: Jacobian recomputed after a rejected step; `evaluations` ignores the 2n Jacobian traces; NaN start cost undetected; short `bounds` silently zipped; no parameter scaling. Done: the Jacobian is computed once per accepted step; central-difference Jacobians add 2n to `evaluations` (a supplied Jacobian adds 0); a non-finite start stops as `Stalled`; short `bounds` panic. Open: parameter scaling.
- [~] `core/dispersion.rs:~85` Sellmeier: negative n² gives NaN; `index_at` ignores `outside_validity`. Done: n² < 0 gives n = 0, k = √(−n²) (an absorbing continuation). Not changed: `index_at` still returns the model value, and callers check `outside_validity`, which the CLI already warns on.
- [~] `design/merit.rs:~24` `chromatic_focal_shift`: doc promises a spot-minimising shift, body returns paraxial BFL spread; returns 0.0 silently with <2 valid wavelengths; no lateral colour. Done: the doc now says it is the paraxial BFL spread; fewer than two foci gives NaN, which the analysis report shows as not computed. Open: lateral colour.
- [x] `sequential.rs:516-545`: TIR in `Diffract` uncounted; evanescent order silently refracts (add a `BlockReason`); `n_cur` not updated after reflect/TIR. Done: a TIR at a diffractive surface is counted and keeps the ray in its medium; an evanescent order blocks the ray with the new `BlockReason::Evanescent`. Reflection already kept `n_cur`, so that part needed no change.
- [x] `core/fresnel.rs:46`: `theta_t` documented complex, imaginary part hard-coded 0. Done: θ₂ = −i·ln(cos θ₂ + i sin θ₂), so the imaginary part is −ln|z|, zero for propagating light and nonzero under TIR.
- [~] Low: `Rng::normal` discards the second sample and is unbounded (add truncated); `glass_code` with n_F = n_C (done: V_d is infinite and the code clamps to 999, so it stays six digits); `intersect_asphere` returns `None` on a negative quadric discriminant and has no damping or failure signal; `glassmap.rs:146` redundant `.min`; `fit_zernike`/`mtf_curves` unvalidated (power of two, n > ~170); rank-deficient fit silently reports raw RMS.

### Phase 2 — wave / waveguide / polarization
- [x] `wave/diffraction.rs`, `bpm.rs`: evanescent components hard-zeroed; use `exp(-k z √(-arg))`. Done: both propagators decay evanescent components as `exp(-k z √(-arg))`, which is 1 at z = 0, so a zero-distance propagation is the identity.
- [~] `angular_spectrum`: add band-limit (Matsushima) and 2N padding; validate `λ`, `n_idx`, `z`, `dx`. Done: validates λ, dx, n and z (panics on bad input). Open: the band-limit and padding.
- [ ] `fresnel_sfft`/`fraunhofer`: 1-D only vs 2-D `angular_spectrum`; add 2-D variants, `z > 0` check, chirp-sampling and far-field checks.
- [ ] `wave/bpm.rs`: no absorbing boundary (periodic wrap-around); `bpm_3d` clones the field and recomputes `h` per step; real-only index profile.
- [ ] `wave/coherence.rs:29-45`: `van_cittert_zernike` ignores `λ`/`z`, wrong DC index (spurious (-1)^k), no quadratic phase.
- [~] `wave/gaussian.rs:123-172`: `mode_match_thin_lens` minimises `w(d_out) − w0_out`, not the waist; `w()` NaN when `q.im == 0`. Done: `mode_match_thin_lens` now matches the waist size (`w0_equivalent`); the waist position follows from the lens, since one focal length cannot set both. Open: `w()` at q.im = 0.
- [x] `wave/fourier.rs`: `psf_from_pupil` loses Strehl and gives NaN for an empty pupil; `mtf_from_psf` clamps to [0,1]. Done: empty pupil panics with a message; `mtf_from_psf` no longer clamps; `strehl_ratio` gives the absolute peak against an ideal pupil, which `psf_from_pupil` does not.
- [ ] `wave/fft.rs`: twiddle recurrence drift, no plan cache, non power of two panics (Bluestein); three duplicate 2-D FFTs (`fft::fft2`, `diffraction::fft2_in_place`, `fourier::fft2_rows_cols`).
- [x] `waveguide/eigen.rs:~108`, `phc.rs:174-200`: QL and Jacobi stop at the iteration cap and return unconverged values; Jacobi tolerance is absolute; no `w.len()` check. Done: QL returns `EigenError` at its cap (`symmetric_eigenvalues`, `graded_lp_modes` and the supercell helpers now return `Result`). Done: the Jacobi solver returns `EigenError` at its cap and converges on a relative test; `pwem_te_bands` returns `Result`.
- [x] `waveguide/fem.rs`: dense, no size guard (vector FEM has one), LU failure only `break`s, no convergence flag, `nx == 0`, stale `sigma` doc. Done: `fem_fundamental_mode` returns `Result` (`FemError`): invalid or oversized domains are refused before allocation (limit 1500 nodes), a singular shift is an error, and the final residual must be below 1e-6. Open: the `sigma` doc still reads loosely.
- [~] `waveguide`: `fiber_lp_modes` fixed 400-point scan misses modes at V ≳ 1000; `graded_lp_modes` fixed R = 4a and `l*l` overflow; `slab_modes` symmetric only and NaN `V` slips through; `pwem_te_bands` even `n_g` / `n_k == 0`. Done: `graded_lp_modes` `l*l` overflow (widened to f64); `fiber_lp_modes` scan density scales with V (the fixed 400 points missed about a sixth of the LP0m modes at V ≈ 1500). Open: the `graded_lp_modes` truncation, the `slab_modes` symmetry, and the other items.
- [ ] `polarization/ray_trace.rs:53`: `PolarizedSystem` holds exactly two elements; document reflection frame/handedness; `rotate_to_frame` (passive) vs `JonesMatrix::rotated` (active).
- [~] `polarization/mueller.rs`: clamp `partial_depolarizer`, validate `diattenuator`; `stokes.rs:61` DOP clamp vs doc; rewrite `crystals.rs` `extraordinary_walkoff` doc. Done: `partial_depolarizer` and `diattenuator` assert their physical range; the DOP clamp is [0, 1]. Open: the `crystals.rs` doc.

### Phase 3 — missing features
- [~] Cooke triplet and double Gauss presets from a verified published prescription (closes the open double-Gauss E2E). Done: `presets::cooke_triplet_50mm` from the published Optiland Cooke-triplet prescription (final optimised values, 50.002 mm), tested in `tests/cooke_triplet_reference.rs` against the tutorial focal length (50.0014 mm here). The glasses are registered at their datasheet d-line indices, so the test is d-line only. Open: double Gauss. The only prescription I found (US 2,784,643) is an OCR-garbled table, so it is not used.
- [ ] Seidel (third-order) aberration sums; OPD-vs-pupil; Strehl; encircled energy; polychromatic spot/MTF; footprint diagram.
- [~] Named tolerance types (tilt, decenter, irregularity, index, Abbe) with a sensitivity table and truncated distributions. Done: `Tolerance::TruncatedNormal` (a sigma inside a hard limit), with a validated limit (`ToleranceError::InvalidTolerance`). Open: the named types (tilt, decenter, irregularity, index, Abbe) and the sensitivity table.
- [~] Thin-film TMM (oblique s/p, complex index) and coating-aware tracing with transmitted power; absorption (`k`) in `Medium`. Done: `thin_film::film_response` for s and p at oblique incidence with complex layer and substrate indices (matches Fresnel for a bare interface, Brewster, AR and power conservation tested). Open: coating-aware tracing and `k` in `Medium`.
- [~] Ring-resonator through/drop transfer function with loss, critical coupling and finesse; coupled-mode phase mismatch Δβ; asymmetric slab; waveguide dispersion engineering. Done: all-pass and add-drop through/drop transmissions (energy-conserving when lossless), critical-coupling condition, finesse. Open: coupled-mode phase mismatch Δβ, asymmetric slab, dispersion engineering.
- [ ] Lu–Chipman decomposition, Mueller realisability check, Jones from Mueller, eigenpolarizations, Jones pupils.
- [ ] Glass substitution / hammer optimisation; Schott temperature coefficients (D0/D1/D2/E0/E1/λtk); vacuum-to-air index correction (~3e-4).
- [ ] Performance: `parallel` feature in geometric/design (Jacobian columns, per-field analysis, glass-map scan, NS rays with per-chunk seeds); cache `zernike_eval` factorials/coefficients; stop rebuilding `TunedParams` per `trace`; `vertex_z` is O(i) per surface.
- [ ] Later: FDTD/RCWA, vector diffraction (Rayleigh–Sommerfeld, Debye–Wolf), Gerchberg–Saxton phase retrieval, Shack–Hartmann.

### Phase 4 — adoption and distribution
- [ ] Make it installable: publish the 9 library crates (root README says `version = "0.1"`, which fails until then); allow the CLI to publish for `cargo install`; cargo-dist prebuilt binaries; maturin wheels (PyPI); `wasm-pack` to npm.
- [ ] README front page: CLI quick start, worked example with real output, add `tpt-optics-cli`/`-wasm`/`-py` to the crate table, link GUIDE/COOKBOOK/examples, badges. Fix the "exactly the same gate" claim.
- [ ] SVG output from the CLI (`draw`/`plot`: layout, spot, ray fans, MTF) and a README gallery; example 03 is named "spot diagram" but prints text.
- [ ] `tpt-optics init --preset singlet|doublet|cooke|double-gauss`, `list-glasses`, `glass --range` table, `cargo-generate` template.
- [ ] CLI: move to clap (help, `--version`, completions, man page); `--format json`; TOML errors with spans and "did you mean"; `assert_cmd` integration tests; `watch` mode. Currently `--flag value` leaves `value` as a positional and typos are accepted.
- [ ] CLI `info` output is stale (omits wasm/py; command list omits `schema`).
- [ ] Python: example notebooks; run `maturin develop` + `pytest` in CI. WASM: run `wasm-pack` + `node tests/smoke.mjs` in CI; deploy `www/` to GitHub Pages (page still not opened in a browser).
- [ ] CI: loop over `crates/tpt-optics/examples/*.rs` instead of the hand-kept 10-line list; `cargo doc -D warnings`; `cargo semver-checks`; `cargo audit`/deny `advisories` (CI checks `licenses sources` only); lychee link check; coverage; criterion regression tracking; dependabot; release-plz; issue/PR templates and CODEOWNERS.
- [ ] Prelude: re-export optimiser, merit-function, tolerance and io types so examples need no deep paths.
- [ ] Empty `crates/tpt-optics-design/tests/`; `verify_system` tests only 3 rays and sequential systems (add TIR and mirror cases); remove or generate hard-coded test counts (`PUBLISHING.md:9` "459 tests", status banner "208 tests").

### Phase 5 — innovation
- [ ] Use the `Dual` path for the LM Jacobian by default in `optimize` (replace finite differences) with glass/asphere/Zernike variables.
- [ ] MCP server: add trace, analyze, optimize and schema tools so an agent can author and iterate `design.toml` (builds on `tpt-optics mcp`).
- [ ] Self-verifying results: report `self_check` + verify invariants (and Kani status) in every CLI result.
- [ ] Browser lens demo on GitHub Pages with live spot diagram and editable prescription.
- [ ] Reproducibility: content-hash design file + library version into every result JSON; include sweep seeds in `manifest.json`.

### Housekeeping
- [ ] Re-tick as `[~]`/`[-]` the items marked `[x]` but deferred/rejected/out of scope (non-sequential `.zmx` writer/parser, binary `.zmx`, STEP, IGES, PCF, Ohara/Hoya/CDGM "via AGF" with no bundled data).
- [ ] Phase 4 exit criterion is ticked while the double-Gauss E2E is open; the "Confirm GitHub org"/"tpt-fem crate name" items are ticked but the migration is not done.
- [ ] If the older "Review follow-up (2026-10-31)" section references a personal plan path, remove it and fix its "nothing was built or run" wording.

# Validation

This page lists the reference checks that the test suite makes, with the
reference value, the tolerance, and the test that makes the check. It is kept
by hand. Each row names a test, so you can run the check with
`cargo test -p <crate> <test name>`, and a missing test will show up as a
broken link to the source.

Use it to judge how far a result can be trusted. The tolerances are the
tolerances the tests assert, not measurement precision.

## Materials and dispersion

| Quantity | Reference | Tolerance | Test |
|---|---|---|---|
| N-BK7 n_d (587.6 nm) | 1.5168 | 2e-4 | `tpt-optics-core` `dispersion` (BK7 Sellmeier) |
| N-BK7 n_F, n_C | 1.52238, 1.51432 (Schott) | 2e-4 | `bk7_single_line_cd_matches_published_indices` (`tpt-optics-io`) |
| N-BK7 V_d | 64.17 | 0.1 | `bk7_abbe_number_and_glass_code_match_the_datasheet` (`tpt-optics-core`) |
| N-BK7 glass code | 517642 | last digit ±1 | same test |
| Fused silica n (1550 nm) | 1.44402 | 5e-5 | `dispersion` (fused silica) |
| Edlén air n (589 nm) | 1.000277 | 1e-6 | `material` (Edlén) |
| Edlén outside its range | NaN | exact | `edlen_outside_its_range_is_nan` |

## Geometric optics

| Quantity | Reference | Tolerance | Test |
|---|---|---|---|
| Singlet EFL | Thick-lens formula | 1e-9 relative | `singlet_efl_matches_thick_lens_formula` (`tpt-optics-geometric`) |
| Plano-convex EFL | R/(n − 1) | 1e-4 relative | `plano_convex_efl_is_r_over_n_minus_one` |
| Achromat F and C foci | Equal at the d-line focal length | 0.02 mm | `achromat_has_the_requested_focal_length_and_equal_f_and_c_foci` |
| Snell reversibility | Round trip returns the ray | proptest | `snell_is_reversible` (`tpt-optics-geometric/tests/properties.rs`) |
| ABCD symplecticity | det = index ratio | proptest | `unit_index_systems_are_symplectic` and `mirror_abcd_is_symplectic` (`tpt-optics-geometric/tests/properties.rs`) |
| Field-ray pupil point | Crosses the first plane at (x, 0) | 1e-12 | `field_ray_crosses_first_plane_at_pupil_point` |
| Focal shift vs curvature | Finite difference | 1e-4 relative | `focal_shift_vs_curvature_matches_finite_difference` |
| Exact gradient (dual numbers) | Central differences | 1e-4 relative | `ray_height_gradient_is_exact_against_central_differences` (`tpt-optics-design`) |
| Analysis focus at the paraxial BFL | Through-focus search | 0.1 mm resolution | `report_matches_the_individual_measurements` (`tpt-optics-design`) |

## Wave optics

| Quantity | Reference | Tolerance | Test |
|---|---|---|---|
| Gaussian beam through a lens | q-parameter ABCD law | suite tolerance | `gaussian_through_lens_matches_q_parameter_abcd_law` |
| Airy disk (uniform pupil) | Analytic | suite tolerance | `vcz_uniform_disk_matches_analytic_airy` |
| Circular MTF | Analytic reference values | suite tolerance | `analytic_circular_mtf_reference_values` |
| Free-space BPM | Analytic evolution | suite tolerance | `free_space_evolution_matches_analytic` |

## Polarisation

| Quantity | Reference | Tolerance | Test |
|---|---|---|---|
| Rotated polariser | Constructor | exact | `rotated_polarizer_matches_constructor` |
| Jones ↔ Mueller consistency | Matrix identities | suite tolerance | `polarization` unit tests |

## Waveguides and photonic crystals

| Quantity | Reference | Tolerance | Test |
|---|---|---|---|
| Slab TE0 n_eff | Transcendental equation | suite tolerance | `fem_slab_mode_matches_transcendental` (`tpt-optics-waveguide`) |
| LP11 cutoff | V = 2.4048 | suite tolerance | `lp11_cutoff_near_2p405` |
| Graded-index ground state | Parabolic oscillator | suite tolerance | `parabolic_ground_state_matches_oscillator` |
| Vector FEM PEC TE10 cutoff | Analytic | 2e-4 relative (16 × 8 cells) | `rectangular_pec_guide_te10_cutoff` |
| Vector FEM dielectric slab TE0 | Transcendental | 1e-4 relative (4 × 64 cells) | `dielectric_slab_te0_matches_exact` |
| Quarter-wave stack stopband | Bragg condition | suite tolerance | `quarter_wave_stack_stopband_brackets_center` |
| Empty lattice folding | Folded light line | suite tolerance | `empty_lattice_matches_folded_light_line` |
| Square-rod TM band gap | (0.324, 0.447) c/a, ε = 8.9, r/a = 0.2 | ±0.02 c/a (grid dependent) | `bulk_gap_is_found_for_square_rods` |
| Missing-rod defect level | 0.393 c/a at Γ | 0.01 c/a | `missing_rod_level_sits_in_tm_gap` |
| Line-defect guided level | Present in the gap | existence | `line_defect_supercell_solves_and_has_a_level_in_the_gap` |
| SMF-28 modes | LP11 cut off at 1550 nm, guided at 1200 nm | existence | `smf28_guides_only_the_fundamental_mode_at_1550` |

## Coupled-mode theory and design

| Quantity | Reference | Tolerance | Test |
|---|---|---|---|
| Directional coupler | Full transfer at L_c | 1e-12 | `coupler_conserves_and_transfers` |
| Grating peak reflectivity | tanh²(κL) | 1e-12 | `spectrum_matches_peak_reflectivity_at_phase_match` |
| Ring FSR | λ²/(n_g L) | ±20 % | `ring_fsr_and_resonances` |
| Tolerance budget | Monte Carlo (400 samples) | 15 % | `linear_budget_agrees_with_monte_carlo_for_small_errors` |
| Parallel tolerance run | Serial, same seeds | exact | `parallel_analysis_matches_the_serial_result_exactly` (feature `parallel`) |
| Sobol first block | 1-D stratification | exact | `sobol_first_block_is_stratified_in_every_dimension` |

## Known gaps

- The FEM initial shift, the photonic-crystal fibre, and the double-Gauss
  reference file are open items in `todo.md`. No check covers them yet.
- The Sobol direction numbers beyond the first pair are checked by
  stratification only, not by the full net property.
- The Cooke triplet and double-Gauss presets are not in the library, so they
  have no reference check.

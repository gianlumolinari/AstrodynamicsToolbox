# Astrodynamics Toolbox (MATLAB)

A modular MATLAB toolbox for **classical astrodynamics**, **trajectory design**, and **dynamical systems analysis in the Circular Restricted Three-Body Problem (CR3BP)**.

The project is under active development and is intended as a research-oriented environment for:

- orbit propagation and perturbation modeling
- Lambert and patched-conics trajectory design
- SPICE-based ephemerides and mission reconstruction
- multiple-gravity-assist analysis
- small-body target screening
- periodic orbit computation in the CR3BP
- family continuation, stability analysis, invariant manifolds, and Poincaré sections

---

## Current Scope

The toolbox currently develops along two main directions:

### 1. Classical astrodynamics
Classical trajectory-analysis capabilities include:

- Cartesian / orbital-element conversions
- two-body propagation
- Lambert solvers
- Hohmann and patched-conics transfers
- J2, drag, and third-body perturbations
- topocentric visibility and ground-station access
- eclipse and conical shadow checks
- multiple-gravity-assist workflows
- SPICE-based mission reconstruction and validation
- small-body / NEO screening examples

### 2. CR3BP and dynamical systems
CR3BP capabilities include:

- equations of motion and zero-velocity curves
- Lagrange point computation
- JPL periodic-orbit seed retrieval
- planar Lyapunov and halo differential correction
- single-shooting and multiple-shooting correction
- natural and pseudo-arclength continuation
- family generation and visualization
- monodromy matrix and eigenspectrum analysis
- stability diagnostics
- stable and unstable manifold generation
- Poincaré sections
- validation against literature benchmarks such as Howell (1984)

---

## Repository Structure

```text
AstrodynamicsToolbox/
├── +astro/                 % Main MATLAB package
├── data/                   % SPICE kernels, constants, auxiliary data
├── examples/
│   ├── classical/
│   │   ├── testing/
│   │   ├── validation/
│   │   └── demos/
│   └── cr3bp/
│       ├── testing/
│       ├── validation/
│       └── demos/
├── external/               % External dependencies (e.g. MICE/SPICE)
├── tests/                  % Automated / regression tests (future expansion)
├── startup.m               % Adds toolbox folders to MATLAB path
└── README.md


## Example Philosophy

Examples are grouped by both **domain** and **purpose**.

### `testing`

These are development-oriented scripts used to:

- verify algorithms
- debug new implementations
- inspect intermediate numerical behavior
- test convergence of correctors and propagators

Typical examples include:

- Lambert solver checks
- perturbation-model tests
- CR3BP differential-correction tests
- JPL seed-retrieval checks

### `validation`

These scripts compare toolbox results against:

- published literature
- tabulated benchmark solutions
- SPICE truth trajectories
- reference mission data

Typical examples include:

- Galileo / JUICE trajectory validation
- JUICE SPICE truth comparison
- Howell (1984) halo-family validation
- CR3BP eigenspectrum benchmark checks

### `demos`

These are more polished examples intended to showcase toolbox capabilities end-to-end.

Typical examples include:

- Earth–Mars mission design workflow
- small-body target-screening demo
- MGA showcase examples
- halo-family continuation
- manifold and Poincaré-section visualizations

---

## Current Example Layout

### Classical examples

#### Testing

Located in `examples/classical/testing`

These cover low-level and intermediate workflows such as:

- orbital element conversion
- two-body propagation
- Hohmann transfer checks
- Lambert validation
- eclipse and visibility utilities
- perturbation modeling
- MGA feasibility and optimization experiments

Current classical testing examples include:

- `test_coe_rv_conversion.m`
- `test_two_body_propagation.m`
- `test_hohmann_transfer.m`
- `test_basic_lambert.m`
- `test_lambert_propagation_check.m`
- `test_compare_universal_vs_izzo.m`
- `test_earth_mars_horizons.m`
- `test_earth_eclipse_check.m`
- `test_ground_station_access.m`
- `test_conical_eclipse_check.m`
- `test_j2_drag_propagation.m`
- `test_topocentric_visibility.m`
- `test_third_body_perturbation.m`
- `test_mga_feasibility_check.m`
- `test_plot_single_flyby.m`
- `test_optimize_mga_dates.m`
- `test_optimize_mga_dates_constrained.m`
- `test_scan_mga_sequence.m`
- `test_optimize_mga_departure_and_tofs.m`

#### Validation

Located in `examples/classical/validation`

These currently include benchmark-style scripts such as:

- Galileo VEEGA reference reconstruction
- JUICE reference trajectory reconstruction
- JUICE SPICE truth comparison

Current classical validation examples include:

- `validate_galileo_veega_reference.m`
- `validate_juice_reference_trajectory.m`
- `validate_juice_spice_truth_comparison.m`

#### Demos

Located in `examples/classical/demos`

These include more polished showcase workflows such as:

- NEO screening
- end-to-end Earth–Mars mission design
- Earth–Venus–Earth–Jupiter MGA example
- Earth–Mars patched-conics analysis
- Earth–Mars window scan
- constrained planetary mission design under `\Delta v` limits

Current classical demo examples include:

- `demo_screen_nea_targets.m`
- `demo_end_to_end_earth_mars_mission_design.m`
- `demo_earth_venus_earth_jupiter.m`
- `demo_planetary_mission_minTOF_under_dv_constraint.m`
- `demo_earth_mars_window_scan.m`
- `demo_earth_mars_patched_conics.m`

---

### CR3BP examples

#### Testing

Located in `examples/cr3bp/testing`

These cover development steps such as:

- Lagrange points and zero-velocity curves
- JPL periodic-orbit seed retrieval
- family plotting
- planar Lyapunov correction
- single-shooting and multiple-shooting periodic correction
- natural parameter continuation

Current CR3BP testing examples include:

- `test_cr3bp_lagrange_and_zvc.m`
- `test_jpl_periodic_orbit_seed.m`
- `test_jpl_full_family_plot.m`
- `test_planar_lyapunov_differential_correction.m`
- `test_single_shooting_periodic.m`
- `test_multiple_shooting_periodic.m`
- `test_natural_parameter_continuation.m`

#### Validation

Located in `examples/cr3bp/validation`

These include benchmark-oriented scripts such as:

- halo eigenspectrum inspection
- Howell (1984) halo-family validation
- Howell-style manifold-surface validation

Current CR3BP validation examples include:

- `validate_halo_eigenspectrum_inspection.m`
- `validate_howell1984_halo_family.m`
- `validate_howell1984_manifold_surface.m`

#### Demos

Located in `examples/cr3bp/demos`

These showcase higher-level workflows such as:

- pseudo-arclength continuation
- family stability analysis
- manifold generation
- Poincaré section analysis
- halo differential correction
- halo family continuation
- halo manifold-tube visualization

Current CR3BP demo examples include:

- `demo_pseudo_arclength_continuation.m`
- `demo_family_stability_analysis.m`
- `demo_manifold_generation.m`
- `demo_poincare_section_manifolds.m`
- `demo_halo_differential_correction.m`
- `demo_halo_family_continuation.m`
- `demo_halo_family_stability_analysis.m`
- `demo_halo_manifold_tubes.m`

---

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/gianlumolinari/AstrodynamicsToolbox.git
cd AstrodynamicsToolbox

Getting Started

1. Clone the repository

```bash
git clone https://github.com/gianlumolinari/AstrodynamicsToolbox.git
cd AstrodynamicsToolbox

2. Open MATLAB in the repository root
Then run:

startup

This adds the toolbox and example folders to the MATLAB path.

Suggested Starting Points

If you are new to the toolbox, these are good entry points.

Classical
	•	examples/classical/testing/test_two_body_propagation.m
	•	examples/classical/testing/test_basic_lambert.m
	•	examples/classical/validation/validate_juice_spice_truth_comparison.m
	•	examples/classical/demos/demo_end_to_end_earth_mars_mission_design.m

CR3BP
	•	examples/cr3bp/testing/test_cr3bp_lagrange_and_zvc.m
	•	examples/cr3bp/testing/test_planar_lyapunov_differential_correction.m
	•	examples/cr3bp/validation/validate_howell1984_halo_family.m
	•	examples/cr3bp/demos/demo_halo_family_continuation.m
	•	examples/cr3bp/demos/demo_poincare_section_manifolds.m

⸻

Validation Strategy

A core design goal of the toolbox is benchmark consistency.

Validation workflows currently include comparison against:
	•	SPICE truth trajectories for mission reconstruction examples
	•	JPL periodic-orbit seeds for CR3BP family studies
	•	literature benchmarks such as Howell (1984) for halo orbits

This makes the toolbox useful not only for development and exploration, but also for reproducible astrodynamics studies.

⸻

Development Status

This is an ongoing project.

The current focus is on building a robust and extensible foundation for:
	•	classical mission design workflows
	•	CR3BP periodic orbit computation and continuation
	•	invariant manifold analysis
	•	benchmark validation against literature and ephemerides

Planned future directions include:
	•	more automated regression tests
	•	additional optimization workflows
	•	improved manifold-surface visualization
	•	expanded small-body targeting tools
	•	higher-fidelity models beyond the CRTBP where appropriate

⸻

Design Principles

The toolbox is being developed with the following priorities:
	•	modularity — reusable functions under the +astro package
	•	transparency — examples written as readable workflows rather than black boxes
	•	validation — benchmark scripts included alongside development and demo cases
	•	research utility — suitable for both learning and more advanced exploratory work

⸻

Notes
	•	Most examples are written as transparent scripts to make the workflow easy to inspect and adapt.
	•	Validation scripts are especially useful as reference cases when extending the toolbox.
	•	The CR3BP branch is being developed with a strong emphasis on literature consistency and benchmark reproducibility.
	•	Some visualizations, especially invariant-manifold plots, may continue to evolve as benchmark cases and plotting conventions are refined.

⸻

Author

Gianluca Molinari
Astrodynamics / trajectory design / CR3BP / mission analysis


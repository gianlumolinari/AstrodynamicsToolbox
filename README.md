# Astrodynamics Toolbox (MATLAB)

A modular MATLAB toolbox for **classical astrodynamics**, **trajectory design**, **CR3BP**, and **dynamical systems analysis**.

Designed to be highly transparent and strictly validated, this toolbox provides a robust foundation for everything from two-body propagation and Lambert / patched-conics design to advanced periodic orbit computation, continuation methods, and SPICE-based ephemeris validation. It is built for both educational learning and advanced research.

---

## 🚀 Getting Started

### 1. Installation
Clone the repository to your local machine:
```bash
git clone https://github.com/gianlumolinari/AstrodynamicsToolbox.git
cd AstrodynamicsToolbox
```

### 2. Initialization
Open MATLAB in the repository root and run the startup script. This will automatically add the toolbox and necessary data folders to your MATLAB path.
```matlab
startup
```

### 3. Explore Interactive Live Scripts
If you are new to the toolbox, the best way to start is by opening the interactive MATLAB Live Scripts (`.mlx`) located in the `liveScripts/` directory:
* `classical_astrodynamics.mlx`: Orbit geometry, propagation, perturbations, Lambert transfers, and MGA workflows.
* `cr3bp_dynamical_systems.mlx`: Zero-velocity curves, periodic orbits, differential correction, continuation, and manifolds.

---

## 🛰️ Core Capabilities

The toolbox is actively developed along two main branches.

### 1. Classical Astrodynamics
* **Core Mechanics:** Two-body propagation, Cartesian / orbital-element conversions, Lambert solvers (Izzo and universal).
* **Mission Design:** Hohmann and patched-conics transfers, multiple-gravity-assist (MGA) workflows, small-body / NEO screening.
* **Perturbations:** J2, atmospheric drag, and third-body perturbation modeling.
* **Geometry & Access:** Topocentric visibility, ground-station access, eclipse, and conical shadow checks.
* **Validation:** SPICE-based ephemerides retrieval, mission reconstruction, and validation.

### 2. CR3BP & Dynamical Systems
* **Fundamentals:** Equations of motion, zero-velocity curves, and Lagrange point computation.
* **Periodic Orbits:** JPL periodic-orbit seed retrieval, single / multiple-shooting correction, planar Lyapunov, and halo differential correction.
* **Continuation & Stability:** Natural and pseudo-arclength continuation, monodromy matrix analysis, eigenspectra, and stability diagnostics.
* **Manifolds:** Stable / unstable manifold generation and Poincaré sections.
* **Validation:** Benchmarked against literature and frozen JPL catalog data.

---

## 💻 Quick Start Examples

### Classical Propagation Example
```matlab
startup

mu = 398600.4418; % Earth gravitational parameter
[r0, v0] = astro.coords.coe2rv(12000, 0.25, 35, 40, 60, 20, mu);
T = astro.maneuvers.orbitalPeriod(12000, mu);
out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, mu), ...
    [0 T], [r0; v0]);

figure
plot3(out.x(:,1), out.x(:,2), out.x(:,3), 'LineWidth', 1.8)
hold on
plot3(out.x(1,1), out.x(1,2), out.x(1,3), 'o', 'MarkerSize', 8, 'LineWidth', 1.5)
plot3(0, 0, 0, 'o', 'MarkerSize', 10, 'LineWidth', 1.5)
grid on
axis equal
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('Two-Body Orbit Propagation')
legend('Trajectory', 'Initial state', 'Earth', 'Location', 'best')
view(3)
```

### CR3BP Propagation Example
```matlab
startup

mu = 0.012150585609624;  % Earth-Moon mass ratio

x0 = [-0.84094370;
       0.00000000;
       0.00000000;
       0.00000000;
       0.25043310;
       0.00000000];

T = 6.27459754;

out = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0 T], x0);

figure
plot3(out.x(:,1), out.x(:,2), out.x(:,3), 'LineWidth', 1.5)
hold on
plot3(-mu, 0, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
plot3(1-mu, 0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
grid on
axis equal
xlabel('x [-]')
ylabel('y [-]')
zlabel('z [-]')
title('Earth-Moon 3:1 Resonant Periodic Orbit')
legend('Trajectory', 'Earth', 'Moon', 'Location', 'best')
view(3)
```

---

## 📁 Repository Structure

The core functionality lives in the `+astro/` package, ensuring modularity and avoiding namespace clashes.

```text
AstrodynamicsToolbox/
├── +astro/                     % Main MATLAB package containing the toolbox source code
│   ├── +bodies/                % Celestial body definitions, physical parameters, and constants
│   ├── +constants/             % General constants used across the toolbox
│   ├── +continuation/          % General continuation algorithms and family-tracking utilities
│   ├── +coords/                % Coordinate, frame, and orbital-element conversions
│   ├── +cr3bp/                 % CR3BP dynamics, differential correction, continuation, and manifolds
│   ├── +ephem/                 % SPICE/Horizons interfaces, kernel loading, and ephemeris retrieval
│   ├── +geometry/              % Geometric tools such as occultation, eclipse, and line-of-sight checks
│   ├── +lambert/               % Lambert problem solvers
│   ├── +maneuvers/             % Classical maneuver analysis and mission-design utilities
│   ├── +manifolds/             % General manifold-related utilities
│   ├── +mga/                   % Multiple-gravity-assist trajectory utilities
│   ├── +opt/                   % Optimization and design support tools
│   ├── +periodic/              % Periodic-orbit helpers and related utilities
│   ├── +perturbations/         % J2, drag, SRP, and third-body models
│   ├── +plot/                  % Plotting and visualization utilities
│   ├── +propagators/           % Numerical propagation routines
│   ├── +smallbody/             % Small-body and asteroid-related utilities
│   ├── +units/                 % Unit-conversion helpers
│   ├── +utils/                 % General-purpose utilities
│   └── +visibility/            % Visibility, access, and observation geometry tools
├── data/                       % SPICE kernels and auxiliary mission/ephemeris data
├── docs/                       % Documentation, notes, derivations, and longer-form writeups
├── examples/                   % Example scripts organised by domain and purpose
├── external/                   % External dependencies, including NAIF MICE/SPICE
├── liveScripts/                % Interactive MATLAB Live Scripts
├── tests/                      % Automated test and validation framework
├── startup.m                   % Initialization script
└── README.md                   
```

### CR3BP Package
The `+astro/+cr3bp` package contains the main dynamical-systems functionality.

```text
+astro/+cr3bp/
├── eomCR3BP.m
├── variationalEOM.m
├── effectivePotential.m
├── jacobiConstant.m
├── lagrangePoints.m
│
├── propagateWithSTM.m
├── monodromyMatrix.m
├── stabilityIndices.m
│
├── differentialCorrectionPlanarLyapunov.m
├── differentialCorrectionHalo.m
├── singleShootingCorrector.m
├── multipleShootingCorrector.m
│
├── continueFamilyNatural.m
├── continueFamilyPseudoArc.m
├── continueHaloPseudoArc.m
├── correctPlanarLyapunovPseudoArc.m
├── correctHaloPseudoArc.m
│
├── manifoldSeeds.m
├── haloManifoldSeeds.m
├── haloManifoldBranch.m
├── howellManifoldBranch.m
├── propagateManifold.m
│
├── propagateToSection.m
├── collectPoincareSection.m
│
├── queryJPLPeriodicOrbits.m
└── parseJPLPeriodicOrbit.m
```

### Examples and Live Scripts
The examples and live scripts are intended to provide interactive, visual entry points into the toolbox.

```text
examples/
├── classical/
│   ├── testing/
│   ├── validation/
│   └── demos/
└── cr3bp/
    ├── testing/
    ├── validation/
    └── demos/

liveScripts/
├── classical_astrodynamics.mlx
└── cr3bp_dynamical_systems.mlx
```

---

## ✅ Validation Strategy

A core design goal of the toolbox is **benchmark consistency and transparency**. 
Most workflows are implemented as readable scripts rather than opaque black boxes, making it easier to verify assumptions, inspect intermediate quantities, and compare against references.

The `examples/` directory is organised into three categories:
* **`testing/`** — development and debugging scripts used to verify implementations and convergence behaviour
* **`validation/`** — comparisons against literature, tabulated benchmarks, or SPICE truth trajectories
* **`demos/`** — polished end-to-end workflows showcasing the toolbox capabilities

---

## 🧪 Test Structure

This folder contains the automated test and validation framework for the Astrodynamics Toolbox.

### Layout

```text
tests/
├── README.md
├── integration/
│   ├── cr3bp/
│   │   ├── testContinuationPseudoArcEarthMoon.m
│   │   ├── testLyapunovFamilyTrendEarthMoon.m
│   │   └── testMultipleShootingCorrectorEarthMoon.m
│   └── lambert/
│       └── test_earth_mars_horizons.m
├── unit/
│   ├── classical/
│   │   ├── test_coe_rv_roundtrip.m
│   │   └── test_two_body_energy.m
│   ├── cr3bp/
│   │   ├── testJplHaloClosureEarthMoon.m
│   │   ├── testJplJacobiConsistency.m
│   │   ├── testJplJacobiDriftEarthMoon.m
│   │   ├── testMonodromyMatrixEarthMoon.m
│   │   └── testSingleShootingCorrectorEarthMoon.m
│   └── lambert/
│       ├── test_basic_lambert.m
│       ├── test_izzo_vs_universal.m
│       └── test_lambert_propagation_check.m
└── validation/
    ├── runCR3BPValidationSuite.m
    └── testJplCatalogData.m
```

### Categories

#### Unit Tests
Unit tests validate individual functions, invariants, or specific numerical properties. They are intended to be fast, deterministic, and easy to diagnose.

**`tests/unit/classical/`**
* `test_coe_rv_roundtrip.m`: verifies consistent conversion between orbital elements and Cartesian state vectors.
* `test_two_body_energy.m`: checks conservation of two-body orbital energy under propagation.

**`tests/unit/lambert/`**
* `test_basic_lambert.m`: basic correctness test for Lambert solutions on a simple transfer case.
* `test_izzo_vs_universal.m`: cross-validates the Izzo Lambert solver against the universal-variable implementation.
* `test_lambert_propagation_check.m`: verifies that a Lambert solution reproduces the expected transfer when propagated.

**`tests/unit/cr3bp/`**
* `testJplJacobiConsistency.m`: confirms that the Jacobi constant recomputed from frozen JPL catalog states matches the catalog value.
* `testJplJacobiDriftEarthMoon.m`: verifies Jacobi conservation during Earth-Moon CR3BP propagation.
* `testJplHaloClosureEarthMoon.m`: propagates Earth-Moon JPL halo states over one period and checks periodic closure.
* `testMonodromyMatrixEarthMoon.m`: validates monodromy matrix properties, including determinant consistency and reciprocal eigenvalue pairing.
* `testSingleShootingCorrectorEarthMoon.m`: validates the generic single-shooting corrector on a robust Earth-Moon L1 Lyapunov case.

#### Integration Tests
Integration tests validate short multi-step workflows involving several functions together.

**`tests/integration/lambert/`**
* `test_earth_mars_horizons.m`: tests an Earth-Mars transfer workflow against ephemeris-based reference conditions.

**`tests/integration/cr3bp/`**
* `testMultipleShootingCorrectorEarthMoon.m`: validates multiple-shooting periodic correction on Earth-Moon halo cases.
* `testContinuationPseudoArcEarthMoon.m`: validates pseudo-arclength continuation for Earth-Moon Lyapunov families.
* `testLyapunovFamilyTrendEarthMoon.m`: compares locally continued Lyapunov-family trends against JPL reference data.

#### Validation Tests
Validation tests focus on external reference datasets and higher-level suite runners.

**`tests/validation/`**
* `testJplCatalogData.m`: checks that JPL periodic-orbit benchmark files exist and contain the expected schema.
* `runCR3BPValidationSuite.m`: runs the main CR3BP validation suite in one command.

### Running the CR3BP Validation Suite
From MATLAB, in the repository root:
```matlab
startup
results = runCR3BPValidationSuite;
```

### Benchmark Data
Several CR3BP tests rely on frozen benchmark files derived from the JPL periodic orbit catalog and stored under:
```text
data/validation/cr3bp/
```
These datasets are used to keep validation deterministic and independent of live API access during testing.

---

## 🗺️ Roadmap & Philosophy

This is an actively evolving project built around **modularity**, **transparency**, and **rigorous validation**.

Planned directions include:
* Automated regression and unit testing
* Expanded optimisation and targeting workflows
* Higher-fidelity perturbation models
* Improved manifold and family visualisation tools
* Additional mission design case studies and benchmark reproductions
* Covariance analysis and propagation

---

## Author

**Gianluca Molinari** *Astrodynamics | Trajectory Design | CR3BP | Mission Analysis*
` ` `
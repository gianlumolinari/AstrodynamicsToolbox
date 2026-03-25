# Astrodynamics Toolbox (MATLAB)
![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg)

A modular MATLAB toolbox for **classical astrodynamics**, **trajectory design**, **CR3BP** and **dynamical systems** analysis. 

Designed to be highly transparent and strictly validated, this toolbox provides a robust foundation for everything from two-body propagation and Lambert/patched-conics design to advanced periodic orbit computation, continuation methods, and SPICE-based ephemeris validation. It is built for both educational learning and advanced research.

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

The toolbox is actively developed along two main branches:

### 1. Classical Astrodynamics
* **Core Mechanics:** Two-body propagation, Cartesian/orbital-element conversions, Lambert solvers (Izzo and universal).
* **Mission Design:** Hohmann and patched-conics transfers, multiple-gravity-assist (MGA) workflows, small-body/NEO screening.
* **Perturbations:** J2, atmospheric drag, and third-body perturbation modeling.
* **Geometry & Access:** Topocentric visibility, ground-station access, eclipse, and conical shadow checks.
* **Validation:** SPICE-based ephemerides retrieval, mission reconstruction, and validation.

### 2. CR3BP & Dynamical Systems
* **Fundamentals:** Equations of motion, zero-velocity curves, and Lagrange point computation.
* **Periodic Orbits:** JPL periodic-orbit seed retrieval, single/multiple-shooting correction, planar Lyapunov, and halo differential correction.
* **Continuation & Stability:** Natural and pseudo-arclength continuation, monodromy matrix analysis, eigenspectrums, and stability diagnostics.
* **Manifolds:** Stable/unstable manifold generation and Poincaré sections.
* **Validation:** Benchmarked against literature (e.g., Howell 1984).

---

## 💻 Quick Start Examples

### Classical Propagation Example
```matlab
startup

mu = 398600.4418; % Earth gravitational parameter
[r0, v0] = astro.coords.coe2rv(12000, 0.25, 35, 40, 60, 20, mu);

out = astro.propagators.propagate( ...
    @(t,x) astro.propagators.eomTwoBody(t, x, mu), ...
    [0 10000], [r0; v0]);

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

mu = 0.012150585609624;   % Earth-Moon mass ratio
L = astro.cr3bp.lagrangePoints(mu);

x0 = [L.L1(1)-0.01; 0; 0; 0; 0.1; 0];
C = astro.cr3bp.jacobiConstant(x0, mu);

out = astro.propagators.propagate( ...
    @(t,x) astro.cr3bp.eomCR3BP(t, x, mu), ...
    [0 5], x0);

figure
plot3(out.x(:,1), out.x(:,2), out.x(:,3), 'LineWidth', 1.8)
hold on
plot3(x0(1), x0(2), x0(3), 'o', 'MarkerSize', 8, 'LineWidth', 1.5)
plot3(-mu, 0, 0, 'o', 'MarkerSize', 10, 'LineWidth', 1.5)
plot3(1-mu, 0, 0, 'o', 'MarkerSize', 8, 'LineWidth', 1.5)
grid on
axis equal
xlabel('x [-]')
ylabel('y [-]')
zlabel('z [-]')
title('CR3BP Trajectory Propagation')
legend('Trajectory', 'Initial state', 'Primary 1', 'Primary 2', 'Location', 'best')
view(3)
```

---

## 📁 Repository Structure

The core functionality lives in the `+astro/` package, ensuring modularity and preventing namespace clashes.

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
│   ├── +lambert/               % Lambert problem solvers (e.g. Izzo and universal-variable methods)
│   ├── +maneuvers/             % Classical maneuver analysis and mission-design utilities
│   ├── +manifolds/             % General manifold-related utilities outside the CR3BP-specific package
│   ├── +mga/                   % Multiple-gravity-assist trajectory utilities
│   ├── +opt/                   % Optimization and design support tools
│   ├── +periodic/              % Periodic-orbit helpers and related utilities
│   ├── +perturbations/         % Perturbation models such as J2, drag, SRP, and third-body effects
│   ├── +plot/                  % Plotting and visualization utilities for trajectories and dynamical structures
│   ├── +propagators/           % Numerical propagation routines, including two-body, Cowell, STM, and high-fidelity models
│   ├── +smallbody/             % Small-body and asteroid-related utilities
│   ├── +units/                 % Unit-conversion helpers
│   ├── +utils/                 % General-purpose utilities (e.g. energy, angular momentum, helper routines)
│   └── +visibility/            % Visibility, access, and observation geometry tools
├── data/                       % SPICE kernels and auxiliary mission/ephemeris data (JUICE, DART, HERA, etc.)
├── docs/                       % Documentation, notes, derivations, and longer-form writeups
├── examples/                   % Example scripts organised by domain and purpose (testing, validation, demos)
├── external/                   % External dependencies, including the NAIF MICE/SPICE MATLAB interface
├── liveScripts/                % Interactive MATLAB Live Scripts for guided exploration of the toolbox
└── startup.m                   % Initialization script that adds the toolbox folders to the MATLAB path

```
### CR3BP package
The `+astro/+cr3bp` package contains the main dynamical systems functionality.

```text
+astro/+cr3bp/
├── eomCR3BP.m                        % CR3BP equations of motion
├── variationalEOM.m                  % CR3BP variational equations for STM propagation
├── effectivePotential.m              % Effective potential Ω(x,y,z)
├── jacobiConstant.m                  % Jacobi integral evaluation
├── lagrangePoints.m                  % L1–L5 computation
│
├── propagateWithSTM.m                % State + STM propagation
├── monodromyMatrix.m                 % Monodromy matrix over one orbit period
├── stabilityIndices.m                % Stability metrics from monodromy eigenstructure
│
├── differentialCorrectionPlanarLyapunov.m   % Planar Lyapunov differential correction
├── differentialCorrectionHalo.m               % Halo differential correction
├── singleShootingCorrector.m                  % Generic single-shooting corrector
├── multipleShootingCorrector.m                % Multiple-shooting corrector
│
├── continueFamilyNatural.m           % Natural-parameter family continuation
├── continueFamilyPseudoArc.m         % Pseudo-arclength continuation
├── continueHaloPseudoArc.m           % Halo-family pseudo-arclength continuation
├── correctPlanarLyapunovPseudoArc.m  % Planar pseudo-arclength corrector step
├── correctHaloPseudoArc.m            % Halo pseudo-arclength corrector step
│
├── manifoldSeeds.m                   % Generic manifold seed generation
├── haloManifoldSeeds.m               % Halo manifold seed rings from stable/unstable subspaces
├── haloManifoldBranch.m              % Ross-style halo manifold branch seeds
├── howellManifoldBranch.m            % Howell-style manifold branch helper
├── propagateManifold.m               % Manifold trajectory propagation
│
├── propagateToSection.m              % Propagate trajectory to a section crossing
├── collectPoincareSection.m          % Collect section crossings from multiple seeds
│
├── queryJPLPeriodicOrbits.m          % Query JPL periodic-orbit database
└── parseJPLPeriodicOrbit.m           % Parse JPL periodic-orbit data into toolbox format
```
### Examples and Live scripts
The examples and live scripts are intended to provide interactive, visual entry points into the toolbox:


```text
examples/
├── classical/
│   ├── testing/                      % Development and debugging scripts
│   ├── validation/                   % SPICE / literature / mission validation workflows
│   └── demos/                        % End-to-end showcase examples
│
└── cr3bp/
    ├── testing/                      % Differential correction, continuation, STM, manifolds
    ├── validation/                   % Benchmark comparisons and reference checks
    └── demos/                        % Halo families, manifold visualisation, Poincaré maps

liveScripts/
├── classical_astrodynamics.mlx       % Interactive introduction to two-body, Lambert, perturbations, MGA
└── cr3bp_dynamical_systems.mlx       % Interactive introduction to CR3BP, periodic orbits, continuation, manifolds
```
## ✅ Validation Strategy

A core design goal of the toolbox is **benchmark consistency and transparency**.  
Most workflows are implemented as readable scripts rather than opaque black boxes, making it easier to verify assumptions, inspect intermediate quantities, and compare against references.

The `examples/` directory is organised into three categories:

1. **`testing/`** — development and debugging scripts used to verify implementations and convergence behaviour  
2. **`validation/`** — comparisons against literature, tabulated benchmarks, or SPICE truth trajectories  
3. **`demos/`** — polished end-to-end workflows showcasing the toolbox capabilities

---

## 🗺️ Roadmap & Philosophy

This is an actively evolving project built around **modularity**, **transparency**, and **rigorous validation**.

Planned directions include:
- automated regression and unit testing
- expanded optimisation and targeting workflows
- higher-fidelity perturbation models
- improved manifold and family visualisation tools
- additional mission design case studies and benchmark reproductions

---

## Author

**Gianluca Molinari**  
Astrodynamics | Trajectory Design | CR3BP | Mission Analysis
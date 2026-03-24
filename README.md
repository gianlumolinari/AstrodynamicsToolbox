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
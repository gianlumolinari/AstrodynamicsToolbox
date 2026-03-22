# astroToolbox

A MATLAB astrodynamics toolbox developed during my PhD for trajectory design, interplanetary mission analysis, Lambert solving, patched-conics studies, and later multi-body dynamics.

## Current capabilities

### Classical astrodynamics
- Planetary body definitions
- COE ↔ Cartesian conversion
- Two-body equations of motion
- Numerical propagation
- Energy and angular momentum diagnostics
- 2D/3D orbit plotting
- Circular velocity, escape velocity
- Hohmann transfer and plane change utilities

### Lambert solvers
- Universal-variable Lambert solver
- Izzo-style Lambert solver
- Propagation-based Lambert validation
- Earth–Mars Lambert examples

### Mission design
- Patched-conics utilities
- Hyperbolic excess, C3, departure and arrival delta-v
- Sphere of influence estimates
- Hyperbolic turning angle and flyby periapsis utilities
- Earth–Mars porkchop plotting
- Constrained planetary mission search (minimum TOF under delta-v limit)

### Ephemerides
- JPL Horizons interface
- SPICE/MICE interface for local planetary ephemerides

## Planned additions
- Total mission-cost porkchops
- Top-N feasible mission reporting
- Transfer geometry visualization for optimal missions
- Flyby / B-plane tools
- J2 and higher-fidelity propagators
- Multi-leg patched-conics mission design
- CR3BP, invariant manifolds, periodic orbits, and continuation tools

## Repository structure

```text
astroToolbox/
│
├── README.md
├── startup.m
├── .gitignore
├── examples/
├── tests/
├── data/
└── +astro/
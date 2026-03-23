# AstrodynamicsToolbox

A modular MATLAB toolbox for astrodynamics, mission design, and dynamical systems analysis.

The project is being developed as a personal and research-oriented toolbox with two main pillars:

1. **Classical astrodynamics**
   - two-body propagation
   - Lambert transfer design
   - patched-conics analysis
   - ephemeris-driven interplanetary trajectories
   - gravity-assist and multiple-gravity-assist workflows
   - perturbation modeling
   - eclipse and visibility analysis

2. **CR3BP and dynamical systems**
   - libration-point analysis
   - retrieval of periodic-orbit seeds from the JPL API
   - periodic-orbit propagation and visualization
   - differential correction and shooting methods
   - family continuation
   - monodromy and stability analysis
   - stable/unstable manifold generation
   - Poincaré sections

The long-term goal is to build a single MATLAB environment that supports both traditional mission design and modern multi-body trajectory design.

---

## Current capabilities

### Classical astrodynamics

- Two-body propagation in Cartesian form
- Lambert solver workflows
- Earth-Mars and interplanetary transfer examples
- Patched-conics mission metrics
  - \(v_\infty\)
  - launch \(C_3\)
  - departure and arrival \(\Delta v\)
- Gravity-assist / MGA prototype workflows
- SPICE-based planetary ephemerides
- JUICE and Galileo-style reference reconstruction examples
- J2 and drag perturbation propagation
- Third-body perturbation examples
- Eclipse analysis
  - cylindrical shadow model
  - conical sunlight / penumbra / umbra classification
- Ground-station visibility and access analysis
- Small-body screening utilities for NEO target selection

### CR3BP and dynamical systems

- CR3BP equations of motion
- Jacobi constant evaluation
- Lagrange point computation
- Zero-velocity / energetics utilities
- JPL periodic-orbit API query support
- Parsing of periodic-orbit family data from JPL
- Family propagation and color-mapped family visualization
- Planar Lyapunov correction experiments
- Single- and multiple-shooting infrastructure
- Natural continuation prototype
- Pseudo-arclength continuation for planar Lyapunov families
- Monodromy matrix computation
- Stability metrics
  - unstable multiplier
  - stable reciprocal multiplier
  - spectral radius
  - stability index
- Stable and unstable manifold seed generation
- Manifold propagation
- Poincaré section generation for manifold crossings

---

## Repository structure

A simplified view of the current layout:

```text
AstrodynamicsToolbox/
│
├── +astro/
│   ├── +bodies/
│   ├── +coords/
│   ├── +cr3bp/
│   ├── +ephem/
│   ├── +geometry/
│   ├── +lambert/
│   ├── +maneuvers/
│   ├── +mga/
│   ├── +propagators/
│   ├── +smallbody/
│   └── +visibility/
│
├── data/
│   └── spice/
│
├── examples/
│   ├── classical/
│   ├── cr3bp/
│   ├── mga/
│   ├── smallbody/
│   └── visibility/
│
├── startup.m
└── README.md
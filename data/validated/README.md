# Validated Orbit Data

This folder stores corrected and validated CR3BP periodic orbits used by examples, manifold generation, and regression tests.

## Structure

- `cr3bp/earth_moon/`
  - validated periodic orbits in normalized Earth-Moon CR3BP coordinates

## Naming convention

`<family>_<libration point>_<qualifier>.mat`

Examples:
- `planar_lyapunov_L1_validated.mat`
- `halo_L1_north_validated.mat`

## Minimum contents of each file

- `state0`
- `T`
- `mu`
- `family`
- `libr`
- `source`
- `closureError`
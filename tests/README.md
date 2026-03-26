Here is the raw text with spaces between the backticks, exactly as requested, so you can easily highlight and copy the whole thing:

# Test Structure

This folder contains the automated test and validation framework for the Astrodynamics Toolbox.

## Layout

``` text
.
├── README.md
├── integration
│   ├── cr3bp
│   │   ├── testContinuationPseudoArcEarthMoon.m
│   │   ├── testLyapunovFamilyTrendEarthMoon.m
│   │   └── testMultipleShootingCorrectorEarthMoon.m
│   └── lambert
│       └── test_earth_mars_horizons.m
├── unit
│   ├── classical
│   │   ├── test_coe_rv_roundtrip.m
│   │   └── test_two_body_energy.m
│   ├── cr3bp
│   │   ├── testJplHaloClosureEarthMoon.m
│   │   ├── testJplJacobiConsistency.m
│   │   ├── testJplJacobiDriftEarthMoon.m
│   │   ├── testMonodromyMatrixEarthMoon.m
│   │   └── testSingleShootingCorrectorEarthMoon.m
│   └── lambert
│       ├── test_basic_lambert.m
│       ├── test_izzo_vs_universal.m
│       └── test_lambert_propagation_check.m
└── validation
    ├── runCR3BPValidationSuite.m
    └── testJplCatalogData.m 
```

## Categories

### Unit Tests
Unit tests validate individual functions, invariants, or tightly scoped numerical properties. They are intended to be fast, deterministic, and easy to diagnose.

#### `tests/unit/classical/`
* **`test_coe_rv_roundtrip.m`**: Verifies consistent conversion between orbital elements and Cartesian state vectors.
* **`test_two_body_energy.m`**: Checks conservation of two-body orbital energy under propagation.

#### `tests/unit/lambert/`
* **`test_basic_lambert.m`**: Basic correctness test for Lambert solutions on a simple transfer case.
* **`test_izzo_vs_universal.m`**: Cross-validates the Izzo Lambert solver against the universal-variable implementation.
* **`test_lambert_propagation_check.m`**: Verifies that a Lambert solution reproduces the expected transfer when propagated.

#### `tests/unit/cr3bp/`
* **`testJplJacobiConsistency.m`**: Confirms that the Jacobi constant recomputed from frozen JPL catalog states matches the catalog value.
* **`testJplJacobiDriftEarthMoon.m`**: Verifies Jacobi conservation during Earth-Moon CR3BP propagation.
* **`testJplHaloClosureEarthMoon.m`**: Propagates Earth-Moon JPL halo states over one period and checks periodic closure.
* **`testMonodromyMatrixEarthMoon.m`**: Validates monodromy matrix properties, including determinant consistency and reciprocal eigenvalue pairing.
* **`testSingleShootingCorrectorEarthMoon.m`**: Validates the generic single-shooting corrector on a robust Earth-Moon L1 Lyapunov case.

---

### Integration Tests
Integration tests validate short multi-step workflows involving several functions together.

#### `tests/integration/lambert/`
* **`test_earth_mars_horizons.m`**: Tests an Earth-Mars transfer workflow against ephemeris-based reference conditions.

#### `tests/integration/cr3bp/`
* **`testMultipleShootingCorrectorEarthMoon.m`**: Validates multiple-shooting periodic correction on Earth-Moon halo cases.
* **`testContinuationPseudoArcEarthMoon.m`**: Validates pseudo-arclength continuation for Earth-Moon Lyapunov families.
* **`testLyapunovFamilyTrendEarthMoon.m`**: Compares locally continued Lyapunov-family trends against frozen JPL reference data.

---

### Validation Tests
Validation tests focus on external reference datasets and higher-level suite runners.

#### `tests/validation/`
* **`testJplCatalogData.m`**: Checks that frozen JPL periodic-orbit benchmark files exist and contain the expected schema.
* **`runCR3BPValidationSuite.m`**: Runs the main CR3BP validation suite in one command.

---

## Running the CR3BP Validation Suite

From MATLAB, in the repository root:

``` matlab
startup
results = runCR3BPValidationSuite;
```

## Benchmark Data

Several CR3BP tests rely on frozen benchmark files derived from the JPL periodic orbit catalog and stored under:

``` text
data/validation/cr3bp/
```

These datasets are used to keep validation deterministic and independent of live API access during testing.
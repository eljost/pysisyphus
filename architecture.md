# Architecture

This document describes the high-level architecture of pysisyphus.

## Bird's Eye View
pysisyphus is an excited states (ES) aware wrapper for several established quantum chemistry (QC) packages, like Gaussian, ORCA, Turbomole [...]. Energies, gradients and Hessian matrices, obtained from QC codes are utilized to optimize stationary points and reaction paths, by means of surface walking, chain of state (COS) methods
and intrinsic reaction coordinate (IRC) calculations.

The overall design is adapted from the well known [Atomic Simulation Environment](https://gitlab.com/ase/ase) (ASE). Central to ASE is the `Atoms` class, encapsulating calculators, that carry out the actual energy and gradient calculations.

Central to pysisyphus is the `Geometry` class, that works in the same way. `Geometry` objects wrap calculators and possibly other objects, to deal with coordinate system, beside Cartesians. For most time pysisyphus was developed in an object-oriented way (see optimizers and calculators). Newer additions (`dynamics` and `intcoords` pacakge) are implemented in a more functional way, to increase pysisyphus value as library.

## Code Map

The most important packages are found at the top, while the remaining packages follow below, sorted alphabetically.

### Important packages & classes

#### `calculators/`
Wrappers for external QC codes and implementations of anayltical 2D potentials.

Analytical potentials derive from `AnaPotBase.py` and are construced from mathematical functions via `sympy`. See `AnaPot.py` for a simple example.

Implementations for external QC codes are mainly derived from the `Calculator` base class. ES aware calculators derive from `OverlapCalculator`, which derives again from `Calculator`.

Calculators implement `get_energy()`, `get_forces()` and possibly `get_hessian()`. These methods all return an `results` dict, containing the appropriate entries. Energy is always present.
ES aware calculators derived from `OverlapCalculator` must implement `prepare_overlap_data()`, which returns energies of all states (ground and ESs), molecular orbital (MO) coefficients and transition density matrices of all ES.
Actual calculations are carried out via calls to the `Calculator.run()` method. Exceptions to this are calculators, having a native python interface (`PySCF` and `xtb-python`).

**Architecture invariant:** Calculation of energies & its derivatives is always requested with Cartesian coordinates, never with internal coordinates.


#### `Geometry.py`
Central class for handling coordinates and delegating calculations of energy and its derivatives to calculators. Defaults to Cartesian coordinates. Other coordinate systems are requesting by setting an appropriate `coord_type`. Implements methods that always return Cartesian quanteties, whereas attributes like `forces` and `hessian` return data in the `coord_type` coordinate system. All quantities are calculated at the current coordinates. Calculations at different coordinates, without storing them in a `Geometry` object are requested with `get_energy_at()` and `get_energy_and_forces_at()`. This is useful for line searches. While the `Geometry` state is not mutated, the state of the calculator probably is.


### Remaining packages & classes

#### `benchmarks`
Offers easy access to geometries and reference energies of several benchmark sets, like Bakers (TS) test set or the S22 set. The `Benchmark` class provides serveral iterators, allowing easy instantiation of `Geometry` objects with pre-set calculators. New test sets are added by extending `data.py` and adding the a new key/value pair to `Benchmark.data_funcs`.
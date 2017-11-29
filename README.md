# pysisphus
pysisphus implements Chain Of States methods like Nudged Elastic Band (NEB) and Simple Zero Temperature String to converge minimum energy paths and provide initial guesses for transition states. In addition pysisyphus features serveral Intrinsic Reaction Coordinate algorithms. The required gradients and/or hessians are calculated by external quantum chemistry codes, that can be called from pysisyphus.

**This software is still work in progress and shouldn't be used for production runs. Use at your own risk.**

## Geometry input and output
All coordinate input read from .xyz or .trj files is expected to be in Ångström. Internally atomic units (Hartree, Bohr, etc.) are used throughout. All coordinate output written is given in Ångström again.

## Implemented methods
### Chain Of States methods
* Initial guess
    * Linear
    * Image Dependent Pair Potential (https://doi.org/10.1063/1.4878664)
* Removal of translation and rotation (http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00360)
#### Nudged elastic band
* Improved tangent (https://doi.org/10.1063/1.1323224)
* Climbing Image (https://doi.org/10.1063/1.1329672)
* Different optimizers

#### Simple Zero Temperature String
* Steepest Descent Optimizer only (https://doi.org/10.1063/1.2720838)
* Splined images
    * Equal spacing
    * Energy dependent spacing, leading to higher resoultion at high energy regions

### Optimizers
* Can be used for single molecules or Chain Of States
* Currently implemented
    * Steepest descent
    * Conjugate Gradient
    * BFGS
    * QuickMin
    * FIRE
    * Wrapper for serveral SciPy optimizer (Experimental, as the step size can't be easily controlled)

### Intrinsic Reaction Coordinate
* Very experimental
* Mass-weighted coordinates
* Initial displacement using the transition vector
### Modified IMK
### Damped velocity verlet
### Gonzales-Schlegel 2
### Euler method
---
## Calculators
* Can be called in parallel are in serial mode
### OpenMolcas
* Wavefunctions from &rasscf and gradients from &alaska
* Reuse of orbitals from previous iterations
* State average calculations to converge excited state MEPs 

### GFN-xTB
* Fast semi-emperical method for ground state calculations (http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00118)

### ORCA
* All groundstate methods that can produce gradients and TDDFT

### Analytical 2D potentials
* All expressions that can be handled by sympy
* Automatic differentation to give gradients and hessian

## Roadmap
* NEB with two climbing images (https://doi.org/10.1063/1.4905209)


tmpdir:
export TMPDIR=[tmpdir]
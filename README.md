# pysisphus
pysisphus implements Chain Of States (COS) methods like Nudged Elastic Band (NEB) and Simple Zero Temperature String (SZTS) to converge minimum energy paths and provide initial guesses for transition states. In addition pysisyphus provides serveral Intrinsic Reaction Coordinate algorithms. The required gradients and/or hessians are calculated by calling external quantum chemistry codes. Everything is done in cartesian coordinates.

**This software is still work in progress and shouldn't be used for production runs. Use at your own risk.**

## Dependencies

	python3.6+
	scipy
	numpy
	matplotlib
	
Usage of pysisyphus within the **anaconda** distribution is recommended.

## Installation

	git clone https://github.com/eljost/pysisyphus
	cd pysisyphus
	python setup.py install|develop
	
or

	pip install pysisyphus
	

## Currently implemented methods
### Chain Of States methods
* Initial guess
    * Linear interpolation
    * Image Dependent Pair Potential (https://doi.org/10.1063/1.4878664)
* Removal of translation and rotation (http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00360)
* Either the first and/or the last image of the path can be fixed

#### Nudged elastic band
* Improved tangent (https://doi.org/10.1063/1.1323224)
* 1 Climbing Image (https://doi.org/10.1063/1.1329672)
* 2 Climbing Image (https://doi.org/10.1063/1.4905209)
* Different optimizers (see below)

#### Simple Zero Temperature String
* Steepest Descent Optimizer only (https://doi.org/10.1063/1.2720838)
* Splined images
    * Equal spacing
    * Energy dependent spacing, leading to higher resoultion at high energy regions

### Optimizers
* Can be used for single molecules or COS-paths
* Currently implemented
    * Steepest Descent
    * Conjugate Gradient
    * BFGS
    * QuickMin
    * FIRE
    * Wrapper for serveral SciPy optimizer that don't require an explicit hessian (Experimental and not recommended, as the step size can't be controlled easily)

### Intrinsic Reaction Coordinate
* Very experimental
* Mass-weighted or non-mass-weighted coordinates
* Initial displacement using the transition vector
### Modified IMK
* http://pubs.acs.org/doi/pdf/10.1021/ja00295a002
### Damped velocity verlet
* http://pubs.acs.org/doi/abs/10.1021/jp012125b
### Gonzales-Schlegel 2
* https://doi.org/10.1063/1.456010
### Euler method
* https://en.wikipedia.org/wiki/Euler_method

---

## Calculators
* Can be called in parallel or in serial mode
### OpenMolcas
* Wavefunctions from &rasscf and gradients from &alaska
* Reuse of orbitals from previous iterations
* State average calculations to converge excited state MEPs 

### GFN-xTB
* Fast semi-emperical method for ground state calculations (http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00118)

### ORCA
* All methods that have gradients
* Recognizes TDDFT optimizations

### Analytical 2D potentials
* All expressions that can be handled by sympy
* Automatic differentation to give gradients and hessian
* Easy visualization of optimization progress with COS methods (AnimPlot class)
---
## Geometry input and output
All coordinate input read from .xyz or .trj files is expected to be in Ångström. Internally atomic units (Hartree, Bohr, etc.) are used throughout. All coordinate output written is given in Ångström again.

## Roadmap
* Get everything to work very smoothly;)

## Special 

tmpdir:
export TMPDIR=[tmpdir]

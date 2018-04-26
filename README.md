# pysisphus
pysisphus implements Chain Of States (COS) methods like Nudged Elastic Band (NEB) and Simple Zero Temperature String (SZTS) to converge minimum energy paths and provide initial guesses for transition states. In addition pysisyphus provides serveral Intrinsic Reaction Coordinate algorithms. The required gradients and/or hessians are calculated by calling external quantum chemistry codes. By default everything is done in cartesian coordinates but an internal coordinates implementation is in progress.

**This software is still work in progress and shouldn't be used for production runs. Use at your own risk.**

## Dependencies

	python3.6+
	scipy
	numpy
	matplotlib
	pytest
	natsort
	pyyaml
	
Usage of pysisyphus within the **anaconda** distribution is recommended.

## Installation
### Install the code
The **recommended way**, as the version on pypi may be outdated:

    cd [root dir where you want to keep pysisyphus]
	git clone https://github.com/eljost/pysisyphus
	cd pysisyphus
	python setup.py develop
	
or **(not recommended)**:

	pip install pysisyphus

### Tell pysisyphus how to call your calculators	
Create a file called `.pysisyphusrc` in your $HOME-directory containing the commands the call your installed calculators, e.g.:

    [orca]
    cmd=/scratch/programme/orca-4.0.1.2/orca
    
    [xtb]
    cmd=xtb
    
    [openmolcas]
    cmd=pymolcas
    
    [gaussian09]
    cmd=/usr/local/gaussian/g09/g09
    
    [wfoverlap]
    cmd=/scratch/wfoverlap_1.0/bin/wfoverlap.x

When the specified path is not absolute but relative (as in the example config for xtb and openmolcas) the corresponding binaries have to be available on the $PATH! The **Turbomole** binaries are hardcoded and don't have to be specified here. 

## Usage
pysisyphus provides three hooks that can be called from the shell (command line). The available commandas can be queried with the `-h` or `--help` arguments:

    # Run optimizations of single molecules or chain of states (NEB, string method)
    pysisrun
    # Visualize the results of chain of states calculations
    pysisplot
    # Manipulate a .trj file or multiple .xyz files
    pysistrj

### pysisrun
Usually called with a YAML-formatted file containing all input for the optimizations. The --restart is quite experimental for now.

    usage: pysisrun [-h] [--clean] [--fclean] [--restart] [yaml]

    positional arguments:
     yaml        Start pysisyphus with input from a YAML file.

    optional arguments:
     -h, --help  show this help message and exit
     --clean     Ask for confirmation before cleaning.
     --fclean    Force cleaning without prior confirmation.
     --restart   Continue a previously crashed/aborted/... pysisphus run.
     
#### Example for a simple xTB NEB run

    cos:
     type: neb          # optimize a NEB
     parallel: 4        # Run 4 calculations in parallel
    opt:
     type: cg           # Use conjugate gradient
     max_cycles: 50     # Optimize for at most 50 cycles
     align: True        # Remove translation and rotation with partial procrustes
     climb: True        # Use climbing image NEB to get a better guess for the TS
    interpol:
     idpp: True         # Interpolate using the image dependent pair potential (recommended)
     between: 15        # Interpolate 15 images for a total of (1 + 15 + 1 = 17) images
    calc:
     type: xtb          # Use xTB calculator
     charge: 0
     mult: 1
    # Interpolate between these two geometries. One can supply also more than two geometries,
    # but then the between value in the interpol section has to be lowered. Supplying a third
    # geoemtry here would lead to a chain of states containing 33 images (1 + 15 + 1 + 15 + 1). 
    xyz: [15_tpss_tzvp_full_scan_acn_013_xtb.xyz, 15_tpss_tzvp_full_scan_acn_024_xtb.xyz]

### pysisplot
    usage: pysisplot [-h] [--until UNTIL]
                     (--saras | --tddft | --params PARAMS | --cosgrad | --energies)

    optional arguments:
        -h, --help       show this help message and exit
        --until UNTIL    Only show until cycle [until].
        --saras          Plot OpenMolcas state average potential energy surfaces
                         over the course of the NEB.
        --tddft          Plot ORCA TDDFT potential energy surfaces over the course
                         of the NEB.
        --params PARAMS  Follow internal coordinates over the course of the NEB. All
                         atom indices have to be 0-based. Use two indices for a
                         bond, three indices for an angle and four indices for a
                         dihedral. The indices for different coordinates have to be
                         separated by ','.
        --cosgrad        Plot image gradients along the path.
        --energies       Plot energies.


### pysistrj
Usually called wit multiple .xyz files or a single .trj file.

    usage: pysistrj [-h] [--idpp] [--xyz XYZ [XYZ ...]]
                (--between BETWEEN | --align ALIGN [ALIGN ...] | --split SPLIT | --reverse REVERSE | --cleantrj CLEANTRJ | --spline SPLINE)

    optional arguments:
        -h, --help            show this help message and exit
        --idpp                Use Image Dependent Pair Potential instead of simple
                              linear interpolation.
        --xyz XYZ [XYZ ...]
        --between BETWEEN     Interpolate additional images.
        --align ALIGN [ALIGN ...]
                              Align geometries onto the first geometry read from
                              multiple .xyz or one .trj file.
        --split SPLIT         Split a supplied .trj file in multiple .xyz files.
        --reverse REVERSE     Reverse a .trj file.
        --cleantrj CLEANTRJ   Clean a .trj file.
        --spline SPLINE       Evenly redistribute geometries along a splined path.
	

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
    * Using only first derivatives
        * Steepest Descent (SD)
        * Conjugate Gradient (CG)
        * QuickMin (QM)
        * FIRE
    * Using second derivatives (Hessian)
        * BFGS
        * Rational Function Optimization (RFO)
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

### Turbomole
* Density functional methods 
* Recognizes TDDFT optimizations

### Gaussian
* (Probably) All methods that support Force run
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

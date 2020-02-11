# pysisphus
[![Documentation Status](https://readthedocs.org/projects/pysisyphus/badge/?version=dev)](https://pysisyphus.readthedocs.io/en/dev/?badge=dev)
![build](https://github.com/eljost/pysisyphus/workflows/Python%20application/badge.svg)
[![codecov](https://codecov.io/gh/eljost/pysisyphus/branch/master/graph/badge.svg)](https://codecov.io/gh/eljost/pysisyphus)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/96281078.svg)](https://zenodo.org/badge/latestdoi/96281078)

pysisyphus is a software-suite for the exploration of potential energy surfaces in ground- and **excited states**. It implements several methods to search for stationary points (minima and first order saddle points) and calculation of minimum energy paths by means of IRC and Chain of States methods like Nudged Elastic Band and Growing String. Furthermore it supports interpolation of geometries in internal coordinates.

The required energies, gradients and hessians are calculated by calling external quantum chemistry codes. pysisyphus can also be used as a library to implement custom quantum chemistry workflows.

**This software is still work in progress. Use at your own risk. Also take a look at the [license](https://github.com/eljost/pysisyphus/blob/master/LICENSE).**

Contrubtions are welcome.

If any issues arise please open an [issue](https://github.com/eljost/pysisyphus/issues) and I'll try to fix it if possible and my time permits it.

## Usage
pysisyphus provides several entry points that can be called from the shell (command line). The available commands of the entry points can be queried with the `-h` or `--help` arguments:

    # Run calculations (Minima optimization, TS search, IRC, NEB, GS, ...)
    pysis
    # Plotting of path-energies, optimization progress, IRC progress, etc ...
    pysisplot
    # Manipulate a .trj file or multiple .xyz files
    pysistrj

### pysis
Called with a YAML-input file. Simple and more complex examples can be found in the `examples` directory of this repository. Some common arguments are given below:

    usage: pysis [-h] [--clean] [--dryrun | --cp CP]
             [yaml]

    positional arguments:
        yaml                  Start pysisyphus with input from a YAML file.

    optional arguments:
        --clean               Ask for confirmation before cleaning.
        --dryrun              Only generate a sample input (if meaningful) for
                              checking.
        --cp CP, --copy CP    Copy .yaml file and corresponding geometries to a new
                              directory. Similar to TURBOMOLEs cpc command.

### pysisplot
Visualization/plotting of the pysisyphus calculations. Some common arguments are given below:

    usage: pysisplot [-h] [--first FIRST] [--last LAST] [--h5 H5]
                 [--orient ORIENT]
                 (--cosgrad | --energies | --aneb | --all_energies | --bare_energies | --afir | --opt | --irc)

    optional arguments:
        -h, --help           show this help message and exit
        --cosgrad, --cg      Plot image gradients along the path.
        --energies, -e       Plot energies.
        --aneb               Plot Adaptive NEB.
        --all_energies, -a   Plot ground and excited state energies from
                             'overlap_data.h5'.
        --bare_energies, -b  Plot ground and excited state energies from
                             'overlap_data.h5'.
        --afir               Plot AFIR and true -energies and -forces from an AFIR
                             calculation.
        --opt                Plot optimization progress.
        --irc                Plot IRC progress.

### pysistrj
Please see `pysistrj --help` for a list of available arguments.

## Available calculators
### Excited state capabilities
|Program|Gradient|Hessian|Exc. states|Comment|
|-|-|-|-|-|
|ORCA|y|y|TD-DFT, TDA||
|Turbomole|y|y|TD-DFT, TDA, ricc2|Tested with V7.2, V7.3 and V7.4. Only ground state hessians.|
|Gaussian16|y|y|tested TD-DFT, TDA|
|PySCF|y|y|tested TD-DFT|
|OpenMOLCAS|y|-|&rasscf|Not derived from OverlapCalculator, so functionality may lag behind for now.|
### Ground states only
|Program|Gradient|Hessian|Comment|
|-|-|-|-|
|MOPAC2016|y|y||
|XTB|y|y|Tested with 6.2.2|
|QCEngine|depends|depends|Capabilities depend on the program harness.|
|Psi4|y|y||
### Pure python calculators & Wrappers
|Program|Gradient|Hessian|Comment|
|-|-|-|-|
|Sympy 2d potentials|y|y|Many analytical potentials with gradients & hessian from sympy (LEPS, Rosenbrock, Cerjan-Miller, ...).|
|Lennard-Jones|y|-|No periodic boundary conditions.|
|**AFIR**|y|-|Wrapper for other pysisyphus calculators|
|**ONIOM**|y|-|Arbitrary levels with multicenter-support in the highest level. Wraps other pysisyphus calculators|
|FakeASE|y|-|Thin bare-bones wrapper for pysisyphus calculators so they can be used with ASE.|

## Available algorithms
### Chain Of State (COS) Methods
| Algorithm |Coordinate system| Comment | Links |
|-----------|---------|-------|-|
|Nudged elastic band (NEB)|Cartesian|DLC planned, Climbing Image variants|
|Adaptive NEB|Cartesian|Not well tested|
|Free-End NEB|Cartesian|Not well tested|
|Simple zero temperature string|Cartesian|Equal-, and energy dependent spacing||
|Growing String method|Cartesian, **Delocalized internals**||
### COS Optimizers
| Algorithm | Comment | Links |
|-----------|---------|-------|
|Steepest Descent (SD)|Backtracking variant||
|Conjugate Gradient (CG)|Backtracking variant||
|QuickMin (QM)|||
|FIRE|||
|BFGS||Removal of translation & rotation is disabled.|

### Transition state optimization
| Algorithm | Comment | Links |
|-----------|---------|-------|
|RS-P-RFO| default choice |https://pubs.acs.org/doi/pdf/10.1021/j100247a015, https://link.springer.com/article/10.1007/s002140050387|
|RS-I-RFO| RFO, Image function |https://link.springer.com/article/10.1007/s002140050387|
|TRIM|Trust region Image function|https://doi.org/10.1016/0009-2614(91)90115-P|
|Dimer method|Hessian free TS search||


### Intrinsic Reaction Coordinate
* Done in mass-weighted cartesian coordinates
* Initial displacement from prescribed energy lowering along the transition vector or by length

| Algorithm | Comment | Links |
|-----------|---------|-------|
|Modified EulerPC| default choice | https://aip.scitation.org/doi/pdf/10.1063/1.3514202, https://doi.org/10.1039/C7CP03722H |
|Modified IMK||http://pubs.acs.org/doi/pdf/10.1021/ja00295a002|
|Damped velocity verlet||http://pubs.acs.org/doi/abs/10.1021/jp012125b|
|Gonzales-Schlegel 2||https://doi.org/10.1063/1.456010|
|LQA||https://aip.scitation.org/doi/pdf/10.1063/1.459634?class=pdf|
|Euler method|Not advisable|https://en.wikipedia.org/wiki/Euler_method|
|Runge-Kutta-4|Not advisable|https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Derivation_of_the_Runge%E2%80%93Kutta_fourth-order_method|

## Additional remarks

    tmpdir:
    export TMPDIR=[tmpdir]

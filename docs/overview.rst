Overview of pysisyphus
**********************

`pysisyphus` is a software-suite for the exploration of potential energy surfaces in ground-
and **excited states**. User input is read from YAML files, but it can also be used as a python
library to set up custom workflows.

Below you can find a screencast of the transition state (TS) search for the famous alanine dipeptide
isomerization at the `xtb` level of theory. It starts with pre-optimizations of the initial
and final geometry, a subsequent growing string calculation in delocalized internal coordinates
to generate a guess for the final TS optimization and with which it concludes.

.. image:: https://asciinema.org/a/300731.png
    :target: https://asciinema.org/a/300731

Besides being able to track excited states over the course of an optimization the idea of `pysisyphus` is to provide tools for handling the whole process of calculating reaction paths, starting from pre-optimizing the given input geometries for chain of states methods and ending with the hessian calculation on the optimized IRC endpoints.

Entry points
============

pysisyphus provides several entry points that can be called from the shell (command line). The available commands of the entry points can be queried with the `-h` or `--help` arguments:

.. code:: bash

    # Run calculations (Minima optimization, TS search, IRC, NEB, GS, ...)
    pysis
    # Plotting of path-energies, optimization progress, IRC progress, etc ...
    pysisplot
    # Manipulate a .trj file or multiple .xyz files
    pysistrj

pysis
-----
Called with a YAML-input file. Simple and more complex examples can be found in the `examples` directory of this repository. Some common arguments are given below:

.. code:: bash

    usage: pysis [-h] [--clean] --cp CP]
             [yaml]

    positional arguments:
        yaml                  Start pysisyphus with input from a YAML file.

    optional arguments:
        --clean               Ask for confirmation before cleaning.
        --cp CP, --copy CP    Copy .yaml file and corresponding geometries to a new
                              directory. Similar to TURBOMOLEs cpc command.

pysisplot
---------
Visualization/plotting of the pysisyphus calculations. Some common arguments are given below:

.. code:: bash

    usage: pysisplot (--cosforces | --cosens | --all_energies | --afir | --opt | --irc | --overlaps)

    optional arguments:
      --cosforces, --cf    Plot image forces along a COS.
      --cosens, --ce       Plot COS energies.
      --all_energies, -a   Plot ground and excited state energies from 'overlap_data.h5'.
      --afir               Plot AFIR and true -energies and -forces from an AFIR calculation.
      --opt                Plot optimization progress.
      --irc                Plot IRC progress.
      --overlaps, -o


pysistrj
--------
Please see `pysistrj --help` for a list of available arguments.

Available calculators
=====================

Excited state capabilities
--------------------------

=========== ======== ======= ================== =======
Program     Gradient Hessian Exc. states        Version
=========== ======== ======= ================== =======
ORCA        y        y       TD-DFT, TDA        4.2.x
Turbomole   y        y       TD-DFT, TDA, ricc2 7.2-7.4
Gaussian16  y        y       tested TD-DFT, TDA 
PySCF       y        y       tested TD-DFT      1.7.5
OpenMOLCAS  y        n       &rasscf            Not yet derived from OverlapCalculator
=========== ======== ======= ================== =======

Ground states capabilities
--------------------------

=========== ======== ======= =========
Program     Gradient Hessian Version
=========== ======== ======= =========
MOPAC2016   y        y       
XTB         y        y       6.3.2
QCEngine    y        n       >= 0.16.0
Psi4        y        y 
=========== ======== ======= =========

Pure python calculators & Wrappers
----------------------------------

============= ======== ======= =========
Program       Gradient Hessian Comment
============= ======== ======= =========
Sympy 2D      y        y       Many analytical potentials (LEPS, Rosenbrock, Cerjan-Miller,
                               Muller-Brown, ...)
Lennard-Jones y        n       **No** periodic boundary conditions
AFIR          y        n       
ONIOM         y        n       Arbitrary number of layers with multicenter-support in the highest layer.
FakeASE       y        n       Wraps `pysisyphus` calculators so they can be used with `ase`.
============= ======== ======= =========

Available algorithms
=====================

Chain Of States Methods
-----------------------

=============================== ====================== =======
Algorithm                       Coordinates            Comment
=============================== ====================== =======
Nudged Elastic Band (NEB)       Cartesian, DLC planned Climbing Image variants, Doubly nudged variant
Adaptive NEB                    Cartesian              Not well tested
Free-End NEB                    Cartesian              Not well tested
Simple Zero-Temperature-String  Cartesian              Equal spacing, energy-dependent spacing
**Growing String Method**       Cartesian, **DLC**
=============================== ====================== =======

Chain Of States Optimizer
--------------------------

================== ==================== =======
Algorithm          Comment              Links
================== ==================== =======
Steepest Descent   Backtracking variant NEB-Optimizers_
Conjugate Gradient Backtracking variant NEB-Optimizers_
QuickMin                                NEB-Optimizers_
FIRE                                    NEB-Optimizers_
BFGS                                    NEB-Optimizers_
================== ==================== =======

.. _NEB-Optimizers: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00360

Transition state optimization
-----------------------------

================== ==================== =======
Algorithm          Comment              Links
================== ==================== =======
RS-P-RFO           default              RFO-Paper_, RS-Paper_
RS-I-RFO                                RFO-Paper_, RS-Paper_
TRIM                                    TRIM-Paper_ 
Dimer method                                    
================== ==================== =======

.. _RFO-Paper: https://pubs.acs.org/doi/pdf/10.1021/j100247a015
.. _RS-Paper: https://link.springer.com/article/10.1007/s002140050387
.. _TRIM-Paper: https://doi.org/10.1016/0009-2614(91)90115-P

Intrinsic Reaction Coordinate integrators
-----------------------------------------

============================= ==================== =======
Algorithm                     Comment              Links
============================= ==================== =======
Damped-Velocity-Verlet                             DVV-Paper_
Euler                         Not recommended
EulerPC                       default              Kaestner-PC_, Euler-PC_
Gonzales-Schlegel 2                                GS2-Paper_
Local Quadratic Approximation                      LQA-Paper_
Modified IMK                                       IMK-Paper_
Runge-Kutta-4                 Not recommended
============================= ==================== =======

.. _Kaestner-PC: https://doi.org/10.1039/C7CP03722H
.. _Euler-PC: https://aip.scitation.org/doi/pdf/10.1063/1.3514202
.. _IMK-Paper: http://pubs.acs.org/doi/pdf/10.1021/ja00295a002
.. _DVV-Paper: http://pubs.acs.org/doi/abs/10.1021/jp012125b
.. _GS2-Paper: https://doi.org/10.1063/1.456010
.. _LQA-Paper: https://aip.scitation.org/doi/pdf/10.1063/1.459634?class=pdf

Additional remarks
==================

`pysisyphus` uses the `tempfile` module from the python stdlib. The location of the temporary
directories can be controlled by setting the **$TMPDIR** environment variable before
executing `pysis`.

.. code:: bash

    export TMPDIR=[tmpdir]

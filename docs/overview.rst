Overview of pysisyphus
**********************

Usage
=====

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

pysisplot
---------
Visualization/plotting of the pysisyphus calculations. Some common arguments are given below:

.. code:: bash

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
PySCF       y        y       tested TD-DFT      1.7.0
OpenMOLCAS  y        n       &rasscf
=========== ======== ======= ================== =======

Ground states capabilities
--------------------------

=========== ======== ======= =========
Program     Gradient Hessian Version
=========== ======== ======= =========
MOPAC2016   y        y       
XTB         y        y       6.2.2
QCEngine    y        n       >= 0.13.0
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


Additional remarks
==================

`pysisyphus` the `tempfile` module from the stdlib. The location of the temporary
directories can be controlled by setting the **$TMPDIR** environment variable before
executing `pysis`.

.. code:: bash

    export TMPDIR=[tmpdir]

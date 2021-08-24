Calculators
***********

Exploring potential energy surfaces (PESs) requires calculation of energies
and its derivatives (gradients, Hessian matrices).

For testing purposes and method development, pysisyphus implements various analytical
2d potentials, as they allow fast evaluations of the aforementioned quantities.
Actual (production) calculations are carried out by wrapping existing quantum chemistry
codes (ORCA, TURBOMOLE, Gaussian etc.) and delegating the required calculations to them.
Pysisyphus generates all necessary inputs, executes the QC code and parses their output
for the requested quantities.

Furthermore pysisyphus provides several "meta"-calculators which wrap other (actual) calculators,
to modify calculated energies and forces. Examples for this are the `Dimer`_ calculator,
used for carrying out transition state searches or the `ONIOM`_ calculator,
allowing  multi-level calculations comprising different levels of theory.

External forces, e.g. a restraining spherical potential or harmonic restraints on
primitive internal coordinates (stretches, bends, torsion) can be applied with
`ExternalPotential`_.

YAML input
==========

Possible keywords for the YAML input can be derived from inspecting the possible arguments
of the `Calculator` base-class (see below) and the possible arguments of the respective calculator,
e.g. ORCA or XTB.

The most commonly used keywords, derived from the `Calculator` baseclass are `mem`,
handling the requested memory per core in MB, `pal`, handling the number of requested CPU cores,
`charge`, the total charge of the system and `mult`, the systems multiplicity.

For (excited state, ES) calculations carried out with calculators derived from `OverlapCalculator`
additional keywords are possible. The most common keywords controlling ES calculations are `track`, activating
ES tracking, `ovlp_type`, selecting the tracking algorithm and `ovlp_with`, handling the selection of the
reference cycle.

An example input highlighting the most important keywords is shown below.

.. code-block:: yaml

    geom:
     [... omitted ...]
    calc:
     type: orca                     # Calculator type, e.g. g09/g16/openmolcas/
                                    # orca/pyscf/turbomole/dftb+/mopac/psi4/xtb

     pal: 1                         # Number of CPU cores
     mem: 1000                      # Memory per core
     charge: 0                      # Charge
     mult: 1                        # Multiplicity
     # Keywords for ES calculations
     track: False                   # Activate ES tracking
     ovlp_type: tden                # Tracking algorithm
     ovlp_with: previous            # Reference cycle selection
     # Additional calculator specific keywords
     [... omitted ...]
    opt:
     [... omitted ...]


Calculator base classes
=======================

.. automodule:: pysisyphus.calculators.Calculator
    :members:
    :undoc-members:

OverlapCalculator base class
============================
.. automodule:: pysisyphus.calculators.OverlapCalculator
    :members:
    :undoc-members:
    :show-inheritance:


Calculators with Excited state capabilities
===========================================

Gaussian09
----------
.. automodule:: pysisyphus.calculators.Gaussian09
    :members:
    :undoc-members:

Gaussian16
----------
.. automodule:: pysisyphus.calculators.Gaussian16
    :members:
    :undoc-members:
    :show-inheritance:

OpenMolcas
----------
.. automodule:: pysisyphus.calculators.OpenMolcas
    :members:
    :undoc-members:
    :show-inheritance:

ORCA 4.2.1 / 5.0.1
----------
.. automodule:: pysisyphus.calculators.ORCA
    :members:
    :undoc-members:
    :show-inheritance:

PySCF 1.7.6
-----------
.. automodule:: pysisyphus.calculators.PySCF
    :members:
    :undoc-members:
    :show-inheritance:

Turbomole 7.x
-------------
Pysisyphus does not implement a wrapper for `define`, so the user has to
manually prepare a directory containing a valid control file.
An automated `define` wrapper, restricted to ground state functionality, is
available via the `QCEngine project <https://github.com/MolSSI/QCEngine>`_,
to which I contributed the Turbomole harness.

Care should be taken to include only the minimum amount of necessary files in the
`control_path` directory, e.g., `(auxbasis, basis, control, coord, mos)` for a
closed-shell calculation using RI. **A gradient file must not be present
in control_path, as well as other subdirectories and files with .out extension.**
The `coord` file, while not strictly required, should be kept too, to facilitate
testing of the setup with standalone Turbomole.

It may be a good idea to pre-converge the calculation in `control_path`,
to see if the setup is correct and actually works. Resulting files like
`energy`, `statistics` can be deleted; `mos` should be kept, as the converged
MOs are reused in pysisyphus.

If an excited-state optimization is desired, care has to be taken, to include
**$exopt [n]** for TD-DFT/TDA or the **geoopt state=([n])** (ricc2)! Tracking
of excited states is currently possible for closed shell `egrad` and `ricc2` calculations.

The current implementation was tested against Turbomole 7.4.1 and QCEngine 0.19.0.
Please see `examples/complex/11_turbomole_gs_tsopt <https://github.com/eljost/pysisyphus/tree/master/examples/complex/11_turbomole_gs_tsopt>`_ for a full example where Turbmole is
utilized in a growing string calculation. The same example, using QCEngine, is found
in `examples/complex/12_qcengine_turbomole_gs_tsopt <https://github.com/eljost/pysisyphus/tree/master/examples/complex/12_qcengine_turbomole_gs_tsopt>`_. MOs are not reused with the
QCEngine calculator, so the native pysisyphus calculator is probably faster.

.. automodule:: pysisyphus.calculators.Turbomole
    :members:
    :undoc-members:
    :show-inheritance:


DFTB+ 20.x
-----------
.. automodule:: pysisyphus.calculators.DFTBp
    :members:
    :undoc-members:
    :show-inheritance:

Calculators with Ground state capabilities
==========================================

MOPAC 2016
----------
.. automodule:: pysisyphus.calculators.MOPAC
    :members:
    :undoc-members:
    :show-inheritance:

Psi4
----
.. automodule:: pysisyphus.calculators.Psi4
    :members:
    :undoc-members:
    :show-inheritance:

QCEngine
--------
.. automodule:: pysisyphus.calculators.QCEngine
    :members:
    :undoc-members:
    :show-inheritance:

XTB 6.x
-------
.. automodule:: pysisyphus.calculators.XTB
    :members:
    :undoc-members:
    :show-inheritance:

Dalton
-------
.. automodule:: pysisyphus.calculators.Dalton
    :members:
    :undoc-members:
    :show-inheritance:

OpenBabel
-------
.. automodule:: pysisyphus.calculators.OBabel
    :members:
    :undoc-members:
    :show-inheritance:



Meta (wrapping) Calculators
===========================

ExternalPotential
-----------------
.. automodule:: pysisyphus.calculators.ExternalPotential.ExternalPotential
    :members:
    :show-inheritance:

Restraint
---------

.. code:: yaml
    
    # General input structure for restraints
    calc:
     type: ext
     # Multiple potentials could be specified here as a list
     potentials:
       # Primitive internal coordinate restraint
       - type: restraint
         # List of restraints; could also be multiple restraints. Each restraint is given as
         # list of 2 or 3 items.
         #
         # The first item always specifies an internal coordinate,
         # whereas the second argument is a force constant (in atomic units; actual units
         # depend on the coordinate). Optionally a reference value (third argument) can be
         # given. If omitted, the initial coordinate value is used as reference value.
         restraints: [[[BOND, 0, 1], 10, 3.0]]
         # The commented out input below would restrain the bond at its initial value.
         #restraints: [[[BOND, 0, 1], 10]]
         # Multiple restraints are specified as given below.
         #restraints: [[[BOND, 0, 1], 10], [[BEND, 0, 1, 2], 1.0]]
    calc:
     type: [actual calculator that is wrapped by ExternalPotential]

.. automodule:: pysisyphus.calculators.ExternalPotential.Restraint
    :members:
    :undoc-members:

HarmonicSphere
--------------

.. automodule:: pysisyphus.calculators.ExternalPotential.HarmonicSphere
    :members:
    :undoc-members:

LogFermi
--------

.. automodule:: pysisyphus.calculators.ExternalPotential.LogFermi
    :members:
    :undoc-members:

AFIR
----
.. automodule:: pysisyphus.calculators.AFIR
    :members:
    :undoc-members:
    :show-inheritance:

ONIOM
-----
.. automodule:: pysisyphus.calculators.ONIOMv2
    :members:
    :undoc-members:
    :show-inheritance:

Dimer
-----
.. automodule:: pysisyphus.calculators.Dimer
    :members:
    :show-inheritance:
    :undoc-members:

Pure Python calculators
=======================

Sympy 2D Potentials
-------------------
.. automodule:: pysisyphus.calculators.AnaPotBase
    :members:
    :undoc-members:
    :show-inheritance:

Lennard-Jones
-------------
.. automodule:: pysisyphus.calculators.LennardJones
    :members:
    :undoc-members:
    :show-inheritance:

FakeASE
-------
.. automodule:: pysisyphus.calculators.FakeASE
    :members:
    :undoc-members:
    :show-inheritance:

TIP3P
-----
.. automodule:: pysisyphus.calculators.TIP3P
    :members:
    :undoc-members:
    :show-inheritance:

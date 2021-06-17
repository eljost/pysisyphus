Calculators
***********

Exploring potential energy surfaces (PESs) requires calculation of energies
and its derivatives (gradients, Hessian matrices). Pysisyphus implements various
analytical 2d potentials for testing purposes, as they allow fast evaluations.
Actual (production) calculations are carried out by wrapping existing quantum chemistry
codes (ORCA, TURBOMOLE, Gaussian etc.) and delegating the required calculations to them.
Pysisyphus generates all necessary inputs, exectes the code and parses the outputs
for the required quantities.

Furthermore pysisyphus provides several "meta"-calculators which wrap other calculators
to modify calculated energies and forces. Examples for this are the `Dimer`_ calculator,
used for carrying out transition state searches or the `ONIOM`_ calculator,
allowing  multi-level calculations comprising different levels of theory.

External forces, e.g. a restraining spherical potential or harmonic restraints on
primitive internal coordinates (stretches, bends, torsion) can be applied with
`ExternalPotential`_.

Calculators with Excited state capabilities
===========================================

Calculators with excited state (ES) capabilities inherit from OverlapCalculator.

OverlapCalculator base class
----------------------------
.. automodule:: pysisyphus.calculators.OverlapCalculator
    :members:
    :undoc-members:
    :show-inheritance:

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

ORCA 4.2.1
----------
.. automodule:: pysisyphus.calculators.ORCA
    :members:
    :undoc-members:
    :show-inheritance:

PySCF 1.7.x
-----------
.. automodule:: pysisyphus.calculators.PySCF
    :members:
    :undoc-members:
    :show-inheritance:

Turbomole 7.x
-------------
For now I have chosen the "easy" way and didnt't try to implement a wrapper for
`define`. That's why the user has to manually prepare a valid TURBOMOLE job directory
beforehand, which is then supplied to the calculator via the `control_path` argument.

If an excited-state optimization is desired care has to be taken to include
**$exopt [n]** for TD-DFT/TDA or the **geoopt state=([n])** (ricc2)! Tracking
of excited states is currently possible for closed shell `egrad` and `ricc2` calculations.

**Right now care has to be taken that no** `gradient` **file is present in the** `control_path`!

An easier, alternative way to use TURBOMOLE in `pysisyphus` is via its `QCEngine` wrapper,
albeit with restricted funtionality (no excited states right now).

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

Calculator base class
=====================

.. automodule:: pysisyphus.calculators.Calculator
    :members:
    :undoc-members:

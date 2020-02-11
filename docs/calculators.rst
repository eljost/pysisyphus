Calculators
***********

Calculator base class
=====================

.. automodule:: pysisyphus.calculators.Calculator
    :members:
    :undoc-members:

Excited state capabilities
==========================

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

PySCF 1.7.0
-----------

.. automodule:: pysisyphus.calculators.PySCF
    :members:
    :undoc-members:
    :show-inheritance:

Turbomole 7.4
-------------

For now I have chosen the "easy" way and didnt't try to implement a wrapper for
`define`. That's why the user has to manually prepare a valid TURBOMOLE job directory
beforehand that is then supplied to the calculator via the `control_path` argument.

If an excited-state optimization is desired care has to be taken to include
**$exopt [n]** for TD-DFT/TDA or the **geoopt state=([n])** (ricc2)! Tracking
of excited states is currently possible for closed shell `egrad` and `ricc2` calculations.

**Right now care has to be taken that no** `gradient` **file is present in the** `control_path`!

An easier, alternative way to use TURBOMOLE in `pysisyphus` is via its `QCEngine` wrapper,
albeit with restricted funtionality (no hessian, no excited states right now).

.. automodule:: pysisyphus.calculators.Turbomole
    :members:
    :undoc-members:
    :show-inheritance:

Ground state capabilities
==========================

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

XTB 6.2
-------

.. automodule:: pysisyphus.calculators.XTB
    :members:
    :undoc-members:
    :show-inheritance:

Pure Python calculators & wrappers
==================================

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

FakeASE
---------------------------------

.. automodule:: pysisyphus.calculators.FakeASE
    :members:
    :undoc-members:
    :show-inheritance:

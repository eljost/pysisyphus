Coordinate Systems
******************

Choice of coordinate system is often crucial for the successful outcome of
geometry optimizations. Pysisyphus supports different coordinate systems.
By default, molecular optimizations are carried out in redundant internal coordinates
(RIC).

Cartesian Coordinates
---------------------

* Strong coupling
* Unambiguously defined, if translation and rotation (TR) are removed
* Redundant set, if TR are not removed

Redundant Internal Coordinates
------------------------------

* Comprising bond stretches, (linear) bends, (improper) torsions and other primitives
* Less coupling, compared to Cartesians
* Easy estimation of diagonal model Hessians
* **Support constraints**
* Require sophisticated setup algorithm
* Internal-Cartesian backtransformation required, to transform steps from internals to Cartesians
* Redundant set

Delocalized Internal Coordinates (DLC)
--------------------------------------

* Complicated linear combinations of primitive internals
* No coupling, coordinates are orthogonal to each other, at least at the geometrey they were defined at
* Non redundant set
* Same comments apply as for RICs

Supported File Formats
======================

All formats can be read, at least to a certain extend, by pysisyphus. Coordinate
files are expected to be in Å, if not otherwise dictated by the respective format.

================ ===== =================================
Suffix           Write   Comment            
================ ===== =================================
.xyz             ✓     Plain XYZ, single geometry.
.trj             ✓     Plain XYZ, multiple geometries.
.pdb             ✓     Protein Data Bank.
.molden          ✗     Restricted to [Geometries] block.
.zmat            ✗     Z-Matrix, see below for an example.
.cjson           ✗     As saved by Avogadro.
================ ===== =================================

Z-Matrix example
----------------

.. code:: text

    C
    N 1 1.35
    H 1 1.0 2 105.0
    H 2 1.4 1 105.0 3 150.0
    H 2 1.4 1 110.0 3 -160.0

Indexing in Z-matrices is 1-based, so the first atom has index 1. Variable substitution
is not supported.

YAML Input
==========

.. code:: yaml

    geom:
     type: cart               # Coordinate system (cart/redund/dlc)
     fn: [input]              # File name or inline input
     define_prims: null       # Additional primitives, to be defined
     constrain_prims: null    # Primitive internals to be constrained
     union: null              # Define same set of primitives at multiple geometries
     isotopes: null           # Specify different isotopes
    preopt:
     geom:                    # geom block in preopt takes same keywords as above
      ...
    endopt:
     geom:                    # geom block in endopt takes same keywords as above
      ...

Inline input
-------------
Inline

Types of Primitive Coordinates
------------------------------
PrimTypes

Define Additional Primitives
----------------------------
DefinePrims

Constraints
-----------
Constraints

Related Literature
==================

1. `The efficient optimization of molecular geometries using redundant internal coordinates <https://doi.org/10.1063/1.1515483>`_
2. `The generation and use of delocalized internal coordinates in geometry optimization <https://doi.org/10.1063/1.471864>`_
3. `The efficient optimization of molecular geometries using redundant internal coordinates <https://doi.org/10.1063/1.1515483>`_
4. `Geometry optimization in redundant internal coordinates <https://doi.org/10.1063/1.462844>`_
5. `Geometry optimization made simple with translation and rotation coordinates <https://doi.org/10.1063/1.4952956>`_



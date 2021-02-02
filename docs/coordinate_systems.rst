Coordinate Systems
******************

Choice of coordinate system is often crucial, for the successful outcome of
geometry optimizations. Pysisyphus supports different coordinate systems.
By default, molecular optimizations are carried out in redundant internal coordinates
(RIC).

Supported Coordinate Systems
----------------------------

Cartesian Coordinates
^^^^^^^^^^^^^^^^^^^^^

* Strong coupling
* Unambiguously defined, if translation and rotation (TR) are removed
* Redundant set, if TR are not removed

Redundant Internal Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Comprising bond stretches, (linear) bends, (improper) torsions and other primitives
* Less coupling, compared to Cartesians
* Easy estimation of diagonal model Hessians
* **Support constraints**
* Require sophisticated setup algorithm
* Internal-Cartesian backtransformation required, to transform steps from internals to Cartesians
* Redundant set

Delocalized Internal Coordinates (DLC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Complicated linear combinations of primitive internals
* No coupling, coordinates are orthogonal to each other, at least at the geometrey they were defined at
* Non redundant set
* More efficient compared to RIC for bigger systems, if initial diagonalization is feasible
* Same comments apply as for RICs

Supported File Formats
----------------------

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
^^^^^^^^^^^^^^^^

.. code:: text

    C
    N 1 1.35
    H 1 1.0 2 105.0
    H 2 1.4 1 105.0 3 150.0
    H 2 1.4 1 110.0 3 -160.0

Indexing in Z-matrices is 1-based, so the first atom has index 1. Variable substitution
is not supported.

YAML Input
----------

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
      ...                     # no 'fn' key here!
    endopt:
     geom:                    # geom block in endopt takes same keywords as above
      ...                     # no 'fn' key here!

Employed coordinates and coordinate systems in `preopt` and `endopt` are similary
controlled by a `geom` block in the respective block. The same keywords are supported,
as for the global `geom` block, except the `fn` key.
block.

Inline input
^^^^^^^^^^^^^
Inline input is supported for XYZ format, coordinates are expected in Å. Take care
of proper indentation.

.. code:: yaml
    
    geom:
     type: redund
     fn: |
      2

      H 0.0 0.0 0.0
      H 0.0 0.0 0.7

Types of Primitive Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pysisyphus implements many different types of (primitive) internal coordinates.
Every coordinate is defined by its type and a set of atom indices,
e.g., 2 indices for a bond, 3 indices for a bend and 4 indices for
a dihedral.

Specification of a type is necessary, as there are many
different types of bonds, bends and dihedrals/out-of-plane coordinates.
One can't assume, that a coordinate comprised of 3 atom indices is always a
regular bend, as it may also be a linear bend or a translational coordinate
(TRANSLATION_X, 13), describin the mean Cartesian X coordinate of 3 atoms.

Atom indices start at 0!

.. code:: python

    # Primitive types
    BOND = 0
    AUX_BOND = 1
    HYDROGEN_BOND = 2
    INTERFRAG_BOND = 3
    AUX_INTERFRAG_BOND = 4
    BEND = 5
    LINEAR_BEND = 6
    LINEAR_BEND_COMPLEMENT = 7
    PROPER_DIHEDRAL = 8
    IMPROPER_DIHEDRAL = 9
    OUT_OF_PLANE = 10
    LINEAR_DISPLACEMENT = 11
    LINEAR_DISPLACEMENT_COMPLEMENT = 12
    TRANSLATION_X = 13
    TRANSLATION_Y = 14
    TRANSLATION_Z = 15
    # Rotational coordinates are not yet fully implemented
    #ROTATION_A = 16
    #ROTATION_B = 17
    #ROTATION_C = 18
    CARTESIAN_X = 19
    CARTESIAN_Y = 20
    CARTESIAN_Z = 21

As some of these types are quite unwieldy, several shortcuts are supported,
that can be used in place of the types above.

.. code:: python

    # Additional shortcuts
    # Using Cartesians in the framework of internal coordinates is mainly
    # useful if one wants to constrain certain atoms.
    "X": [PT.CARTESIAN_X],
    "Y": [PT.CARTESIAN_Y],
    "Z": [PT.CARTESIAN_Z],
    "XY": [PT.CARTESIAN_X, PT.CARTESIAN_Y],
    "XZ": [PT.CARTESIAN_X, PT.CARTESIAN_Z],
    "YZ": [PT.CARTESIAN_Y, PT.CARTESIAN_Z],
    "XYZ": [PT.CARTESIAN_X, PT.CARTESIAN_Y, PT.CARTESIAN_Z],
    "ATOM": [PT.CARTESIAN_X, PT.CARTESIAN_Y, PT.CARTESIAN_Z],
    # Primitive aliases
    "B": [PT.BOND],
    "A": [PT.BEND],
    "D": [PT.PROPER_DIHEDRAL],
    "DIHEDRAL": [PT.PROPER_DIHEDRAL],

Define Additional Primitives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pysisyphus tries its best, to automatically come up with a reasonable set
of internal coordinates, but sometimes the algorithm misses an important one.
Especially at transition state guesses, where increased atom
distances are common, bonds may be missed.

In such cases, additional coordinate can be requested explicitly. If additional
coordinates are requested, **a nested list is expected [[coord0], [coord1], ...]**.

.. code:: yaml

    # General structure (list of coordinate lists)
    define_prims: [[PrimType or Shortcut], *[atom indices], ...]

    # Examples

    # Additional bond between atoms 4 and 7 (0-based indexing).
    # All three lines below result in the same bond; the latter two use shortcuts.
    define_prims: [[0, 4, 7]]
    define_prims: [[B, 4, 7]]
    define_prims: [[BOND, 4, 7]]

    # Wrong specification (forgot outer list/brackets):
    define_prims: [0, 4, 7]

    # Also define an additional dihedral, beside the bond
    define_prims: [[0, 4, 7], ["D", 0, 1, 2, 3]]

Constraints
^^^^^^^^^^^
Constraints are currently only supported in conjunction with RIC (`coord_type: redund`).
It is not (yet) possible to modify the value of the specified coordinate via YAML
input; the internal coordinate is constrained at its initial value. The same syntax
as for `define_prims` is used. If the coordinate of the requested constraint is not
already defined, it will be defined subsequently. There is no need to also add the
constrained coordinate to `define_prims`.

.. code:: yaml

    # General structure (nested list of coordinates)
    constrain_prims: [[[PrimType or Shortcut], *[atom indices]], ...]

    # Examples

    # Constrain Cartesian coordinate of atom 0.
    # All three lines result in the same constraint.
    constrain_prims: [[0]]  # Special case: one atom index yields Cartesians of atom
    constrain_prims: [[XYZ, 0]]
    constrain_prims: [[ATOM, 0]]

    # Constrain only Cartesian X and Y component of atom 0.
    constrain_prims: [[XY, 0]]

    # Constraint bond between atoms 4 and 7 (0-based indexing).
    # All three lines below result in the same constraint; the latter two use shortcuts.
    constrain_prims: [[0, 4, 7]]
    constrain_prims: [[B, 4, 7]]
    constrain_prims: [[BOND, 4, 7]]

Constraining the Cartesian coordinates (X, Y and Z) of one atom does not affect
the final energy of an optimization. **But constraining more than one atome does.**

Isotopes
^^^^^^^^
Different isotope masses can be requested. The system works similar to Gaussians system.
An list of pairs is expected, where the first number specifies the atom and the
second number is either an integer or a float. If it is an integer, the isotope
mass closest to this integer is looked up in an internal database. Floats are used as is.

.. code:: yaml

    # General structure (nested list of coordinates)
    isotopes: [[[atom index], [new mass, integer/float], ...]

    # Modify the mass of atom with index 2 (hydrogen in this case)
    # Both lines give identical results (deuterium).
    # In the second line, the mass is given directly.
    isotopes: [[2, 2]]
    isotopes: [[2, 2.014101778]]

Different isotope masses affect calculated frequencies and IRCs. Atoms can be fixed
in IRC calculations by specifying a very high mass.

.. code:: yaml

    # Fix atom 0 in IRC calculation.
    isotopes: [[0, 1e9]]


Related Literature
------------------

1. `The efficient optimization of molecular geometries using redundant internal coordinates <https://doi.org/10.1063/1.1515483>`_
2. `The generation and use of delocalized internal coordinates in geometry optimization <https://doi.org/10.1063/1.471864>`_
3. `The efficient optimization of molecular geometries using redundant internal coordinates <https://doi.org/10.1063/1.1515483>`_
4. `Geometry optimization in redundant internal coordinates <https://doi.org/10.1063/1.462844>`_
5. `Geometry optimization made simple with translation and rotation coordinates <https://doi.org/10.1063/1.4952956>`_



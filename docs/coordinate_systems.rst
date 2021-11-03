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
* :code:`type: cart`

Redundant Internal Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Comprising bond stretches, (linear) bends, (improper) torsions and other primitives
* Less coupling, compared to Cartesians
* Easy estimation of diagonal model Hessians
* **Support constraints**
* Require sophisticated setup algorithm
* Iterative internal-Cartesian backtransformation, which may fail
* Usually highly redundant set
* :code:`type: redund`

Delocalized Internal Coordinates (DLC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Complicated linear combinations of primitive internals
* No coupling, coordinates are orthogonal to each other, at least at the geometrey they were defined at
* Non redundant set
* More efficient compared to RIC for bigger systems (if initial DLC generation is feasible)
* Same comments apply, as for RICs
* :code:`type: dlc`

Translation & Rotation Internal Coordinates (TRIC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Especially suited to optimize solvated/non-covalently bound systems
* Translation and rotation coordinates are assigned to every fragment
* Avoids error-prone assignment of interfragment coordinates
* See `10.1063/1.4952956 <https://doi.org/10.1063/1.4952956>`_ for a full discussion
* By default the B-Matrix is recalculated in every step of the internal-Cartesian
  backtransformation when TRIC is enabled
* :code:`type: tric`

Supported File Formats
----------------------

All formats can be read, at least to a certain extent. Coordinate files are expected to be in Å,
if not otherwise dictated by the respective format.

================ ===== =================================
Suffix           Write   Comment            
================ ===== =================================
.xyz             ✓     Plain XYZ, single geometry.
.trj             ✓     Plain XYZ, multiple geometries.
.pdb             ✓     Protein Data Bank.
.molden          ✗     Restricted to [Geometries] block.
.zmat            ✗     Z-Matrix, see below for an example.
.cjson           ✗     As saved by Avogadro.
.crd             ✗     CHARMM card format
.sdf             ✗     Structure-data file
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
is not supported in the current Z-matrix parser.

YAML Input
----------

See below for an explanation of possible inputs in the `geom` section. Input related
to internal coordinates is given mostly in the `coord_kwargs` subgroup, whereas
atom-related input (`isotops`, `freeze_atoms`) is given one level above. See the
example below:

.. code:: yaml

    geom:
     type: cart                     # Coordinate system (cart/redund/dlc/tric)
     fn: [input]                    # File name or inline input
     union: False                   # Define same set of primitives at multiple geometries
     isotopes: null                 # Specify different isotopes
     freeze_atoms: null             # Freeze Cartesians of certain atoms
     coord_kwargs:                  # Keywords that are passed to the internal coordinate class
      define_prims: null            # Additionally define these primitives
      constrain_prims: null         # Primitive internals to be constrained
      freeze_atoms_exclude: False   # Whether to set up internal coordinates for frozen atoms
    preopt:
     geom:                          # geom block in preopt takes same keywords as above
      ...                           # no 'fn' key here!
    endopt:
     geom:                          # geom block in endopt takes same keywords as above
      ...                           # no 'fn' key here!

Employed coordinates and coordinate systems in `preopt` and `endopt` are similary
controlled by a `geom` block. Same keywords are supported, as for the `geom` block,
at the top level, except the `fn` key.

Inline input
^^^^^^^^^^^^^
Inline coordinates are supported for XYZ format, and expected in Å. Take care
of proper indentation. The example below would yield RIC for the hydrogen molecule
(1 bond).

.. code:: yaml
    
    geom:
     type: redund
     fn: |
      2

      H 0.0 0.0 0.0
      H 0.0 0.0 0.7

Types of Primitive Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pysisyphus implements many different (primitive) internal coordinates.
Every coordinate is defined by its type and a set of atom indices,
e.g., 2 indices for a bond, 3 indices for a bend and 4 indices for
a dihedral.

Specification of a type is necessary, as there are many
different kinds of bonds, bends and dihedrals/out-of-plane.
One can't just assume, that a coordinate comprised of 3 atom indices is always a
regular bend, as it may also be a linear bend or a translational coordinate
(TRANSLATION_X, 14), describin the mean Cartesian X coordinate of 3 atoms.

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
    # TRANSLATION = 13  # Dummy coordinate
    TRANSLATION_X = 14
    TRANSLATION_Y = 15
    TRANSLATION_Z = 16
    # ROTATION = 17  # Dummy coordinate
    ROTATION_A = 18
    ROTATION_B = 19
    ROTATION_C = 20
    # CARTESIAN = 21  # Dummy coordinate
    CARTESIAN_X = 22
    CARTESIAN_Y = 23
    CARTESIAN_Z = 24
    BONDED_FRAGMENT = 25
    DUMMY_TORSION = 26
    DISTANCE_FUNCTION = 27
    # atan2 based coordinates
    BEND2 = 28
    TORSION2 = 29

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
    "A2": [PT.BEND2],
    "D": [PT.PROPER_DIHEDRAL],
    "DIHEDRAL": [PT.PROPER_DIHEDRAL],
    "TORSION": [PT.PROPER_DIHEDRAL],
    "D2": [PT.PROPER_DIHEDRAL2],
    "DIHEDRAL2": [PT.PROPER_DIHEDRAL2],
    "TORSION2": [PT.PROPER_DIHEDRAL2],
    # Translation & Rotation coordinates
    "TRANSLATION": [PT.TRANSLATION_X, PT.TRANSLATION_Y, PT.TRANSLATION_Z],
    "ROTATION": [PT.ROTATION_A, PT.ROTATION_B, PT.ROTATION_C],
    # Miscellaneous
    "DIST_FUNC": [PT.DISTANCE_FUNCTION]

Define Additional Primitives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pysisyphus tries its best, to automatically come up with a reasonable set
of internal coordinates, but sometimes the algorithm misses an important one.
Especially at transition state guesses, where increased atom
distances are common, bonds may be missed.

In such cases, additional coordinates can be requested explicitly. If additional
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

Freeze Atoms
^^^^^^^^^^^^

All three Cartesian coordinates (X, Y, Z) of certain atoms can be frozen, so
they always remain at their initial value. By setting `freeze_atoms_exclude`
in `coord_kwargs`, frozen atoms can be excluded from the internal coordinate setup.
By default frozen atoms are included.

.. code:: yaml

    geom:
     type: [type]
     fn: [fn]
     freeze_atoms: [*atom indices]
     coord_kwargs:
      freeze_atoms_exclude: False

    # Example; fully freeze Cartesians of first and second atom.
    geom:
     type: cart
     fn: input.xyz
     freeze_atoms: [0, 1]
     coord_kwargs:
      freeze_atoms: True

Constraints
^^^^^^^^^^^
**Constraints beyond frozen atoms are currently only supported in conjunction with
RIC (`type: redund`).**
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
    # Both lines result in the same constraint.
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

Harmonic restraints to selected primitive internals can be specified in the `calc:`
section (see the :ref:`Restraint` documentation).

Isotopes
^^^^^^^^
Different isotope masses can be requested. The system works similar to Gaussians system.
A list of pairs is expected, where the first number specifies the atom and the
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
in IRC calculations by specifying a very high mass. **YAML does not recognize 1e9 as
float**, take care to add a dot (**1.e9**).

.. code:: yaml

    # Fix atom 0 in IRC calculation.
    isotopes: [[0, 1.e9]]

Geometry & RedundantCoords
--------------------------

The central class in pysisyphus, handling coordinates and delegating calculations
to (external QC) codes, is the `Geometry` class, similar to ASE's `Atoms`.

.. automodule:: pysisyphus.Geometry
    :members:
    :undoc-members:

Basic infrastructure to work with internal coordinates is provided by `RedundantCoords`.

.. automodule:: pysisyphus.intcoords.RedundantCoords
    :members:
    :undoc-members:

Related Literature
------------------

1. `The efficient optimization of molecular geometries using redundant internal coordinates <https://doi.org/10.1063/1.1515483>`_
2. `The generation and use of delocalized internal coordinates in geometry optimization <https://doi.org/10.1063/1.471864>`_
3. `Geometry optimization in redundant internal coordinates <https://doi.org/10.1063/1.462844>`_
4. `Geometry optimization made simple with translation and rotation coordinates <https://doi.org/10.1063/1.4952956>`_



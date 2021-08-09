import itertools as it

from pysisyphus.helpers_pure import OrderedEnum
from pysisyphus.intcoords import (
    Bend,
    CartesianX,
    CartesianY,
    CartesianZ,
    LinearBend,
    LinearDisplacement,
    OutOfPlane,
    RotationA,
    RotationB,
    RotationC,
    Stretch,
    TranslationX,
    TranslationY,
    TranslationZ,
    Torsion,
)


class PrimTypes(OrderedEnum):
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
    TRANSLATION = 13
    TRANSLATION_X = 14
    TRANSLATION_Y = 15
    TRANSLATION_Z = 16
    ROTATION = 17
    ROTATION_A = 18
    ROTATION_B = 19
    ROTATION_C = 20
    CARTESIAN = 21
    CARTESIAN_X = 22
    CARTESIAN_Y = 23
    CARTESIAN_Z = 24


# Alias for easier access
PT = PrimTypes


PrimTypeShortcuts = {
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
    "TORSION": [PT.PROPER_DIHEDRAL],
    # Translation & Rotation coordinates
    "TRANSLATION": [PT.TRANSLATION_X, PT.TRANSLATION_Y, PT.TRANSLATION_Z],
    "ROTATION": [PT.ROTATION_A, PT.ROTATION_B, PT.ROTATION_C],
}

# The tuples below can be used to decide whether a given type belongs
# to a certain class of primitive.
Bonds = (
    PT.BOND,
    PT.AUX_BOND,
    PT.HYDROGEN_BOND,
    PT.INTERFRAG_BOND,
    PT.AUX_INTERFRAG_BOND,
)
Bends = (PT.BEND, )
LinearBends = (
    PT.LINEAR_BEND,
    PT.LINEAR_BEND_COMPLEMENT,
    PT.LINEAR_DISPLACEMENT,
    PT.LINEAR_DISPLACEMENT_COMPLEMENT,
)
Dihedrals = (PT.PROPER_DIHEDRAL, PT.IMPROPER_DIHEDRAL)
OutOfPlanes = (PT.OUT_OF_PLANE, )
Cartesians = (PT.CARTESIAN_X, PT.CARTESIAN_Y, PT.CARTESIAN_Z)
Rotations = (PT.ROTATION_A, PT.ROTATION_B, PT.ROTATION_C)
Translations = (PT.TRANSLATION_X, PT.TRANSLATION_Y, PT.TRANSLATION_Z)


def get_rot_coord(cls):
    def func(indices, ref_coords3d):
        return cls(indices, ref_coords3d=ref_coords3d)

    return func


# Maps primitive types to their classes
PrimMap = {
    PT.BOND: Stretch,
    PT.AUX_BOND: Stretch,
    PT.HYDROGEN_BOND: Stretch,
    PT.INTERFRAG_BOND: Stretch,
    PT.AUX_INTERFRAG_BOND: Stretch,
    PT.BEND: Bend,
    PT.LINEAR_BEND: LinearBend,
    PT.LINEAR_BEND_COMPLEMENT: lambda indices: LinearBend(indices, complement=True),
    PT.PROPER_DIHEDRAL: lambda indices: Torsion(indices, periodic=True),
    PT.IMPROPER_DIHEDRAL: lambda indices: Torsion(indices, periodic=True),
    PT.OUT_OF_PLANE: OutOfPlane,
    PT.LINEAR_DISPLACEMENT: LinearDisplacement,
    PT.LINEAR_DISPLACEMENT_COMPLEMENT: lambda indices: LinearDisplacement(
        indices, complement=True
    ),
    PT.TRANSLATION_X: TranslationX,
    PT.TRANSLATION_Y: TranslationY,
    PT.TRANSLATION_Z: TranslationZ,
    PT.ROTATION_A: get_rot_coord(RotationA),
    PT.ROTATION_B: get_rot_coord(RotationB),
    PT.ROTATION_C: get_rot_coord(RotationC),
    PT.CARTESIAN_X: CartesianX,
    PT.CARTESIAN_Y: CartesianY,
    PT.CARTESIAN_Z: CartesianZ,
}


def normalize_prim_input(prim_inp):
    """Normalize input for define_prims and constrain_prims

    The intcoords.RedundantCoords constructor expects lists of integer lists
    (tuples) for arguments like 'define_prims' and 'constrain_prims'. The first item
    of every list determines the type of primitive coordinate. Currently
    there are about 20 different types and it is hard to remember all of
    them.

    So we also allow a more human friendly input, that is normalized here.
    The most common primitives are:

    0: BOND
    5: BEND
    8: PROPER_DIHEDRAL

    This function maps inputs like ["BOND", 1, 2] to [PrimTypes.BOND, 1, 2] etc.

    Always returns a list of tuples, as some prim_inps expand to multiple
    coordinates, e.g., XYZ or ATOM.
    """
    if prim_inp is None:
        return []

    prim_type, *indices = prim_inp

    # Nothing to do
    if isinstance(prim_type, PrimTypes):
        return [prim_inp]

    # First check if we got something like an integer
    try:
        return [tuple([PrimTypes(int(prim_type))] + indices)]
    # Raised when prim_type is, e.g., "BOND"
    except ValueError:
        pass

    # Check if we got a PrimType name
    try:
        prim_type_ = getattr(PrimTypes, str(prim_type).upper())
        return [tuple([prim_type_] + indices)]
    except AttributeError:
        pass

    # Check if we got a shortcut, e.g, X/Y/Z/XYZ/ATOM etc.
    try:
        prim_types_ = PrimTypeShortcuts[str(prim_type).upper()]
        return [tuple([prim_type_] + indices) for prim_type_ in prim_types_]
    except KeyError as error:
        print(f"Could not normalize 'prim_inp'={prim_inp}!")
        raise error


def normalize_prim_inputs(prim_inps):
    # Flatten list of tuples
    return list(it.chain(*[normalize_prim_input(pi) for pi in prim_inps]))


def prims_from_prim_inputs(prim_inps):
    norm_prim_inps = normalize_prim_inputs(prim_inps)
    prims = [PrimMap[prim_type](indices) for prim_type, *indices in norm_prim_inps]
    return prims

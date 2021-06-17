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
    TRANSLATION_X = 13
    TRANSLATION_Y = 14
    TRANSLATION_Z = 15
    ROTATION_A = 16
    ROTATION_B = 17
    ROTATION_C = 18
    CARTESIAN_X = 19
    CARTESIAN_Y = 20
    CARTESIAN_Z = 21


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
}


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
    # PT.ROTATION_A: lambda indices, ref_coords3d: RotationA(indices, ref_coords3d=ref_coords3d),
    # PT.ROTATION_B: lambda indices, ref_coords3d: RotationB(indices, ref_coords3d),
    # PT.ROTATION_C: lambda indices, ref_coords3d: RotationC(indices, ref_coords3d),
    PT.ROTATION_A: get_rot_coord(RotationA),
    PT.ROTATION_B: get_rot_coord(RotationB),
    PT.ROTATION_C: get_rot_coord(RotationC),
    PT.CARTESIAN_X: CartesianX,
    PT.CARTESIAN_Y: CartesianY,
    PT.CARTESIAN_Z: CartesianZ,
}

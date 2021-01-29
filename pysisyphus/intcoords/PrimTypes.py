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


# Maps primitive types to their classes
PrimMap = {
    PrimTypes.BOND: Stretch,
    PrimTypes.AUX_BOND: Stretch,
    PrimTypes.HYDROGEN_BOND: Stretch,
    PrimTypes.INTERFRAG_BOND: Stretch,
    PrimTypes.AUX_INTERFRAG_BOND: Stretch,
    PrimTypes.BEND: Bend,
    PrimTypes.LINEAR_BEND: LinearBend,
    PrimTypes.LINEAR_BEND_COMPLEMENT: lambda indices: LinearBend(
        indices, complement=True
    ),
    PrimTypes.PROPER_DIHEDRAL: lambda indices: Torsion(indices, periodic=True),
    PrimTypes.IMPROPER_DIHEDRAL: lambda indices: Torsion(indices, periodic=True),
    PrimTypes.OUT_OF_PLANE: OutOfPlane,
    PrimTypes.LINEAR_DISPLACEMENT: LinearDisplacement,
    PrimTypes.LINEAR_DISPLACEMENT_COMPLEMENT: lambda indices: LinearDisplacement(
        indices, complement=True
    ),
    PrimTypes.TRANSLATION_X: TranslationX,
    PrimTypes.TRANSLATION_Y: TranslationY,
    PrimTypes.TRANSLATION_Z: TranslationZ,
    PrimTypes.ROTATION_A: RotationA,
    PrimTypes.ROTATION_B: RotationB,
    PrimTypes.ROTATION_C: RotationC,
    PrimTypes.CARTESIAN_X: CartesianX,
    PrimTypes.CARTESIAN_Y: CartesianY,
    PrimTypes.CARTESIAN_Z: CartesianZ,
}

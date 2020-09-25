from pysisyphus.helpers_pure import OrderedEnum
from pysisyphus.intcoords import Stretch, Bend, LinearBend, Torsion


class PrimTypes(OrderedEnum):
    BOND = 0
    HYDROGEN_BOND = 1
    INTERFRAG_BOND = 2
    AUX_INTERFRAG_BOND = 3
    BEND = 4
    LINEAR_BEND = 5
    LINEAR_BEND_COMPLEMENT = 6
    PROPER_DIHEDRAL = 7
    IMPROPER_DIHEDRAL = 8


# Maps primitive types to their classes
PrimMap = {
    PrimTypes.BOND: Stretch,
    PrimTypes.HYDROGEN_BOND: Stretch,
    PrimTypes.INTERFRAG_BOND: Stretch,
    PrimTypes.AUX_INTERFRAG_BOND: Stretch,
    PrimTypes.BEND: Bend,
    PrimTypes.LINEAR_BEND: LinearBend,
    PrimTypes.LINEAR_BEND_COMPLEMENT: lambda indices: LinearBend(indices),
    PrimTypes.PROPER_DIHEDRAL: lambda indices: Torsion(indices, periodic=True),
    PrimTypes.IMPROPER_DIHEDRAL: lambda indices: Torsion(indices, periodic=True),
}

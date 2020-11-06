from pysisyphus.helpers import geom_loader
from pysisyphus.calculators import DFTBp


def test_dftbp_forces():
    geom = geom_loader("lib:h2o.xyz")
    calc = DFTBp()
    geom.set_calculator(calc)

    forces = geom.forces
    energy = geom.energy

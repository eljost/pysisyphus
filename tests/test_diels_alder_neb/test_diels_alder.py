import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.testing import using

@pytest.mark.skip
@pytest.mark.parametrize(
    "opt_cls",
    [
        pytest.param(QuickMin, marks=using("pyscf")),
        pytest.param(LBFGS, marks=using("pyscf")),
    ]
)
def test_diels_alder_neb(opt_cls):
    geoms = geom_loader("diels_alder_interpolated.trj")
    for i, geom in enumerate(geoms):
        calc_kwargs = {
            "basis": "sto3g",
            "pal": 2,
            "calc_number": i,
        }
        calc = PySCF(**calc_kwargs)
        geom.set_calculator(calc)
    neb = NEB(geoms)

    opt_kwargs = {
        "dump": True,
        "align": True,
        "max_cycles": 15,
    }
    opt = opt_cls(neb, **opt_kwargs)
    opt.run()

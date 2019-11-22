import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.Psi4 import Psi4
from pysisyphus.calculators.PySCF import PySCF


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    [
        pytest.param(Gaussian16, {"route": "HF/STO-3G"}),
        pytest.param(ORCA, {"keywords": "HF STO-3G tightscf"}),
        pytest.param(Psi4, {"method": "scf", "basis": "sto-3g",
                            "to_set": {"scf_type": "pk"}}
        ),
        pytest.param(PySCF, {"basis": "sto-3g"}),
])
def test_cytosin_opt(calc_cls, calc_kwargs):
    geom = geom_from_library("cytosin.xyz", coord_type="redund")
    calc = calc_cls(**calc_kwargs, mem=1000)
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau_tight",
        "overachieve_factor": 2.,
        # "trust_radius": 0.3,
        # "trust_max": 0.3,
        "line_search": True,
        "gdiis": True,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    # gau tight
    assert geom.energy == pytest.approx(-387.54925361)

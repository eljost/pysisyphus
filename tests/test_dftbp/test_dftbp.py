import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.calculators import DFTBp
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@using("dftbp")
@pytest.mark.parametrize(
    "slakos", [
        (None),
        ("slakos"),
    ]
)
def test_dftbp_opt(slakos, this_dir):
    geom = geom_loader("lib:h2o.xyz", coord_type="redund")
    if slakos:
        slakos = this_dir / slakos
    calc_kwargs = {
        "slakos": slakos,
        "parameter": "mio-ext",
    }

    calc = DFTBp(**calc_kwargs)
    geom.set_calculator(calc)

    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-4.07793793)


def test_dftb_energy():
    geom = geom_loader("lib:h2o.xyz", coord_type="redund")
    calc = DFTBp(parameter="mio-ext")
    geom.set_calculator(calc)

    assert geom.energy == pytest.approx(-4.0777507524)

import numpy as np
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


@using("dftbp")
def test_dftb_energy():
    geom = geom_loader("lib:h2o.xyz")
    calc = DFTBp(parameter="mio-ext")
    geom.set_calculator(calc)

    assert geom.energy == pytest.approx(-4.0777507524)


@using("dftbp")
def test_missing_parameters():
    with pytest.raises(AssertionError):
        calc = DFTBp(parameter="these-are-missing")


@using("dftbp")
def test_charge():
    geom = geom_loader("lib:h2o.xyz")
    calc = DFTBp(parameter="mio-ext", charge=2)
    geom.set_calculator(calc)

    assert geom.energy == pytest.approx(-2.6298838274)


@using("dftbp")
def test_triplet():
    geom = geom_loader("lib:h2o.xyz")
    with pytest.raises(AssertionError):
        calc = DFTBp(parameter="mio-ext", mult=3)


@using("dftbp")
@pytest.mark.parametrize(
    "root, ref_energy, ref_norm", [
        (None, -19.203183051, 0.078868392),
        (1, -19.203183051, 0.38002410032),
    ]
)
def test_cytosin_es_forces(root, ref_energy, ref_norm):
    geom = geom_loader("lib:cytosin.xyz")
    calc_kwargs = {
        "parameter": "mio-ext",
        "root": root,
        "track": bool(root),
    }
    calc = DFTBp(**calc_kwargs)
    geom.set_calculator(calc)

    forces = geom.forces
    energy = geom.energy
    norm = np.linalg.norm(forces)
    assert norm == pytest.approx(ref_norm)
    print("energy", energy)


@using("dftbp")
def test_cytosin_s1_opt():
    geom = geom_loader("lib:cytosin.xyz", coord_type="redund")
    calc_kwargs = {
        "parameter": "mio-ext",
        "root": 1,
        "track": True,
    }
    calc = DFTBp(**calc_kwargs)
    geom.set_calculator(calc)

    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-19.13130711)

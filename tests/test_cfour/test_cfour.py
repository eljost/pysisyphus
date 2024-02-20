import pytest

from pysisyphus.calculators import CFOUR
from pysisyphus.helpers import geom_loader
from pysisyphus.irc import EulerPC
from pysisyphus.testing import using


@pytest.fixture
def water():
    geom = geom_loader("lib:h2o.xyz", coord_type="redund")
    cfour_input = {
        "calc": "SCF",
        "basis": "STO-3G",
    }
    calc = CFOUR(cfour_input, wavefunction_dump=False)
    geom.set_calculator(calc)
    return geom


@pytest.fixture
def ch2oo_ts(this_dir):
    geom = geom_loader(this_dir / "ch2oo_ts_bohr.xyz", coord_type="cart")
    cfour_input = {
        "calc": "SCF",
        "basis": "PVDZ",
        "scf_conv": 12,
        "lineq_conv": 12,
    }
    calc = CFOUR(cfour_input, wavefunction_dump=False)
    geom.set_calculator(calc)
    return geom


@using("cfour")
def test_energy(water):
    energy = water.energy
    assert energy == pytest.approx(-74.96070)


@using("cfour")
def test_irc(ch2oo_ts, this_dir):
    irc = EulerPC(
        ch2oo_ts,
        hessian_init=(this_dir / "FCMFINAL"),
        corr_func="scipy",
        rms_grad_thresh=1e-2,
    )
    irc.run()
    assert irc.forward_is_converged
    assert irc.backward_is_converged

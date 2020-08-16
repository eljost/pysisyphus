import pytest

import pysisyphus
from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "cdds, fallback, no_mwfn, no_jmol", [
        (None    , None  , False, False),
        ("calc"  , None  , True , False),
        ("render", "calc", False, True),
        ("render", None  , True , True),
    ]
)
def test_cdds_fallback(cdds, fallback, no_mwfn, no_jmol, monkeypatch):
    # Disable Mwfn/Jmol as requested
    def mock_available(calculator):
        if (calculator == "mwfn") and no_mwfn:
            return False
        if (calculator == "jmol") and no_jmol:
            return False
        return True
    monkeypatch.setattr(pysisyphus.calculators.OverlapCalculator,
                        "available", mock_available)

    calc_kwargs = {
        "cdds": cdds,
    }
    calc = OverlapCalculator(**calc_kwargs)

    assert calc.cdds == fallback


@pytest.fixture
def water():
    geom = geom_loader("lib:h2o.xyz")
    init_logging()
    calc_kwargs = {
        "xc": "pbe0",
        "method": "tddft",
        "basis": "sto3g",
        "nstates": 2,
        "root": 1,
        "track": True,
        "cdds": "calc",
    }
    calc = PySCF(**calc_kwargs)
    geom.set_calculator(calc)
    return geom


@using("pyscf")
def test_mwfn_crash_fallback(water, monkeypatch):
    calc = water.calculator
    calc.cdds = "calc"

    # Mock method to ensure the CDD calculation always crashes.
    def mock_calc_cdd_cube(*args):
        raise Exception("Mocked Multiwfn crash!")
    monkeypatch.setattr(OverlapCalculator, "calc_cdd_cube", mock_calc_cdd_cube)

    energy = water.energy
    # Force recalculation
    water.clear()
    energy_ = water.energy

    # Check that CDD calculation was disabled, after calc_cdds_crashed
    assert calc.cdds == None

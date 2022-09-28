import h5py
import numpy as np
import pytest

import pysisyphus
from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import describe
from pysisyphus.init_logging import init_logging
from pysisyphus.testing import using
from pysisyphus.wavefunction.excited_states import norm_ci_coeffs, tden_overlaps


@pytest.mark.parametrize(
    "cdds, fallback, mwfn, jmol",
    [
        (None, None, True, True),
        ("calc", None, False, True),
        ("render", "calc", True, False),
        ("render", None, False, False),
    ],
)
def test_cdds_fallback(cdds, fallback, mwfn, jmol, monkeypatch):
    # Disable Mwfn/Jmol as requested
    def mock_get_cmd(calculator):
        return {
            "jmol": jmol,
            "mwfn": mwfn,
        }[calculator]

    monkeypatch.setattr(
        pysisyphus.calculators.OverlapCalculator, "get_cmd", mock_get_cmd
    )

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
        # OverlapCalculator specific
        "track": True,
        "cdds": "calc",
        "ovlp_type": "tden",
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


@pytest.mark.parametrize(
    "h5_fn",
    [
        "cytosin_orca_overlap_data.h5",
        "cytosin_trip_orca_overlap_data.h5",
    ],
)
def test_tden_self_overlap(h5_fn, this_dir):
    with h5py.File(this_dir / h5_fn, "r") as handle:
        mo_coeffs = handle["mo_coeffs"][:]
        ci_coeffs = handle["ci_coeffs"][:]

    calc = OverlapCalculator()

    def tden_self_overlap(C, Xa):
        S_AO = calc.get_sao_from_mo_coeffs(C)
        Xa, _ = norm_ci_coeffs(Xa, np.zeros_like(Xa), restricted_norm=1.0)
        return tden_overlaps(C, Xa, C, Xa, S_AO)

    for Ca, Xa in zip(mo_coeffs, ci_coeffs):
        ovlps = tden_self_overlap(Ca.T, Xa)
        I = np.eye(ovlps.shape[0])
        np.testing.assert_allclose(ovlps, I, atol=5e-5)

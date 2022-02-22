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

    def tden_self_overlap(mo_coeffs, ci_coeffs):
        ao_ovlp = calc.get_sao_from_mo_coeffs(mo_coeffs)

        ci_norm = np.linalg.norm(ci_coeffs, axis=(1, 2))
        mo_norm = calc.get_mo_norms(mo_coeffs, ao_ovlp)

        ci_coeffs = ci_coeffs / ci_norm[:, None, None]
        mo_coeffs = calc.renorm_mos(mo_coeffs, ao_ovlp)

        ci_norm = np.linalg.norm(ci_coeffs, axis=(1, 2))
        mo_norm = calc.get_mo_norms(mo_coeffs, ao_ovlp)
        # print(f"norm(CI): {ci_norm}")
        # print(f"norm(MOs): {describe(mo_norm)}")

        overlaps = calc.tden_overlaps(
            mo_coeffs, ci_coeffs, mo_coeffs, ci_coeffs, ao_ovlp
        )
        return overlaps

    for i, (moc, cic) in enumerate(zip(mo_coeffs, ci_coeffs)):
        ovlps = tden_self_overlap(moc, cic)
        I = np.eye(ovlps.shape[0])
        np.testing.assert_allclose(ovlps, I, atol=5e-5)

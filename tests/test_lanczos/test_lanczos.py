import logging
import sys

import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.modefollow import geom_lanczos
from pysisyphus.testing import using


logger = logging.getLogger("lanczos")
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)


def test_anapot_lanczos():
    geom = AnaPot.get_geom((-0.5767, 1.6810, 0.0))
    guess = (0.4, 0.3, 0.0)
    w_min, v_min = geom_lanczos(geom, guess=guess, dx=1e-5, dl=1e-5, logger=logger)

    H = geom.hessian
    w, v = np.linalg.eigh(H)
    w_ref = w[0]
    v_ref = v[:, 0]
    assert w_min == pytest.approx(w_ref, abs=1e-4)
    assert any([np.allclose(v, v_ref, atol=5e-2) for v in (v_min, -v_min)])


@using("pyscf")
@pytest.mark.parametrize(
    "coord_type, guess",
    [
        ("cart", (0.4, 0.0, 0.0, 0.0, -0.3, 0.0, 0.0, 0.0, 0.0)),
        ("redund", (-0.1, 0.1, 0.2, 0.4, -0.3, 0.1)),
        pytest.param(
            "dlc", (0.2, 0.6, -0.1), marks=pytest.mark.skip
        ),  # test seems really flaky in the CI
    ],
)
def test_hcn_iso_lanczos(coord_type, guess):
    geom = geom_loader("lib:hcn_iso_hf_sto3g_ts_opt.xyz", coord_type=coord_type)
    calc = PySCF(pal=2, basis="sto3g")
    geom.set_calculator(calc)

    w_min, v_min = geom_lanczos(geom, guess=guess, logger=logger)

    # Reference values
    H = geom.hessian
    w, v = np.linalg.eigh(H)
    w_ref = w[0]
    v_ref = v[:, 0]

    assert w_min == pytest.approx(w_ref, abs=1e-2)
    assert any([np.allclose(v, v_ref, atol=5e-2) for v in (v_min, -v_min)])

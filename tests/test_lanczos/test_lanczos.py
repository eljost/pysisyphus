import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.modefollow import geom_lanczos
from pysisyphus.testing import using


def test_anapot_lanczos():
    geom = AnaPot.get_geom((-0.5767, 1.6810, 0.))
    w_min, v_min = geom_lanczos(geom, dx=1e-5)

    H = geom.hessian
    w, v = np.linalg.eigh(H)
    w_ref = w[0]
    v_ref = v[:,0]
    assert w_min == pytest.approx(w_ref, abs=1e-3)
    assert any([np.allclose(v, v_ref, atol=5e-2) for v in (v_min, -v_min)])


@using("pyscf")
def test_hcn_iso_lanczos():
    geom = geom_loader("lib:hcn_iso_hf_sto3g_ts_opt.xyz")
    calc = PySCF(pal=2, basis="sto3g")
    geom.set_calculator(calc)

    r_prev = np.zeros_like(geom.coords)
    r_prev[0] = 0.4
    r_prev[4] = -0.3
    w_min, v_min = geom_lanczos(geom, r_prev=r_prev)

    # Reference values
    H = geom.hessian
    w, v = np.linalg.eigh(H)
    w_ref = w[0]
    v_ref = v[:,0]

    assert w_min == pytest.approx(w_ref, abs=1e-2)
    assert any([np.allclose(v, v_ref, atol=5e-2) for v in (v_min, -v_min)])

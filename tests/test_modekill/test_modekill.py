import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.irc import ModeKill
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.testing import using


@using("xtb")
def test_modekill_xtb(this_dir):
    fn = this_dir / "shaked.geom_000.xyz"
    geom = geom_loader(fn)
    calc = XTB(pal=4)
    geom.set_calculator(calc)
    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)
    np.testing.assert_allclose(nus[[0, 1]], (-199.6709792, -94.054831))

    modekill = ModeKill(geom, kill_inds=[0, ])
    modekill.run()

    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)

    assert nus[0] == pytest.approx(-96.20299777)


@using("pyscf")
def test_modekill_pyscf(this_dir):
    fn = this_dir / "ethane_shaked.xyz"
    geom = geom_loader(fn)
    calc = PySCF(basis="sto3g", xc="bp86", pal=2)
    geom.set_calculator(calc)

    w, v = np.linalg.eigh(geom.eckart_projection(geom.mw_hessian))
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(-266.2978391005)

    modekill = ModeKill(geom, kill_inds=[0, ])
    modekill.run()
    assert modekill.converged

    w, v = np.linalg.eigh(geom.eckart_projection(geom.mw_hessian))
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(324.4358155, abs=1e-4)

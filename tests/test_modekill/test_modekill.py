import numpy as np
import pytest

from pysisyphus.calculators import Gaussian16, XTB
from pysisyphus.irc import ModeKill
from pysisyphus.helpers import geom_loader, eigval_to_wavenumber
from pysisyphus.testing import using


@using("xtb")
def test_modekill_xtb(this_dir):
    fn = this_dir / "shaked.geom_000.xyz"
    geom = geom_loader(fn)
    calc = XTB(pal=4)
    geom.set_calculator(calc)
    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(-199.5586146) 
    assert nus[1] == pytest.approx(-94.01553723) 

    modekill = ModeKill(geom, kill_inds=[0, ])
    modekill.run()

    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)

    assert nus[0] == pytest.approx(-96.31027310)


@using("gaussian16")
def test_modekill_g16(this_dir):
    fn = this_dir / "ethane_shaked.xyz"
    geom = geom_loader(fn)
    calc = Gaussian16(route="bp86 sto-3g")
    geom.set_calculator(calc)

    w, v = np.linalg.eigh(geom.eckart_projection(geom.mw_hessian))
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(-265.4967) 

    modekill = ModeKill(geom, kill_inds=[0,])
    modekill.run()
    assert modekill.converged

    w, v = np.linalg.eigh(geom.eckart_projection(geom.mw_hessian))
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(-6.475704754831899e-05)

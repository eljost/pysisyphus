#!/usr/bin/env python3

import numpy as np
import pytest

from pysisyphus.calculators import Gaussian16
from pysisyphus.calculators.XTB import XTB
from pysisyphus.irc import ModeKill
from pysisyphus.helpers import geom_from_xyz_file, eigval_to_wavenumber


def test_modekill_xtb(datadir):
    fn = datadir / "shaked.geom_000.xyz"
    geom = geom_from_xyz_file(fn)
    calc = XTB(pal=4)
    geom.set_calculator(calc)
    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(-199.026) 
    assert nus[1] == pytest.approx(-98.2082) 

    modekill = ModeKill(geom, kill_inds=[0, ])
    modekill.run()

    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)

    assert nus[0] == pytest.approx(-98.2557)


def test_modekill_g16(datadir):
    fn = datadir / "ethane_shaked.xyz"
    geom = geom_from_xyz_file(fn)
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

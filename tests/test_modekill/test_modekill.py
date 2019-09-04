#!/usr/bin/env python3

import numpy as np
from pytest import approx

from pysisyphus.calculators.XTB import XTB
from pysisyphus.irc import ModeKill
from pysisyphus.helpers import geom_from_xyz_file, eigval_to_wavenumber


def test_modekill(datadir):
    fn = datadir / "shaked.geom_000.xyz"
    geom = geom_from_xyz_file(fn)
    calc = XTB(pal=4)
    geom.set_calculator(calc)
    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)
    assert nus[0] == approx(-199.026) 
    assert nus[1] == approx(-98.2082) 

    modekill = ModeKill(geom, kill_inds=[0, ])
    modekill.run()

    w, v = np.linalg.eigh(geom.mw_hessian)
    nus = eigval_to_wavenumber(w)

    assert nus[0] == approx(-98.2557)

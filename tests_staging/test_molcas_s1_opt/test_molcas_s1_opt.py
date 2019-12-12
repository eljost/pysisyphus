#!/usr/bin/env python3

import os
from pathlib import Path
import pytest

from pysisyphus.helpers import geom_from_xyz_file
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.calculators.OpenMolcas import OpenMolcas
from pysisyphus.init_logging import init_logging


@pytest.mark.skip
def test_molcas_s1_opt():
    """
    Optimization of the S1 of trans-Butadien
    Spent 23.6 s preparing the first cycle.
    cycle   max(force)    rms(force)     max(step)     rms(step)       s/cycle
        1     0.072748     0.033504     0.040000     0.019036         26.0
        2     0.020289     0.008667     0.040000     0.017887         25.1
        3     0.011788     0.005995     0.016634     0.007007         25.2
        4     0.010659     0.005249     0.005894     0.002998         25.3
        5     0.007637     0.003809     0.020000     0.009861         24.9
    Number of cycles exceeded!
    """
    path = Path(os.path.dirname(os.path.abspath(__file__)))
    xyz_fn = path / "trans_butadien.xyz"
    inporb_fn = path / "butadien_vdzp.RasOrb"

    geom = geom_from_xyz_file(xyz_fn)
    kwargs = {
        "basis": "ano-rcc-vdzp",
        "inporb": inporb_fn,
        "charge": 0,
        "mult": 1,
        "roots": 5,
        "mdrlxroot": 2,
    }
    calc = OpenMolcas(**kwargs)
    geom.set_calculator(calc)
    opt_kwargs = {
        "dump": True,
        "max_cycles": 5,
    }
    opt = BFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.max_forces[-1] == pytest.approx(0.0076367199999)
    assert opt.rms_forces[-1] == pytest.approx(0.0038088424813)


if __name__ == "__main__":
    init_logging()
    test_molcas_s1_opt()

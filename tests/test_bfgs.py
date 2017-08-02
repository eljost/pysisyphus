#!/usr/bin/env python3

import pytest

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.BFGS import BFGS

from qchelper.geometry import parse_xyz_file

@pytest.mark.skip(reason="just implementing this")
def test_bfgs():
    ethan_xyz = "xyz_files/ethan.xyz"
    atoms, coords = parse_xyz_file(ethan_xyz)
    ethan_geom = Geometry(atoms, coords.flatten())
    ethan_geom.set_calculator(ORCA())

    kwargs = {
        "max_cycles":28,
    }
    opt = BFGS(ethan_geom, **kwargs)
    opt.run()

#!/usr/bin/env python3

import pytest

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.FIRE import FIRE

from qchelper.geometry import parse_xyz_file


@pytest.mark.skip(reason="Takes too long.")
def test_ethan():
    ethan_xyz = "xyz_files/ethan.xyz"
    atoms, coords = parse_xyz_file(ethan_xyz)
    ethan_geom = Geometry(atoms, coords.flatten())
    ethan_geom.set_calculator(ORCA())

    kwargs = {
        "max_cycles":28,
    }
    opt = FIRE(ethan_geom, **kwargs)
    opt.run()
    assert(opt.is_converged)

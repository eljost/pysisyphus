#!/usr/bin/env python

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient

import numpy as np
import pytest

np.set_printoptions(suppress=True, precision=4)

def get_geom():
    init_logging("./")
    geom = geom_from_library("h2.xyz")
    kwargs = {
        "keywords": "BP86 def2-SV(P)",
        "charge": 0,
        "mult": 1,
        "blocks": "%tddft nroots 2 iroot 1 end",
        "track": True,
    }
    orca = ORCA(**kwargs)
    geom.set_calculator(orca)

    return geom


@pytest.mark.skip
def test_orca_track():
    geom = get_geom()
    forces = geom.forces
    wfow = geom.calculator.wfow
    res = wfow.compare(wfow)
    np.testing.assert_allclose(res, [0, 1, 2])


@pytest.mark.skip
def test_orca_tracked_opt():
    geom = get_geom()
    opt = ConjugateGradient(geom)
    opt.run()
    assert opt.cur_cycle == 15
    assert opt.converged


if __name__ == "__main__":
    test_orca_track()
    #test_orca_tracked_opt()

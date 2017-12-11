#!/usr/bin/env python3

import os
from pathlib import Path
import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.IDPP import idpp_interpolate
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.calculators.XTB import XTB
from pysisyphus.cos.NEB import NEB


def prepare_opt(idpp=True):
    path = Path(os.path.dirname(os.path.abspath(__file__)))
    fns = ("hcn.xyz", "hcn_iso_ts.xyz", "nhc.xyz")
    geoms = [geom_from_library(fn) for fn in fns]

    calc_kwargs = {
        "charge": 0,
        "mult": 1,
    }
    if idpp:
        geoms = idpp_interpolate(geoms, 5)
        neb = NEB(geoms)
    else:
        neb = NEB(geoms)
        neb.interpolate(5)
    for image in neb.images:
        image.set_calculator(XTB(**calc_kwargs))
    return neb

@pytest.mark.skip
def test_xtb_hcn_iso():
    neb = prepare_opt()
    opt_kwargs = {
        "dump": True,
        "align": True,
        "dont_skip": False,
    }
    #opt = BFGS(neb, **opt_kwargs)
    opt = ConjugateGradient(neb, **opt_kwargs)
    opt.run()

    assert (opt.is_converged)
    assert (opt.cur_cycle == 22)


@pytest.mark.skip
def test_xtb_hcn_climb_iso():
    neb = prepare_opt()
    opt_kwargs = {
        "dump": True,
        "align": True,
        "climb": True,
        "max_cycles": 100,
    }
    opt = ConjugateGradient(neb, **opt_kwargs)
    #opt = SteepestDescent(neb, **opt_kwargs)
    opt.run()

    assert (opt.is_converged)
    assert (opt.cur_cycle == 21)


if __name__ == "__main__":
    test_xtb_hcn_iso()
    ##test_xtb_hcn_climb_iso()

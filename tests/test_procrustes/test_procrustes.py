#!/usr/bin/env python3

from pysisyphus.calculators.IDPP import idpp_interpolate
from pysisyphus.calculators.XTB import XTB
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.SteepestDescent import SteepestDescent

import pytest


def get_geoms():
    fns = ("xtbopt.anfang.xyz", "xtbopt.ende2.xyz")
    geoms = [geom_from_library(fn) for fn in fns]
    geoms = idpp_interpolate(geoms, 4)
    return geoms


def get_cos():
    geoms = get_geoms()
    cos = ChainOfStates(geoms)
    return cos


def get_neb():
    geoms = get_geoms()
    for g in geoms:
        g.set_calculator(XTB())
    neb = NEB(geoms)
    return neb


@pytest.mark.skip
def test_procrustes():
    cos = get_cos()

    """
    # Save original coordinates
    with open("org.trj", "w") as handle:
        handle.write(cos.as_xyz())
    """

    opt = Optimizer(cos)
    opt.procrustes()
    opt.write_cycle_to_file()

    """
    # Save rotated coordinates
    with open("rot.trj", "w") as handle:
        handle.write(cos.as_xyz())
    """


@pytest.mark.skip
def test_fit_rigid():
    neb = get_neb()
    opt = BFGS(neb, align=True)
    opt.run()
    assert(opt.cur_cycle == 33)


if __name__ == "__main__":
    #test_procrustes()
    test_fit_rigid()

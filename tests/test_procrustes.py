#!/usr/bin/env python3

from pysisyphus.calculators.IDPP import idpp_interpolate
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.Optimizer import Optimizer


def test_procrustes():
    fns = ("xtbopt.anfang.xyz", "xtbopt.ende2.xyz")
    geoms = [geom_from_library(fn) for fn in fns]
    geoms = idpp_interpolate(geoms, 4)
    cos = ChainOfStates(geoms)

    """
    # Save original coordinates
    with open("org.trj", "w") as handle:
        handle.write(cos.as_xyz())
    """

    opt = Optimizer(cos)
    opt.procrustes()

    """
    # Save rotated coordinates
    with open("rot.trj", "w") as handle:
        handle.write(cos.as_xyz())
    """


if __name__ == "__main__":
    test_procrustes()

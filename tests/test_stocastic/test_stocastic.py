#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.stocastic.Kick import Kick
from pysisyphus.stocastic.FragmentKick import FragmentKick

import numpy as np


np.set_printoptions(suppress=True, precision=2)


def test_kick():
    geom = geom_from_library("benzene_and_chlorine.xyz")
    kick = Kick(geom, cycle_size=10, radius=1.25, seed=1532002565, cycles=5)
    kick.run()


def test_fragment_kick():
    geom = geom_from_library("benzene_and_chlorine.xyz")
    benz_frag = range(12)
    chlorine_frag = (12, 13)
    fragments = (benz_frag, chlorine_frag)
    kwargs = {
        "cycle_size": 5,
        "radius": 1.25,
        "cycles": 2,
        "seed": 1532002565,
    }
    fkick = FragmentKick(geom, fragments, **kwargs)
    fkick.run()

def test_toluene():
    geom = geom_from_library("toluene_and_cl2.xyz")
    toluene_frag = range(15)
    cl2_frag = (15, 16)
    fragments = (toluene_frag, cl2_frag)
    # kwargs = {
        # "cycle_size": 10,
        # "radius": 1.25,
        # "cycles": 2,
        # "seed": 1532002565,
    # }
    # fkick = FragmentKick(geom, fragments, **kwargs)
    # fkick.run()
    kwargs = {
        "cycle_size": 50,
        "radius": 1.25,
        "cycles": 2,
        "seed": 1532002565,
    }
    fkick = FragmentKick(geom, fragments, **kwargs)
    fkick.run_mod()


if __name__ == "__main__":
    # test_kick()
    # test_fragment_kick()
    test_toluene()

#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import pytest

from pysisyphus.helpers import geom_from_library, geoms_from_trj
from pysisyphus.stocastic.Kick import Kick
from pysisyphus.stocastic.FragmentKick import FragmentKick

np.set_printoptions(suppress=True, precision=4)
THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))

# Input fragments have to be reasonable


@pytest.mark.skip
def test_kick():
    geom = geom_from_library("benzene_and_chlorine.xyz")
    kick = Kick(geom, cycle_size=10, radius=1.25, seed=1532002565, cycles=5)
    kick.run()


@pytest.mark.skip
def test_benz_chlorine_fragment_kick():
    geom = geom_from_library("benzene_and_chlorine.xyz")
    benz_frag = range(12)
    chlorine_frag = (12, 13)
    fragments = (benz_frag, chlorine_frag)
    kwargs = {
        "cycle_size": 30,
        "cycles": 10,
        # "cycle_size": 10,
        # "cycles": 5,
        "radius": 4.5,
        "seed": 1532002565,
        "fragments": fragments,
        "fix_fragments": (0, ),
        "rmsd_thresh": .2,
    }
    fkick = FragmentKick(geom, **kwargs)
    fkick.run()


@pytest.mark.skip
def test_benz_no_plus_fragment_kick():
    geom = geom_from_library("benzene_and_no.xyz")
    benz_frag = range(12)
    no_frag = (12, 13)
    fragments = (benz_frag, no_frag)
    kwargs = {
        "cycle_size": 10,
        "cycles": 3,
        "radius": 4.5,
        "seed": 1532002565,
        "fragments": fragments,
        "fix_fragments": (0, ),
        "rmsd_thresh": .2,
        "calc_kwargs": {
            "charge": 1,
        },
    }
    fkick = FragmentKick(geom, **kwargs)
    fkick.run()


@pytest.mark.skip
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
        "cycle_size": 75,
        "radius": 3,
        "cycles": 3,
        "seed": 1532002565,
    }
    fkick = FragmentKick(geom, fragments, **kwargs)
    fkick.run()


def test_atoms_are_too_close():
    trj_fn = THIS_DIR / "test_reject.trj"
    geoms = geoms_from_trj(trj_fn)
    kick = Kick(geoms[0])
    reject = [kick.atoms_are_too_close(geom, factor=.7) for geom in geoms]
    inds = [i for i, r in enumerate(reject) if r]
    assert inds == [6, 9, 14, 19, 20, 23, 28]


if __name__ == "__main__":
    # test_kick()
    # test_benz_chlorine_fragment_kick()
    test_benz_no_plus_fragment_kick()
    # test_toluene()
    # test_atoms_are_too_close()

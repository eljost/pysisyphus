#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import pytest

from pysisyphus.helpers import geom_loader, geom_from_library, geoms_from_trj
from pysisyphus.stocastic import *
from pysisyphus.testing import using


def test_get_fragments():
    geom = geom_loader("lib:benzene_and_chlorine.xyz")
    # Only define chlorine fragment and determine second fragment automatically
    fragments = ((12, 13), )
    kwargs = {
        "fragments": fragments,
    }
    fkick = FragmentKick(geom, **kwargs)
    _, benz_frag = fkick.fragments
    assert benz_frag.tolist() == list(range(12))


@using("xtb")
def test_kick():
    geom = geom_loader("lib:benzene_and_chlorine.xyz")
    stoc_kwargs = {
        "cycle_size": 10,
        "radius": 1.25,
        "seed": 1532002565,
        "max_cycles": 5,
    }
    stoc = Kick(geom, **stoc_kwargs)
    stoc.run()

    assert stoc.cur_cycle == 4
    assert len(stoc.new_geoms) == 9
    assert min(stoc.new_energies) == pytest.approx(-24.9688182)


@using("xtb")
def test_benz_chlorine_fragment_kick():
    geom = geom_loader("lib:benzene_and_chlorine.xyz")
    fragments = (list(range(12)), )
    stoc_kwargs = {
        "cycle_size": 15,
        "max_cycles": 3,
        "radius": 4.5,
        "seed": 1532002565,
        "fragments": fragments,
        "fix_fragments": (0, ),
        "rmsd_thresh": .2,
    }
    stoc = FragmentKick(geom, **stoc_kwargs)
    stoc.run()

    assert stoc.cur_cycle == 2
    assert len(stoc.new_geoms) == 5
    assert min(stoc.new_energies) == pytest.approx(-24.9911479)


@pytest.mark.skip
def test_benz_no_plus_fragment_kick():
    geom = geom_from_library("benzene_and_no.xyz")
    benz_frag = range(12)
    no_frag = (12, 13)
    fragments = (benz_frag, no_frag)
    kwargs = {
        "cycle_size": 10,
        "cycles": 3,
        "radius": 2.5,
        "seed": 1532002565,
        "fragments": fragments,
        "fix_fragments": (0, ),
        "rmsd_thresh": .2,
        "calc_kwargs": {
            "charge": 1,
        },
        "random_displacement": True,
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
        # "cycle_size": 15,
        "radius": 3,
        "cycles": 15,
        # "cycles": 5,
        "seed": 1532002565,
    }
    fkick = FragmentKick(geom, fragments, **kwargs)
    fkick.run()


def test_atoms_are_too_close():
    THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))
    trj_fn = THIS_DIR / "test_reject.trj"
    geoms = geoms_from_trj(trj_fn)
    kick = Kick(geoms[0])
    reject = [kick.atoms_are_too_close(geom, factor=.7) for geom in geoms]
    inds = [i for i, r in enumerate(reject) if r]
    assert inds == [6, 9, 14, 19, 20, 23, 28]




if __name__ == "__main__":
    # test_kick()
    # test_benz_chlorine_fragment_kick()
    # test_benz_no_plus_fragment_kick()
    test_toluene()
    # test_atoms_are_too_close()
    # test_get_fragments()

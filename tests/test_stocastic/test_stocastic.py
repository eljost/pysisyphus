#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import pytest

from pysisyphus.helpers import (geom_from_library, geoms_from_trj,
                                geom_from_xyz_file)
from pysisyphus.stocastic.Kick import Kick
from pysisyphus.stocastic.FragmentKick import FragmentKick
from pysisyphus.stocastic.align import matched_rmsd

np.set_printoptions(suppress=True, precision=4)
THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


@pytest.mark.skip
def test_kick():
    geom = geom_from_library("benzene_and_chlorine.xyz")
    kick = Kick(geom, cycle_size=10, radius=1.25, seed=1532002565, cycles=5)
    kick.run()


@pytest.mark.skip
def test_fragment_kick():
    geom = geom_from_library("benzene_and_chlorine.xyz")
    benz_frag = range(12)
    chlorine_frag = (12, 13)
    fragments = (benz_frag, chlorine_frag)
    kwargs = {
        "cycle_size": 100,
        "radius": 3.5,  #1.25,
        "cycles": 25,
        "seed": 1532002565,
        "fix_fragments": (0, ),
        "energy_thresh": 1e-4,
    }
    fkick = FragmentKick(geom, fragments, **kwargs)
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


def test_match_rmsd():
    geom1 = geom_from_xyz_file("eins.xyz")
    # Calling with the identical geometries should return RMSD of 0.
    min_rmsd, (geom1_matched, geom2_matched) = matched_rmsd(geom1, geom1)
    np.testing.assert_allclose(min_rmsd, 0.0, atol=1e-10)
    np.testing.assert_allclose(geom1_matched.coords, geom2_matched.coords)

    geom2 = geom_from_xyz_file("zwei.xyz")
    min_rmsd, _ = matched_rmsd(geom1, geom2)
    np.testing.assert_allclose(min_rmsd, 0.057049, atol=1e-5)


if __name__ == "__main__":
    # test_kick()
    # test_fragment_kick()
    # test_toluene()
    # test_atoms_are_too_close()
    test_match_rmsd()

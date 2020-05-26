import numpy as np
import pytest

from pysisyphus.helpers import geom_loader
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


def test_atoms_are_too_close():
    geoms = geom_loader("lib:test_reject.trj")
    stoc = Kick(geoms[0])
    reject = [stoc.atoms_are_too_close(geom, factor=.7) for geom in geoms]
    inds = [i for i, r in enumerate(reject) if r]
    assert inds == [6, 9, 14, 19, 20, 23, 28]


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
    stoc_kwargs = {
        "cycle_size": 15,
        "max_cycles": 3,
        "radius": 4.5,
        "seed": 1532002565,
        "fragments": (list(range(12)), ),
        "rmsd_thresh": .2,
    }
    stoc = FragmentKick(geom, **stoc_kwargs)
    stoc.run()

    assert stoc.cur_cycle == 2
    assert len(stoc.new_geoms) == 5
    assert min(stoc.new_energies) == pytest.approx(-24.9911479)


@using("xtb")
def test_benz_no_plus_fragment_kick():
    geom = geom_loader("lib:benzene_and_no.xyz")
    stoc_kwargs = {
        "cycle_size": 10,
        "max_cycles": 3,
        "radius": 2.5,
        "seed": 1532002565,
        "fragments": (list(range(12)), ),
        "rmsd_thresh": .2,
        "calc_kwargs": {
            "charge": 1,
        },
        "random_displacement": True,
    }

    stoc = FragmentKick(geom, **stoc_kwargs)
    stoc.run()

    assert stoc.cur_cycle == 2
    assert len(stoc.new_geoms) == 17
    assert min(stoc.new_energies) == pytest.approx(-22.37273674)


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

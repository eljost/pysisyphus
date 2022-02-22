import pytest

from pysisyphus.calculators.Rastrigin import Rastrigin
from pysisyphus.calculators.PySCF import PySCF
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


@pytest.mark.skip_ci
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


@pytest.mark.skip_ci
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


@pytest.mark.skip_ci
@using("xtb")
def test_benz_no_plus_fragment_kick():
    geom = geom_loader("lib:benzene_and_no.xyz")
    stoc_kwargs = {
        "cycle_size": 10,
        "max_cycles": 3,
        "radius": 5.,
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
    assert len(stoc.new_geoms) == 11
    assert min(stoc.new_energies) == pytest.approx(-22.3727376)


@pytest.mark.skip_ci
@using("xtb")
def test_toluene():
    geom = geom_loader("lib:toluene_and_cl2.xyz")
    stoc_kwargs = {
        "cycle_size": 10,
        "max_cycles": 3,
        "radius": 3,
        "seed": 1532002565,
        "fragments": (list(range(15)), ),
    }
    stoc = FragmentKick(geom, **stoc_kwargs)
    stoc.run()

    assert stoc.cur_cycle == 2
    assert len(stoc.new_geoms) == 4
    assert min(stoc.new_energies) == pytest.approx(-28.1439149)


def test_rastrigin():
    calc = Rastrigin()

    def calc_getter(calc_number):
        calc.calc_counter = 0
        calc.calc_number = calc_number
        return calc

    geom = calc.get_minima()[0]
    stoc_kwargs = {
        "seed": 20180325,
        "calc_getter": calc_getter,
        "cycle_size": 25,
        "max_cycles": 10,
        "radius": 4,
    }
    stoc = Kick(geom, **stoc_kwargs)
    geoms = stoc.run()

    assert stoc.cur_cycle == 9
    assert len(stoc.new_geoms) == 14
    assert min(stoc.new_energies) == pytest.approx(0., abs=1.5e-8)

    # import numpy as np
    # import matplotlib.pyplot as plt
    # calc.plot()
    # ax = calc.ax
    # xs, ys = np.array([geom.coords[:2] for geom in geoms]).T
    # ax.scatter(xs, ys, s=35, marker="X", color="red", zorder=10)
    # plt.show()


# Takes 30 min in the CI which is way too much for a rarely used feature
@pytest.mark.skip_ci
@using("pyscf")
def test_pyscf_stocastic():
    geom = geom_loader("lib:benzene_and_no.xyz")

    def calc_getter(calc_number):
        calc_kwargs = {
            "charge": +1,
            "mult": 1,
            "pal": 2,
            "basis": "321g",
            "calc_number": calc_number,
        }
        calc = PySCF(**calc_kwargs)
        return calc

    stoc_kwargs = {
        "calc_getter": calc_getter,
        "cycle_size": 3,
        "max_cycles": 3,
        "radius": 5.,
        "seed": 20180325,
        "fragments": (list(range(12)), ),
        "rmsd_thresh": .2,
        "random_displacement": True,
    }

    stoc = FragmentKick(geom, **stoc_kwargs)
    stoc.run()

    assert stoc.cur_cycle == 2
    assert len(stoc.new_geoms) == 4
    assert min(stoc.new_energies) == pytest.approx(-357.605594464)

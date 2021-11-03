import matplotlib.pyplot as plt
import numpy as np

import pytest

from pysisyphus.calculators.FourWellAnaPot import FourWellAnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.cos.GrowingNT import GrowingNT
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


def plot_gnt(calc, gnt):
    calc.plot_geoms(gnt.images)
    calc.plot_geoms(gnt.sp_images)

    # fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    # ens = np.array(gnt.all_energies)
    # xs = np.arange(ens.size)
    # ens -= ens.min()
    # max_ind = ens.argmax()
    # norms = np.array(gnt.all_real_forces)
    # norms = np.linalg.norm(norms, axis=1)
    # ax0.plot(ens, "o-")
    # ax0.set_title("Energies")
    # ax1.plot(norms, "o-")
    # ax1.set_title("norm(forces)")
    # ax1.set_xlabel("Image")
    # for ax in (ax0, ax1):
        # ax.axvline(xs[max_ind], c="k", ls="--")
    plt.show()


@pytest.mark.parametrize(
    "between, r_update",
    [
        (9, False),
        (21, False),
    ],
)
def test_mb_growingnt(between, r_update):
    geoms = MullerBrownPot().get_minima()
    geoms = geoms[1], geoms[0]
    geom0, geom1 = geoms

    gnt_kwargs = {
        "between": between,
        "rms_thresh": 0.08,
        "final_geom": geom1,
        "r_update_thresh": 0.9,
        # "r_update": r_update,
        "r_update": True,
    }
    gnt = GrowingNT(geom0, **gnt_kwargs)
    opt_kwargs = {
        "max_cycles": 100,
    }
    opt = PreconLBFGS(gnt, **opt_kwargs)
    opt.run()

    # plot_gnt(geoms[0].calculator, gnt)

    assert opt.is_converged
    assert gnt.images[-1].energy == pytest.approx(-146.69951721)


def test_four_well_growingnt():
    geoms = FourWellAnaPot().get_minima()
    geoms = geoms[0], geoms[-1]
    geom0, geom1 = geoms

    gnt_kwargs = {
        "between": 23,
        "rms_thresh": 0.005,
        "final_geom": geom1,
    }
    gnt = GrowingNT(geom0, **gnt_kwargs)
    opt_kwargs = {
        "max_cycles": 100,
    }
    opt = PreconLBFGS(gnt, **opt_kwargs)
    opt.run()

    # plot_gnt(geoms[0].calculator, gnt)

    assert opt.is_converged
    assert gnt.images[-1].energy == pytest.approx(-6.76245257)


import pytest
from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@pytest.mark.skip
@using("xtb")
@pytest.mark.parametrize(
    "bonds",
    [
        # [[8, 7, 1], [8, 0, -1], [0, 5, -1], [7, 6, -1]],
        # [[8, 7, 1], [8, 0, -1], [0, 5, -1]],
        # [[18, 3, 1], [18, 15, -0.5]],
        [[18, 3, 1], [18, 15, -1]],
    ],
)
def test_biaryl_growingnt(bonds, this_dir):
    fn = "aligned.trj"
    fn = "left_irc_start_right.trj"
    geom0, geom1 = geom_loader(fn)
    geom0.set_calculator(XTB(pal=6, quiet=True))

    gnt_kwargs = {
        "bonds": bonds,
        # "final_geom": geom1,
        "rms_thresh": 0.0015,
        "step_len": 0.75,
        # "between": 10,
        "update_r": True,
    }
    gnt = GrowingNT(geom0, **gnt_kwargs)
    opt_kwargs = {
        "max_cycles": 1000,
        "dump": True,
    }
    # opt = PreconLBFGS(gnt, **opt_kwargs)
    from pysisyphus.optimizers.LBFGS import LBFGS

    opt = LBFGS(gnt, max_step=0.1, **opt_kwargs)
    opt.run()


@pytest.mark.skip
@using("xtb")
@pytest.mark.parametrize(
    "bonds",
    [
        [[11, 3, 1], [0, 10, 1]],
    ],
)
def test_diels_alder_growingnt(bonds, this_dir):
    geoms = geom_loader(
        "/home/johannes/Code/pysisyphus/pysisyphus/geom_library/birkholz_rx/02_hcn_original.trj"
    )
    geoms = geom_loader("lib:diels_alder_interpolated.trj")
    geom0 = geoms[0]
    geom1 = geoms[-1]
    geom0.set_calculator(XTB(pal=6, quiet=True))

    gnt_kwargs = {
        # "bonds": bonds,
        "final_geom": geom1,
        "between": 20,
        # "step_len": 0.05,
        # "rms_thresh": 0.001,
    }
    gnt = GrowingNT(geom0, **gnt_kwargs)
    opt_kwargs = {
        "max_cycles": 100,
        "dump": True,
        "max_step": 0.1,
        "line_search": False,
    }
    opt = PreconLBFGS(gnt, **opt_kwargs)
    # from pysisyphus.optimizers.LBFGS import LBFGS
    # opt = LBFGS(gnt, max_step=0.1, **opt_kwargs)

    opt.run()


@pytest.mark.skip
@using("xtb")
def test_hcn_growingnt():
    geoms = geom_loader("lib:birkholz_rx/02_hcn_original.trj")
    geom0 = geoms[0]
    geom1 = geoms[-1]
    geom0.set_calculator(XTB(pal=1, quiet=True))

    # from pysisyphus.intcoords import Bend
    # indices = (1, 0, 2)
    # _, r = Bend._calculate(geom0.coords3d, indices, gradient=True)

    gnt_kwargs = {
        # "step_len": 0.2,
        "final_geom": geom1,
        "between": 18,
        # "r": r,
        # "rms_thresh": 0.001,
    }
    gnt = GrowingNT(geom0, **gnt_kwargs)
    opt_kwargs = {
        "max_cycles": 100,
        "dump": True,
        "max_step": 0.1,
        "line_search": False,
    }
    opt = PreconLBFGS(gnt, **opt_kwargs)
    # from pysisyphus.optimizers.LBFGS import LBFGS
    # opt = LBFGS(gnt, **opt_kwargs)

    opt.run()

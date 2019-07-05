#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.StabilizedQNMethod import StabilizedQNMethod
from pysisyphus.calculators.XTB import XTB
from pysisyphus.calculators.Gaussian16 import Gaussian16


def test_sqnm():
    # geom = AnaPot.get_geom((0, 4, 0))
    geom = AnaPot.get_geom((-0.8, 1.73, 0))

    opt_kwargs = {
        "max_cycles": 15,
        "eps": 1e-4,
        "E_thresh": 1e-4,
        "alpha": 0.5,
        "hist_max": 10,
        # "thresh": "gau_tight",
        "trust_radius": 0.1,
        # "bio": False,
        "dump": True,
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)
    # opt = SQNM3(geom, **opt_kwargs)
    opt.run()
    c = np.array(opt.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(c[:,0], c[:,1], "o-")

    plt.show()


def test_sqnm_bio_mode():
    atoms = "h h".split()
    coords = ((0, 0, 0), (0, 0, 1))
    geom = Geometry(atoms, coords)
    xtb = XTB()
    geom.set_calculator(xtb)

    opt = StabilizedQNMethod(geom)
    cur_grad = geom.gradient
    stretch_grad, rem_grad = opt.bio_mode(cur_grad)
    # There is only one bond in H2, so stretch_gradient == cur_grad and the
    # remainder should be the zero vector.
    np.testing.assert_allclose(stretch_grad, cur_grad)
    np.testing.assert_allclose(np.linalg.norm(rem_grad), 0.)
    return


def test_sqnm_xtb():
    geom = geom_from_library("split.image_021.xyz")
    # geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    xtb = XTB()
    geom.set_calculator(xtb)

    opt_kwargs = {
        "max_cycles": 100,
        "hist_max": 10,
        "dump": True,
        "trust_radius": 0.1,
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)

    # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    # opt = RFOptimizer(geom)

    opt.run()


def test_mecoome_sqnm_xtb():
    geom = geom_from_library("mecoome_split.image_010.xyz")
    xtb = XTB()
    geom.set_calculator(xtb)

    opt_kwargs = {
        # "alpha": 0.5,
        # "alpha_stretch": 0.5,
        "alpha": 0.5,
        "alpha_stretch": 0.5,
        "hist_max": 5,
        "dump": True,
        "trust_radius": 0.5,
        "bio": False,
        # "thresh": "gau",
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)

    # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    # opt = RFOptimizer(geom)

    opt.run()


def test_mecoome_sqnm_g16():
    geom = geom_from_library("mecoome_split.image_010.xyz")
    calc = Gaussian16("PM6", pal=4)
    geom.set_calculator(calc)

    opt_kwargs = {
        # "alpha": 0.5,
        # "alpha_stretch": 0.5,
        "alpha": 0.5,
        "alpha_stretch": 0.5,
        "hist_max": 5,
        "dump": True,
        "trust_radius": 0.5,
        # "thresh": "gau",
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)

    # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    # opt = RFOptimizer(geom)

    opt.run()


# def test_abnr_xtb():
    # from pysisyphus.optimizers.ABNR import ABNR
    # geom = geom_from_library("split.image_021.xyz")
    # # geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    # xtb = XTB()
    # geom.set_calculator(xtb)

    # opt_kwargs = {
        # "max_cycles": 15,
        # # "hist_max": 10,
        # # "dump": True,
        # # "trust_radius": 0.1,
    # }
    # opt = ABNR(geom, **opt_kwargs)

    # # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    # # opt = RFOptimizer(geom)

    # opt.run()


if __name__ == "__main__":
    # test_sqnm()
    # test_sqnm_bio_mode()
    # test_sqnm_xtb()
    test_mecoome_sqnm_xtb()
    # test_mecoome_sqnm_g16()
    # test_abnr_xtb()

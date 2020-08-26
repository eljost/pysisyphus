import numpy as np

from pysisyphus.calculators import XTB
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.closures import bfgs_multiply

from BFGS import BFGS


def test_bfgs_multiply():
    size = 10
    s_list = list()
    y_list = list()
    forces = np.ones(size)
    Hf = bfgs_multiply(s_list, y_list, forces)
    np.testing.assert_allclose(Hf, -forces)


def test_double_damped_neb():
    geoms = geom_loader("interpolated.trj")
    # geoms = geom_loader("cycle_034.input.trj")
    for i, geom in enumerate(geoms):
        calc_kwargs = {
            "pal": 2,
            "calc_number": i,
            "charge": 0,
            "mult": 1,
        }
        calc = XTB(**calc_kwargs)
        geom.set_calculator(calc)
    cos = NEB(geoms)

    # quickmin_kwargs = {
        # "max_cycles": 35,
        # "dump": True,
    # }
    # opt = QuickMin(cos, **opt_kwargs)

    opt_kwargs = {
        "max_cycles": 10,
        "max_step": 0.1,
        # "update": "damped",
        # "update": "bfgs",
        "update": "double",
        # "align": True,  # does not work!
        "dump": True,
    }
    opt = BFGS(cos, **opt_kwargs)

    # opt_kwargs = {
        # # "align": True,
        # "max_step": 0.1,
    # }
    # opt = LBFGS(cos, **opt_kwargs)

    opt.run()

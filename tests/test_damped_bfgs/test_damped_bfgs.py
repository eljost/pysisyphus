import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.calculators import XTB
from pysisyphus.testing import using


def test_bfgs_multiply_empty_lists():
    size = 10
    s_list = list()
    y_list = list()
    forces = np.ones(size)
    Hf = bfgs_multiply(s_list, y_list, forces)
    np.testing.assert_allclose(Hf, forces)


@pytest.mark.skip_ci
@using("xtb")
@pytest.mark.parametrize(
    "opt_cls, _opt_kwargs, ref_cycle", [
        # The first two should yield identical results
        (BFGS, {"update": "bfgs", }, 25),
        # Deactivate gamma_mult to make it comparable to BFGS
        (LBFGS, {"keep_last": 25, "double_damp": False, "gamma_mult": False, }, 25),
        # s is nearly never damped, so this is not much different to undamped BFGS
        (BFGS, {"update": "damped", }, 25),
        # y gets damped often, so double damping makes a difference
        (BFGS, {"update": "double", }, 16),
        (LBFGS, {"keep_last": 10, "double_damp": True, "gamma_mult": True,
                 "align": False, }, 15),
        # Geometries are already fairly aligned, so this doesn't do much
        # (LBFGS, {"keep_last": 10, "double_damp": True, "gamma_mult": True,
                 # "align": True, }, 15),
    ]
)
def test_double_damped_neb(opt_cls, _opt_kwargs, ref_cycle, this_dir):
    geoms = geom_loader(this_dir / "neb_input.trj")
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

    opt_kwargs = {
        "max_step": 0.1,
        "rms_force": 0.005,
        "rms_force_only": True,
    }
    opt_kwargs.update(_opt_kwargs)

    opt = opt_cls(cos, **opt_kwargs)
    opt.run()

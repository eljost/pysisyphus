import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.DimerMethod import DimerMethod
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


def test_dimer():
    coords = (-0.2, 1.1, 0)
    geom = Geometry(("X", ), coords)
    N_init = np.array((0.3, 0.7, 0.))

    # New implementation
    dimer_kwargs = {
        "calculator": AnaPot(),
        "geometry": geom,
        "N_init": N_init,
    }
    dimer = DimerMethod(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt = PreconLBFGS(geom, precon=False, line_search=None, max_step_element=0.25,
                      thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 9
    assert geom.energy == pytest.approx(2.80910484)

    # AnaPot().plot_opt(opt)

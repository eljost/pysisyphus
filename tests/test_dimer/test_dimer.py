import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


@pytest.mark.parametrize(
    "rotation_method, ref_cycle",
    [
        ("direct", 11),
        ("fourier", 9),
    ]
)
def test_dimer(rotation_method, ref_cycle):
    coords = (-0.2, 1.1, 0)
    geom = Geometry(("X", ), coords)
    N_init = np.array((0.3, 0.7, 0.))

    # New implementation
    dimer_kwargs = {
        "rotation_method": rotation_method,
        "calculator": AnaPot(),
        "N_init": N_init,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "precon": False,
        "line_search": None,
        "max_step_element": 0.25,
        "thresh": "gau_tight",
        "max_cycles": 15,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(2.80910484)

    # AnaPot().plot_opt(opt)

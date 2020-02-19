import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library, do_final_hessian
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


init_logging()


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


def test_dimer_hcn():
    geom = geom_from_library("baker_ts/01_hcn.xyz")
    ref_energy = -92.24604
    N_init = " 0.5858  0.      0.0543 " \
             "-0.7697 -0.      0.061 " \
             "0.2027  0.     -0.1295".split()

    calc = PySCF("321g", pal=2)

    dimer_kwargs = {
        "rotation_method": "fourier",
        "calculator": calc,
        "N_init": N_init,
    }
    dimer = Dimer(**dimer_kwargs)
    geom.set_calculator(dimer)

    opt_kwargs = {
        "precon": True,
        "max_step_element": 0.25,
        "max_cycles": 15,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 8
    assert geom.energy == pytest.approx(ref_energy)

    # do_final_hessian(geom, save_hessian=False)

import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.cos.NEB import NEB
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.plot import plot_cos_energies, plot_cos_forces, plot_opt


def test_multi_run():
    # NEB
    geoms = AnaPot().get_path(10)
    cos = NEB(geoms)
    neb_kwargs = {
        "dump": True,
        "h5_group_name": "neb",
    }
    neb_opt = SteepestDescent(cos, **neb_kwargs)
    neb_opt.run()

    # GrowingString
    geoms = AnaPot().get_path(10)
    gsm = GrowingString(geoms, calc_getter=AnaPot, perp_thresh=0.25)
    gsm_kwargs = {
        "dump": True,
        "h5_group_name": "gsm",
        "stop_in_when_full": 0,
    }
    gsm_opt = StringOptimizer(gsm, **gsm_kwargs)
    gsm_opt.run()
    # calc = geoms[0].calculator
    # calc.anim_opt(opt, show=True)

    # Simple Optimization
    geom = AnaPot().get_path(10)[5]
    opt_kwargs = {
        "dump": True,
        "h5_group_name": "opt",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    # h5_fn = "optimization.h5"
    # h5_group = "neb"
    # plot_cos_energies(h5_fn, h5_group)

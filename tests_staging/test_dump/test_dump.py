import h5py
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.cos.NEB import NEB
from pysisyphus.interpolate import interpolate
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.run import run_from_dict


def test_dump_neb():
    # NEB
    geoms = AnaPot().get_path(10)
    cos = NEB(geoms)
    opt_kwargs = {
        "dump": True,
        # "rms_force": 1.0,
        "h5_group_name": "neb",
    }
    opt = SteepestDescent(cos, **opt_kwargs)
    # opt = SteepestDescent(geoms[7], **opt_kwargs)
    opt.run()

    # GrowingString
    geoms = AnaPot().get_path(10)
    gsm = GrowingString(geoms, calc_getter=AnaPot, perp_thresh=0.25)
    opt_kwargs = {
        "dump": True,
        "h5_group_name": "gsm",
        "gamma": 10.,
        "stop_in_when_full": 0,
    }
    opt = StringOptimizer(gsm, **opt_kwargs)
    opt.run()
    calc = geoms[0].calculator
    # calc.anim_opt(opt, show=True)

    # # Simple Optimization
    # geom = AnaPot().get_path(10)[5]
    # opt_kwargs = {
        # "dump": True,
        # "h5_group_name": "opt",
    # }
    # opt = RFOptimizer(geom, **opt_kwargs)
    # opt.run()
    return

    hei_coords, *_ = cos.get_splined_hei()
    ts_geom = geoms[0].copy()
    ts_geom.coords = hei_coords
    ts_geom.set_calculator(AnaPot())
    ts_opt_kwargs = {
        # "dump": True,
        "thresh": "gau_tight",
    }
    ts_opt = RSIRFOptimizer(ts_geom, **ts_opt_kwargs)
    ts_opt.run()

    # calc = geoms[0].calculator
    # calc.anim_opt(opt, show=True)


@pytest.mark.skip
def test_run():
    run_dict = {
        "preopt": {
            "max_cycles": 5,
            "coord_type": "cart",
        },
        "interpol": {
            "type": "linear",
            "between": 12,
            "align": False,
        },
        "cos": {
            "type": "neb",
        },
        "opt": {
            "type": "sd",
            "align": False,
            "rms_force": 1.0,
            "dump": False,
        },
        "tsopt": {
            "type": "rsirfo",
            "thresh": "gau_tight",
        },
        "calc": {
            "type": "anapot",
        },
        # Inline xyz
        "xyz": """
         1
         
         X -1.2 1.4 0.0
         1
         
         X 2.0 4.0 0.0
         """,
         "coord_type": "cart",
    }
    results = run_from_dict(run_dict)


def test_h5():
    fn = "test.h5"
    size = 5
    max_ = 5
    shape = (max_, size)
    maxshape = (None, size)
    with h5py.File(fn, "w") as f:
        opt = f.create_group("opt")
        test = opt.create_dataset("test", shape)
        # test.resize((10, size))
        test_ = f["opt/test"][:]

        unlim = f.create_dataset("unlim", shape, maxshape=maxshape)
        unlim.resize((10, size))
        unlim_ = f["unlim"][:]
    print(test_)
    print(unlim_)


def test_opt_dump():
    geom = AnaPot().get_geom((1, 1, 0.))
    calc = geom.calculator
    opt_kwargs = {
        "dump": True,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    calc.plot_opt(opt)#, show=True)


def test_h5_sizes():
    fn = "hess.h5"
    cycles1 = 10
    cycles2 = 100
    coord_size = 300
    hess_shape = (coord_size, coord_size)
    # maxshape = (None, size)
    with h5py.File(fn, "w") as f:
        hess1 = f.create_dataset("hess1", (cycles1, *hess_shape), dtype=np.float64)
        hess2 = f.create_dataset("hess2", (cycles2, *hess_shape), dtype=np.float64, chunks=True)
        import pdb; pdb.set_trace()
        pass

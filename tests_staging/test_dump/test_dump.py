import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.interpolate import interpolate
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.run import run_from_dict


def test_dump():
    geoms = AnaPot().get_path(10)
    cos = NEB(geoms)
    opt_kwargs = {
        "dump": True,
        "rms_force": 1.0,
        # "max_cycles": 5,
    }
    opt = SteepestDescent(cos, **opt_kwargs)
    # opt = SteepestDescent(geoms[7], **opt_kwargs)
    opt.run()
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


import h5py
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

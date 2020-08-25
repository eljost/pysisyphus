from pysisyphus.calculators import XTB
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_loader
from BFGS import BFGS
from pysisyphus.optimizers.QuickMin import QuickMin


def test_neb():
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
        "max_cycles": 35,
        "max_step": 0.1,
        # "update": "damped",
        # "update": "bfgs",
        "update": "double",
        "dump": True,
    }
    opt = BFGS(cos, **opt_kwargs)
    opt.run()

#!/usr/bin/env python3

from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.helpers import geom_from_library, geom_from_xyz_file, do_final_hessian
# from pysisyphus.optimizers.ANCOptimizer import ANCOptimizer
from pysisyphus.optimizers.NCOptimizer import NCOptimizer


# geom = geom_from_library("azetidine_guess.xyz")
# calc = Gaussian16("HF 4-31G", pal=4)
# geom = geom_from_xyz_file("guess.xyz")
#geom = geom_from_xyz_file("guess2.xyz")
geom = geom_from_xyz_file("guess3.xyz")
calc = Gaussian16("PM6", pal=4)
# from pysisyphus.calculators.XTB import XTB
# calc = XTB(pal=4)
geom.set_calculator(calc)

# opt = ANCOptimizer(geom, dump=True)
opt_kwargs = {
    "dump": True,
    "hessian_init": "calc",
    "freeze_modes": 200,
    "max_cycles": 20,
    "prefix": "frozen_"
}
opt = NCOptimizer(geom, **opt_kwargs)
opt.run()
do_final_hessian(geom)

# from pysisyphus.Geometry import Geometry
# from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
# geom = Geometry(geom.atoms, geom.coords, coord_type="redund", define_prims=((20, 19),))
# geom.set_calculator(calc)
# tsopt = RSIRFOptimizer(geom, hessian_recalc=5, trust_max=0.3)
# tsopt.run()
# do_final_hessian(geom)

# from pysisyphus.irc.EulerPC import EulerPC
# geom = Geometry(geom.atoms, geom.cart_coords)
# geom.set_calculator(calc)
# irc = EulerPC(geom)
# irc.run()

#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.tsoptimizers.PRFOptimizer import PRFOptimizer
from pysisyphus.optimizers.RSRFOptimizer import RSRFOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.calculators.XTB import XTB


geom = geom_from_library("hcn_iso_ts.xyz")
# geom = geom_from_library("codein.xyz")
xtb = XTB()
geom.set_calculator(xtb)

# opt_kwargs = {
    # "max_cycles": 5,
    # # "thresh": "gau",
    # "max_size": 0.4,
    # # "recalc_hess": 2,
    # "max_micro_cycles": 15,
    # # "trust_radius": 0.06,
    # "trust_radius": 0.3,
# }
# opt = RSPRFOptimizer(geom, **opt_kwargs)
# opt.run()

min_opt_kwargs = {
    # "max_cycles": 3,
    "recalc_hess": 10,
    # "hess_update": "flowchart",
    "hess_update": "flowchart",
    # "thresh": "gau_tight",
    # "max_micro_cycles": 1,
    "trust_radius": 0.1,
}
opt = RSRFOptimizer(geom, **min_opt_kwargs)
# opt = RFOptimizer(geom, **min_opt_kwargs)
opt.run()

with open("opt.xyz", "w") as handle:
    handle.write(geom.as_xyz())

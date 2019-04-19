#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.tsoptimizers.PRFOptimizer import PRFOptimizer
from pysisyphus.calculators.XTB import XTB

geom = geom_from_library("hcn_iso_ts.xyz")
xtb = XTB()
geom.set_calculator(xtb)

opt_kwargs = {
    # "max_cycles": 8,
    # "thresh": "gau",
    "max_size": 0.4,
    # "recalc_hess": 2,
}
opt = PRFOptimizer(geom, **opt_kwargs)
opt.run()

with open("opt.xyz", "w") as handle:
    handle.write(geom.as_xyz())

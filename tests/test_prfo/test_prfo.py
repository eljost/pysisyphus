#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.tsoptimizers.PRFOptimizer import PRFOptimizer
from pysisyphus.optimizers.RSRFOptimizer import RSRFOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.calculators.XTB import XTB


geom = geom_from_library("hcn_iso_ts.xyz")
# geom = geom_from_library("codein.xyz")
xtb = XTB(pal=1)
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
    "recalc_hess": 1,
    # "hess_update": "flowchart",
    "hess_update": "flowchart",
    "thresh": "gau",
    # "max_micro_cycles": 1,
    "trust_radius": 0.1,
}
opt = RSRFOptimizer(geom, **min_opt_kwargs)
# opt = RFOptimizer(geom, **min_opt_kwargs)
opt.run()

with open("opt.xyz", "w") as handle:
    handle.write(geom.as_xyz())

H = geom.hessian
f = geom.forces
import numpy as np
import pdb; pdb.set_trace()
step = np.linalg.inv(H).dot(f)
step_rms = np.sqrt(np.mean(step**2))
step_max = np.abs(step).max()
step_norm = np.linalg.norm(step)
print(f"rms(step)={step_rms:.6f}, max(step)={step_max:.6f}, norm(step)={step_norm:.6f}")

#!/usr/bin/env python3

from pysisyphus.helpers import get_baker_opt_ts_geoms
from pysisyphus.color import green, red
from pysisyphus.calculators import Gaussian16
from pysisyphus.irc import EulerPC


colors = {
    True: green,
    False: red,
}

geoms, meta_data = get_baker_opt_ts_geoms()
only_fn = "01_hcn_opt_ts.xyz"
meta_data = {
    only_fn: meta_data[only_fn]
}

for i, (xyz_fn, (charge, mult)) in enumerate(meta_data.items()):
    geom = geoms[xyz_fn]
    print(f"{i:02d}: {xyz_fn} ({charge},{mult})")

    calc_kwargs = {
        "route": "HF/3-21G",
        "pal": 4,
        "charge": charge,
        "mult": mult,
    }
    geom.set_calculator(Gaussian16(**calc_kwargs))

    irc_kwargs = {
        "step_length": .2,
        "hessian_recalc": 5,
    }
    irc = EulerPC(geom, **irc_kwargs)
    irc.run()
    col = colors[irc.converged]
    print()
    print(col(f"Converged: {irc.converged}"))

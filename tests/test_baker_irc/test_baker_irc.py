#!/usr/bin/env python3

import itertools as it
from pathlib import Path
import shutil

from pysisyphus.helpers import get_baker_opt_ts_geoms
from pysisyphus.color import green, red
from pysisyphus.calculators import Gaussian16
from pysisyphus.irc import EulerPC


colors = {
    True: green,
    False: red,
}

geoms, meta_data = get_baker_opt_ts_geoms()
# only_fn = "01_hcn_opt_ts.xyz"
# only_fn = "16_h2po4_anion_opt_ts.xyz"
# meta_data = {
    # only_fn: meta_data[only_fn]
# }

works = (
    "01_hcn_opt_ts.xyz",
    "03_h2co_opt_ts.xyz",
    "04_ch3o_opt_ts.xyz",
    "05_cyclopropyl_opt_ts.xyz",
    "06_bicyclobutane_opt_ts.xyz",
    "07_bicyclobutane_opt_ts.xyz",
    "08_formyloxyethyl_opt_ts.xyz",
    "13_hf_abstraction_opt_ts.xyz",
    "14_vinyl_alcohol_opt_ts.xyz",
    "16_h2po4_anion_opt_ts.xyz",
    "18_silyene_insertion_opt_ts.xyz",
    "19_hnccs_opt_ts.xyz",
    "21_acrolein_rot_opt_ts.xyz",
    "20_hconh3_cation_opt_ts.xyz",
    "23_hcn_h2_opt_ts.xyz",
    "24_h2cnh_opt_ts.xyz",
    "25_hcnh2_opt_ts.xyz",
)

fails = (
    # "05_cyclopropyl_opt_ts.xyz",
)

use = (
    works,
    fails,
)
use_meta_data = dict()
for key in it.chain(*use):
    use_meta_data[key] = meta_data[key]

for i, (xyz_fn, (charge, mult)) in enumerate(use_meta_data.items()):
    geom = geoms[xyz_fn]
    print(f"@{i:02d}: {xyz_fn} ({charge},{mult})")

    calc_kwargs = {
        "route": "HF/3-21G",
        "pal": 4,
        "charge": charge,
        "mult": mult,
    }
    # Otherwise SCF converges to wrong solution with no imaginary frequencies
    if xyz_fn == "05_cyclopropyl_opt_ts.xyz":
        calc_kwargs["route"] += " scf=qc"
    geom.set_calculator(Gaussian16(**calc_kwargs))

    irc_kwargs = {
        # "hessian_recalc": 5,
        # "displ_energy": 1e-3,
        # "step_length": .4,
        "hessian_recalc": 10,
        # "displ": "length",
        # "backward": False,
        # "forward": False,
    }
    irc = EulerPC(geom, **irc_kwargs)
    irc.run()
    col = colors[irc.converged]
    print()
    print(col(f"@Converged: {irc.converged}"))
    src = "finished_irc.trj"
    dst = Path(xyz_fn).stem + "_finished_irc.trj"
    shutil.copy(src, dst)
    print(f"@Copied '{src}' to '{dst}'")

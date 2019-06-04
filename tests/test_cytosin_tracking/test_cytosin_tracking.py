#!/usr/bin/env python3

import itertools as it
import os
from pathlib import Path

from pysisyphus.helpers import geom_from_library, highlight_text
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


def get_calc(ovlp_type, ovlp_with):
    calc_kwargs = {
        "control_path": THIS_DIR / "control_pbe0_clean",
        "track": True,
        "ovlp_type": ovlp_type,
        "ovlp_with": ovlp_with,
        # "cdds": "render",
        # "orient": "reset;center {-3.5575373 2.2096293 1.4794437E-5}; " \
                    # "rotate z -112.54; rotate y 41.57; rotate z 99.02;",
        "adapt_args": [0.5, 0.3, 0.6],
        "pal": 4,
    }
    calc = Turbomole(**calc_kwargs)
    return calc


def run():
    ovlp_types = "wf tden nto_org nto".split()
    # ovlp_types = ("nto", )
    ovlp_withs = "adapt first previous".split()
    for i, (ovlp_type, ovlp_with) in enumerate(it.product(ovlp_types, ovlp_withs)):
        # ovlp_type = "wf"
        # ovlp_with = "adapt"
        print(
            highlight_text(f"i={i:02d}, ovlp_type={ovlp_type}, ovlp_with={ovlp_with}")
        )
        geom = geom_from_library("cytosin.xyz", coord_type="redund")
        calc = get_calc(ovlp_type, ovlp_with)
        geom.set_calculator(calc)
        opt = RFOptimizer(geom)
        opt.run()
        assert calc.root_flips[2]  # == True
        assert all([flipped == False for i, flipped in enumerate(calc.root_flips)
                    if i != 2]
        )
        assert calc.root == 2
        assert opt.cur_cycle == 4
        print()


if __name__ == "__main__":
    run()

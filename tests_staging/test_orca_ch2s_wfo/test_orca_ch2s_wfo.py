#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.ORCA import ORCA


np.set_printoptions(precision=4, suppress=True)
GEOM = geom_from_library("ch2s_bp86def2sv_opt.xyz")


def self_compare_base(geom, wfo_basis):
    calc_kwargs = {
        "keywords": "BP86 def2-SVP",
        "track": False,
        "wfo_basis": wfo_basis,
        "wfo_charge": 0,
    }

    orca = ORCA(**calc_kwargs)
    orca.run_calculation(geom.atoms, geom.coords)
    #wfow = turbo.wfow
    #wfow.compare(wfow)
    #turbo.root = 0
    # turbo.root = 0
    # turbo.track_root(geom.atoms, geom.coords)
    #wfow.track()


def test_bp86():
    self_compare_base(GEOM, "def2svp")


if __name__ == "__main__":
    test_bp86()

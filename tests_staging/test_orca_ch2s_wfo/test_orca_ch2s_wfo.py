#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.calculators.ORCA import ORCA


np.set_printoptions(precision=4, suppress=True)
init_logging()
GEOM = geom_from_library("ch2s_bp86def2sv_opt.xyz")
#GEOM = geom_from_library("h2.xyz")


def self_compare_base(geom, wfo_basis):

    calc_kwargs = {
        "keywords": "B3LYP def2-TZVP",
        "blocks": "%tddft nroots 2 maxdim 5 end",
        "track": True,
        "pal": 4,
    }

    orca = ORCA(**calc_kwargs)
    orca.run_calculation(geom.atoms, geom.coords)
    # cis = "calculator_0.000.orca.cis"
    # gbw = "calculator_0.000.orca.gbw"
    # out = "calculator_0.000.orca.out"
    # orca.cis = cis
    # orca.gbw = gbw
    # orca.out = out

    orca.store_wfo_data(geom.atoms, geom.coords)
    wfow = orca.wfow
    res = wfow.compare(wfow)
    print(res)

    # mo_coeffs = orca.parse_gbw(gbw)
    # print(mo_coeffs)
    # print(mo_coeffs.shape)
    #mos = fake_turbo_mos(coeffs)
    #orca.parse_mo_numbers(out)

    # coeffs = orca.parse_cis(cis)
    # print(coeffs)
    # print(coeffs[0])
    # print(coeffs.shape)

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

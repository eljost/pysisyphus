#!/usr/bin/env python3

from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.findiffs import verify_geom


def test_wilson_by_findiff():
    geom = geom_loader("lib:h2o2_hf_321g_opt.xyz")
    # geom = geom_loader("lib:baker_ts/19_hnccs.xyz")
    assert verify_geom(geom)

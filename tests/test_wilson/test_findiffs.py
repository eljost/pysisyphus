#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.intcoords.findiffs import verify_geom


def test_wilson_by_findiff():
    geom = geom_from_library("h2o2_hf_321g_opt.xyz")
    assert verify_geom(geom)

#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.intcoords.findiffs import test_geom


def run():
    geom = geom_from_library("h2o2_hf_321g_opt.xyz")
    test_geom(geom)


if __name__ == "__main__":
    run()

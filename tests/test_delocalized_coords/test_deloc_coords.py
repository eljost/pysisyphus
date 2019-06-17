#!/usr/bin/env python3


from pysisyphus.calculators.XTB import XTB
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_from_library


def run():
    # geom = geom_from_library("split.image_021.xyz")
    # geom.set_calculator(XTB())
    # opt = RFOptimizer(geom)
    # opt.run()

    # geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    # geom.set_calculator(XTB())
    # opt = RFOptimizer(geom)
    # opt.run()

    xyz_fn = "fluorethylene.xyz"
    # xyz_fn = "split.image_021.xyz"
    geom = geom_from_library(xyz_fn, coord_type="dlc")
    # int_ = geom.internal
    # int_.B
    # int_.set_delocalized_vectors()
    geom.set_calculator(XTB())
    opt_kwargs = {
        # "max_cycles": 1,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()


if __name__ == "__main__":
    run()

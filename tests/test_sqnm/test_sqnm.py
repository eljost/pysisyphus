#!/usr/bin/env python3

import matplotlib.pyplot as plt

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.StabilizedQNMethod import StabilizedQNMethod


def test_sqnm():
    # geom = AnaPot.get_geom((0, 4, 0))
    geom = AnaPot.get_geom((-0.8, 1.73, 0))
    # geom = AnaPot.get_geom((-1, 1.05, 0))
    # calc = geom.calculator 
    # calc.plot()
    # plt.show()

    opt_kwargs = {
        "max_cycles": 15,
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)
    opt.run()


if __name__ == "__main__":
    test_sqnm()

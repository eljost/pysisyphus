import numpy as np

from pysisyphus.calculators.AnaPotBase import AnaPotBase


class AnaPot(AnaPotBase):
    def __init__(self, **kwargs):
        kwargs_ = {
            "V_str": "4 + 4.5*x - 4*y + x**2 + 2*y**2-2*x*y + x**4 - 2*x**2*y",
            "xlim": (-2, 2.5),
            "ylim": (0, 5),
            "levels": np.linspace(-3, 4, 80),
            "minima": ((-1.05274, 1.02776, 0), (1.94101, 3.85427, 0)),
            "saddles": ((0.6117313, 1.4929732, 0.0),),
        }
        kwargs_.update(kwargs)
        super().__init__(**kwargs_)

    def __str__(self):
        return "AnaPot calculator"

from pysisyphus.calculators.AnaPotBase import AnaPotBase


class Quadratic1dPot(AnaPotBase):
    def __init__(self, **kwargs):
        kwargs_ = {
            "V_str": "-x**2",
            "saddles": ((0.0, 0.0, 0.0),),
            "xlim": (-1, 1),
        }
        kwargs_.update(kwargs)
        super().__init__(**kwargs_)

    def __str__(self):
        return "Quadratic1dPot calculator"


class Cubic1dPot(AnaPotBase):
    def __init__(self, **kwargs):
        kwargs_ = {
            "V_str": "0.7 * x**3 - 2.0 * x**2 + 1.0 * x",
            "saddles": ((0.2959976785, 0.0, 0.0),),
            "xlim": (-2, 2),
        }
        kwargs_.update(kwargs)
        super().__init__(**kwargs_)

    def __str__(self):
        return "Cubic1dPot calculator"

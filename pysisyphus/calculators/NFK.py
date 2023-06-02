# [1] https://doi.org/10.1016/j.cplett.2004.07.079
#     The reaction pathway of a potential energy surface as curve
#     with induced tangent
#     Hirsch, Quapp, 2004


import numpy as np

from pysisyphus.calculators.AnaPotBase import AnaPotBase


class NFK(AnaPotBase):
    def __init__(self, C=0.06, **kwargs):
        """NFK surface of Hirsch and Quapp.

        See Fig. 4 in [1] and the appendix.
        """

        # Note that C is substituted in
        V_str = f"{C} * (x**2 + y**2)**2 + x*y - 9*exp(-(x-3)**2-y**2) - 9*exp(-(x+3)**2 - y**2)"
        kwargs_ = {
            "V_str": V_str,
            "xlim": (-3, 3),
            "ylim": (-3, 3),
            "levels": np.linspace(-7, 3.5, 70),
            "minima": ((2.71268103, -0.15093971, 0.0), (-2.7126810, 0.1509397, 0.0)),
            "saddles": ((0.0, 0.0, 0.0),),
        }
        kwargs_.update(kwargs)
        super().__init__(**kwargs_)

    def __str__(self):
        return "NFK calculator"


class ModNFK(NFK):
    def __init__(self, C=0.03, **kwargs):
        kwargs["minima"] = (
            (-2.84721930, 0.15702574, 0.0),
            (2.84721930, -0.15702574, 0.0),
        )
        super().__init__(C, **kwargs)

    def __str__(self):
        return "ModNFK calculator"

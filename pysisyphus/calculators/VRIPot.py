# [1] https://doi.org/10.1016/j.cplett.2021.138970
#     The dynamical significance of valley-ridge inflection points
#     Víctor J.García-Garrido, Stephen Wiggins

import numpy as np

from pysisyphus.calculators.AnaPotBase import AnaPotBase


class VRIPot(AnaPotBase):
    def __init__(self, **kwargs):
        VTS = 0.5
        xs = 1.0
        xi = 0.3265
        A = 1.0
        B = 1.0
        C = 1.0
        lim = 1.5
        lims = (-lim, lim)
        kwargs_ = {
            "V_str": (
                f"{VTS} / {xs}**4 * x**2 * (x**2 - 2*{xs}**2)"
                f"+ {A}*y**2 * ({xi}-x) "
                f"+ y**4 * ({B} + {C}*x)"
            ),
            "xlim": lims,
            "ylim": lims,
            "levels": np.linspace(-3 * VTS, 2 * VTS, 65),
            # "minima": ((-1.05274, 1.02776, 0), (1.94101, 3.85427, 0)),
            # "saddles": ((0.6117313, 1.4929732, 0.0),),
        }
        _ = kwargs_["V_str"]
        kwargs_.update(kwargs)
        super().__init__(**kwargs_)

    def __str__(self):
        return "VRIPot calculator"

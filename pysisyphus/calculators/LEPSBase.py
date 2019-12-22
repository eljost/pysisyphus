#!/usr/bin/env python3

from pysisyphus.calculators.AnaPotBase import AnaPotBase
from pysisyphus.calculators.LEPSExpr import LEPSExpr


class LEPSBase(AnaPotBase):

    def __init__(self, pot_type="leps"):
        leps_expr = LEPSExpr()
        V_expr, xlim, ylim, levels = leps_expr.get_expr(pot_type)
        V_str = str(V_expr)

        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim, levels=levels)

    def __str__(self):
        return "LEPSBase calculator"


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    choices = "leps harmonic tot dimer".split()
    for c in choices:
        lp = LEPSBase(pot_type=c)
        lp.plot()
        plt.show()

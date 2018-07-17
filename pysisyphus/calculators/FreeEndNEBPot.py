import numpy as np
from sympy import atan, symbols

from pysisyphus.calculators.AnaPotBase import AnaPotBase
from pysisyphus.calculators.LEPSExpr import LEPSExpr


class FreeEndNEBPot(AnaPotBase):

    def __init__(self):
        leps_expr = LEPSExpr()
        V_expr, xlim, ylim, levels = leps_expr.get_expr("harmonic")

        self.x0 = 1.93
        x = symbols("x")
        
        V_expr = V_expr - 2*atan(5*(x-self.x0)) - 2*x
        V_str = str(V_expr)
        xlim = (0.5, 2.25)
        ylim = (0.0, 1.6)
        levels = np.linspace(-10, 0, 125)

        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim, levels=levels)

    def __str__(self):
        return "FreeEndNEBPot calculator"


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    fep = FreeEndNEBPot()
    fep.plot()
    plt.show()

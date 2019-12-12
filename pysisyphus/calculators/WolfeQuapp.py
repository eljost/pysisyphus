from pysisyphus.calculators.AnaPotBase import AnaPotBase


class WolfeQuapp(AnaPotBase):

    def __init__(self):
        V_str = "x**4 + y**4 - 2*x**2 - 4*y**2 + x*y + 0.3*x + 0.1*y"
        xlim = (-2.0, 2.0)
        ylim = (-2.0, 2.0)
        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim)

    def __str__(self):
        return "WolfeQuapp calculator"

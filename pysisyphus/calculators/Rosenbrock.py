import numpy as np

from pysisyphus.calculators.AnaPotBase import AnaPotBase

class Rosenbrock(AnaPotBase):

    def __init__(self): 
        V_str = "(1-x)**2 + 100*(y - x**2)**2"
        xlim = (-2.5, 2.5)
        ylim = (-1.5, 3.5)
        levels = np.logspace(-5, 10, 50, base=2)
        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim, levels=levels)

    def __str__(self):
        return "Rosenbrock calculator"

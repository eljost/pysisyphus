import numpy as np

from pysisyphus.calculators.AnaPotBase import AnaPotBase


class Rastrigin(AnaPotBase):
    """
        http://www.sfu.ca/~ssurjano/rastr.html
    """

    def __init__(self): 
        # d == 2
        # f(x) = 10*d + sum_(i+1)^d (x_i**2 - 10*cos(2*pi*x_i))
        V_str = "20 + x**2 - 10*cos(2*pi*x) + y**2 - 10*cos(2*pi*y)"
        xlim = (-5.12, 5.12)
        ylim = xlim
        levels = np.linspace(0, 50, 50)
        minima = ((0., 0., 0.), )
        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim, levels=levels, minima=minima)

    def __str__(self):
        return "Rastrigin calculator"

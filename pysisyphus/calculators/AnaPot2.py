import numpy as np
import sympy as sym
from sympy import atan, exp, tan, sin, pi


from pysisyphus.calculators.AnaPotBase import AnaPotBase
from pysisyphus.calculators.Calculator import Calculator

# https://www.wolframalpha.com/input/?i=plot+arccot(-exp(y)*cot(x/2-pi/4))+-+2*exp(-(y-sin(x))^2/2)
# [1] http://aip.scitation.org/doi/abs/10.1063/1.461606
# https://www.wolframalpha.com/input/?i=derivative+of+(arccot(-exp(y)*cot(x/2-pi/4))+-+2*exp(-(y-sin(x))^2/2))

class AnaPot2_(Calculator):

    def __init__(self): 
        super().__init__()

    def get_energy(self, atoms, coords):
        x, y, z = coords
        cot = 1 / np.tan(x/2 - np.pi/4)
        arccot = np.arctan(
                1 / (-np.exp(y) * cot)
        )
        energy = (
            arccot - 2 * np.exp(-(y-np.sin(x))**2 / 2)
        )
        return {"energy": energy}

    def __str__(self):
        return "AnaPot2 calculator"


class AnaPot2(AnaPotBase):
    """We can't use sympify as it replaces 1/tan by cot and this isn't
    supported by numpy when we call lambdify."""

    def __init__(self): 
        x, y = sym.symbols("x y")
        V_str = atan(exp(-y)/tan(x/2 + pi/4)) - 2*exp(-(y - sin(x))**2/2)
        # xlim = (-np.pi/2, np.pi)
        xlim = (-np.pi, np.pi)
        ylim = (-2, 2)
        levels = np.linspace(-2, 1, 40)
        super().__init__(V_str=V_str, xlim=xlim, ylim=ylim, levels=levels,
                         use_sympify=False)

    def __str__(self):
        return "AnaPot2 calculator"

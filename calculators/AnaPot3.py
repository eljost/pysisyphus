import numpy as np
from sympy import symbols, diff, lambdify

from pysisyphus.calculators.Calculator import Calculator

# [1] http://www.cims.nyu.edu/~eve2/string_jcp_simplified07.pdf

class AnaPot3(Calculator):

    def __init__(self): 
        super(AnaPot3, self).__init__()
        x, y = symbols("x y")
        V = (1 - x**2 - y**2)**2 + (y**2) / (x**2 + y**2)
        dVdx = diff(V, x)
        dVdy = diff(V, y)
        self.V = lambdify((x, y), V, "numpy")
        self.dVdx = lambdify((x, y), dVdx, "numpy")
        self.dVdy = lambdify((x, y), dVdy, "numpy")

    def get_energy(self, atoms, coords):
        x, y, z = coords
        energy = self.V(x, y)
        return {"energy": energy}

    def get_forces(self, atoms, coords):
        x, y, z = coords
        dVdx = self.dVdx(x, y)
        dVdy = self.dVdy(x, y)
        dVdz = np.zeros_like(dVdx)
        forces = -np.array((dVdx, dVdy, dVdz))
        results = self.get_energy(atoms, coords)
        results["forces"] = forces
        return results

    def __str__(self):
        return "AnaPot3 calculator"

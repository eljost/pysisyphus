import numpy as np
from sympy import symbols, diff, lambdify, sympify

from pysisyphus.calculators.Calculator import Calculator

class AnaPotBase(Calculator):

    def __init__(self, V_str): 
        super(AnaPotBase, self).__init__()
        x, y = symbols("x y")
        V = sympify(V_str)
        dVdx = diff(V, x)
        dVdy = diff(V, y)
        self.V = lambdify((x, y), V, "numpy")
        self.dVdx = lambdify((x, y), dVdx, "numpy")
        self.dVdy = lambdify((x, y), dVdy, "numpy")

        dVdxdx = diff(V, x, x)
        dVdxdy = diff(V, x, y)
        dVdydy = diff(V, y, y)

        self.dVdxdx = lambdify((x, y), dVdxdx, "numpy")
        self.dVdxdy = lambdify((x, y), dVdxdy, "numpy")
        self.dVdydy = lambdify((x, y), dVdydy, "numpy")

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

    def get_hessian(self, atoms, coords):
        x, y, z = coords
        dVdxdx = self.dVdxdx(x, y)
        dVdxdy = self.dVdxdy(x, y)
        dVdydy = self.dVdydy(x, y)
        hessian = np.array(((dVdxdx, dVdxdy, 0),
                            (dVdxdy, dVdydy, 0),
                            (0, 0, 0))
        )
        results = self.get_forces(atoms, coords)
        results["hessian"] = hessian
        return results

#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, diff, lambdify, sympify

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.Geometry import Geometry

class AnaPotBase(Calculator):

    def __init__(self, V_str, xlim=(-1,1), ylim=(-1,1), levels=None,
                 use_sympify=True):
        super(AnaPotBase, self).__init__()
        self.xlim = xlim
        self.ylim = ylim
        self.levels = levels
        x, y = symbols("x y")
        if use_sympify:
            V = sympify(V_str)
        else:
            V = V_str
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

        self.fake_atoms = ("X", ) # X, dummy atom

        self.analytical_2d = True

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

    def plot(self, levels=None, show=False, **figkwargs):
        self.fig, self.ax = plt.subplots(**figkwargs)
        x = np.linspace(*self.xlim, 100)
        y = np.linspace(*self.ylim, 100)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, 0)
        pot_coords = np.stack((X, Y, Z))
        pot = self.get_energy(self.fake_atoms, pot_coords)["energy"]

        if levels is None:
            if self.levels is None:
                levels = np.linspace(pot.min(), pot.max(), 35)
            else:
                levels = self.levels

        # Draw the contourlines of the potential
        contours = self.ax.contour(X, Y, pot, levels)
        self.fig.colorbar(contours)

        if show:
            plt.show()

    def plot_eigenvalue_structure(self, grid=50, levels=None):
        self.plot(levels=levels)
        xs = np.linspace(*self.xlim, grid)
        ys = np.linspace(*self.ylim, grid)
        X, Y = np.meshgrid(xs, ys)
        z = list()
        for x_, y_ in zip(X.flatten(), Y.flatten()):
            H = self.get_hessian(self.fake_atoms,  (x_, y_, 0))["hessian"]
            w, v = np.linalg.eigh(H)
            z.append(
                1 if (w < 0).any() else 0
            )
        Z = np.array(z).reshape(X.shape)
        self.ax.contourf(X, Y, Z, cmap=cm.Reds)#, alpha=0.5)


    @classmethod
    def get_geom(cls, coords, atoms=("X", )):
        geom = Geometry(atoms, coords)
        geom.set_calculator(cls())
        return geom

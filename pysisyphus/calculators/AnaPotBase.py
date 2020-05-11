#!/usr/bin/env python3

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, diff, lambdify, sympify

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.Geometry import Geometry
from pysisyphus.interpolate import interpolate
from pysisyphus.plotters.AnimPlot import AnimPlot


class AnaPotBase(Calculator):

    def __init__(self, V_str, xlim=(-1,1), ylim=(-1,1), levels=None,
                 use_sympify=True, minima=None):
        super(AnaPotBase, self).__init__()
        self.xlim = xlim
        self.ylim = ylim
        self.levels = levels
        if minima is None:
            minima = list()
        self.minima = np.array(minima, dtype=float)

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

    def plot_eigenvalue_structure(self, grid=50, levels=None, show=False):
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
        if show:
            plt.show()

    def plot_opt(self, opt, enum=True):
        coords = np.array(opt.coords)
        self.plot()
        xs, ys = coords.T[:2]
        self.ax.plot(xs, ys, "o-")
        if enum:
            for i, (x, y) in enumerate(zip(xs, ys)):
                self.ax.annotate(i, (x, y))
        plt.show()

    def anim_opt(self, opt, energy_profile=False, colorbar=False, figsize=(8, 6),
                 show=False):
        try:
            min_ = self.levels.min()
            max_ = self.levels.max()
            num = self.levels.size
            levels = (min_, max_, num)
        except TypeError:
            levels = None

        anim = AnimPlot(self.__class__(),
                        opt,
                        xlim=self.xlim,
                        ylim=self.ylim,
                        levels=levels,
                        energy_profile=energy_profile,
                        colorbar=colorbar,
                        figsize=figsize
        )
        anim.animate()
        if show:
            plt.show()
        return anim

    @classmethod
    def get_geom(cls, coords, atoms=("X", ), V_str=None):
        geom = Geometry(atoms, coords)
        if V_str:
            geom.set_calculator(cls(V_str=V_str))
        else:
            geom.set_calculator(cls())
        return geom

    def get_path(self, num, minima_inds=None):
        between = num - 2

        inds = 0, 1
        if minima_inds is not None:
            inds = minima_inds
        initial_ind, final_ind = inds

        initial_geom = self.get_geom(self.minima[initial_ind])
        final_geom = self.get_geom(self.minima[final_ind])
        geoms = interpolate(initial_geom, final_geom, between=between)
        for geom in geoms:
            # Creating new instances can be really slow when the sympy calls
            # need some time. For now we just reuse the current calculator...
            # geom.set_calculator(self.__class__())
            geom.set_calculator(self)
        return geoms

    def get_minima(self):
        return [self.get_geom(coords) for coords in self.minima]

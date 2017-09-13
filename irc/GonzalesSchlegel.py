#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import newton

from pysisyphus.Geometry import Geometry
from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter

# [1] An improved algorithm for reaction path following
# http://aip.scitation.org/doi/pdf/10.1063/1.456010

class GonzalesSchlegel(IRC):

    def __init__(self, geometry, **kwargs):
        super(GonzalesSchlegel, self).__init__(geometry, **kwargs)

        self.pivot_coords = list()
        self.micro_coords = list()

        micro_header = "# |dx| |tangent|".split()
        micro_fmts = ["d", ".2E", ".3E"]
        self.micro_formatter = TableFormatter(micro_header, micro_fmts, 10)

    def bfgs_update(self, grad_diffs, coord_diffs):
        """BFGS update of the hessian as decribed in Eq. (16) of [1]."""
        y = grad_diffs[:,None]
        yT = grad_diffs[None,:]
        s = coord_diffs[:,None]
        sT = coord_diffs[None,:]
        first_term = y.dot(yT) / yT.dot(s)
        second_term = (self.hessian.dot(s).dot(sT).dot(self.hessian)
                       / sT.dot(self.hessian).dot(s)
        )
        dH = first_term - second_term
        return self.hessian + dH

    def micro_step(self):
        """Constrained optimization on a hypersphere."""
        eye = np.eye(self.displacement.size)

        gradient = -self.geometry.forces
        gradient_diff = gradient - self.last_gradient
        coords_diff = self.geometry.coords - self.last_coords
        self.last_gradient = gradient
        # Without copy we would only store the reference...
        self.last_coords = copy.copy(self.geometry.coords)

        self.hessian = self.bfgs_update(gradient_diff, coords_diff)
        eigvals, eigvecs = np.linalg.eig(self.hessian)
        hessian_inv = np.linalg.inv(self.hessian)

        def lambda_func(lambda_):
            # Eq. (11) in [1]
            # (H - λI)^-1
            hmlinv = np.linalg.inv(self.hessian - eye*lambda_)
            # (g - λp)
            glp = gradient - self.displacement*lambda_
            tmp = self.displacement - hmlinv.dot(glp)
            return tmp.dot(tmp) - 0.25*(self.step_length**2)

        smallest_eigval = np.sort(eigvals)[0]
        # Initial guess for λ.
        # λ must be smaller then the smallest eigenvector
        lambda_ = np.sort(eigvals)[0]
        lambda_ *= 1.5 if (lambda_ < 0) else 0.5
        # Find the root with scipy
        lambda_ = newton(lambda_func, lambda_)

        # Calculate dx from optimized lambda
        dx = -np.dot(
                np.linalg.inv(self.hessian-lambda_*eye),
                gradient-lambda_*self.displacement
        )
        self.displacement += dx
        self.geometry.coords += dx

        displ_norm = np.linalg.norm(self.displacement)
        tangent = gradient - gradient.dot(self.displacement)/displ_norm * gradient

        return dx, tangent
        

    def step(self):
        gradient0 = -self.geometry.forces
        gradient0_norm = np.linalg.norm(gradient0)
        self.irc_energies.append(self.geometry.energy)
        # For the BFGS update in the first micro step we use the original
        # point and the initial guess to calculate gradient and
        # coordinate differences.
        self.last_gradient = gradient0
        self.last_coords = self.geometry.coords

        # Take a step against the gradient to the pivot point x*_k+1.
        pivot_step = 0.5*self.step_length * gradient0/gradient0_norm
        pivot_coords = self.geometry.coords - pivot_step
        self.pivot_coords.append(pivot_coords)

        # Make initial guess for x'_k+1. Here we take another half
        # step from the pivot point.
        self.geometry.coords = pivot_coords - pivot_step
        # Initial displacement p' from the pivot point
        self.displacement = self.geometry.coords - self.pivot_coords[-1]

        these_micro_coords = list()
        i = 0
        print(self.micro_formatter.header)
        while True:
            if i == 20:
                logging.warning("Max micro cycles exceeded!")
                break
            #print(f"Micro {i}")
            dx, tangent = self.micro_step()
            these_micro_coords.append(self.geometry.coords)
            norm_dx = np.linalg.norm(dx)
            #print("norm(dx) = {:.2E}".format(norm_dx))
            norm_tangent = np.linalg.norm(tangent)
            #print("norm(tangent) = {:.2E}".format(norm_tangent))
            print(self.micro_formatter.line(i, norm_dx, norm_tangent))

            if (np.linalg.norm(dx) <= 1e-3):
                break
            i += 1

        self.micro_coords.append(np.array(these_micro_coords))

    def postprocess(self):
        self.pivot_coords = np.array(self.pivot_coords)
        self.micro_coords = np.array(self.micro_coords)

    def show2d(self):
        fig, ax = plt.subplots(figsize=(8,8))

        xlim = (-1.75, 1.25)
        ylim = (-0.5, 2.25)
        levels=(-150, -15, 40)
        x = np.linspace(*xlim, 100)
        y = np.linspace(*ylim, 100)
        X, Y = np.meshgrid(x, y)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y))
        pot = self.geometry.calculator.get_energy(fake_atoms, pot_coords)["energy"]
        levels = np.linspace(*levels)
        contours = ax.contour(X, Y, pot, levels)

        # Pivot points
        ax.plot(*zip(*self.pivot_coords), "bo", ls="-", label="pivot")
        # Constrained optmizations
        for mc in self.micro_coords:
            ax.plot(*zip(*mc), "yo", ls="-")
            for i, m in enumerate(mc):
                ax.text(*m, str(i))
        ax.plot(*zip(*self.coords), "ro", ls="-")
        plt.legend()
        plt.show()

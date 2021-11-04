# [1] https://pubs.acs.org/doi/pdf/10.1021/ct1000268
#     Updated Branching Plane for Finding Conical Intersections without Coupling
#     Derivative Vectors
#     Maeda, Ohno, Morokuma, 2010


from copy import deepcopy
from dataclasses import dataclass
from math import sqrt

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers_pure import hash_args


def update_y(x, x_prev, y_prev):
    """Update approximate coupling derivative vector y."""

    y_prev_x = y_prev.dot(x)
    x_prev_x = x_prev.dot(x)
    y = (y_prev_x * x_prev - x_prev_x * y_prev) / sqrt(y_prev_x ** 2 + x_prev_x ** 2)
    return y


def get_P(x, y):
    """Projector to project out components in branching plane."""

    I = np.eye(x.size)
    x_ = x / np.linalg.norm(x)
    y_ = y / np.linalg.norm(y)
    P = I - np.outer(x_, x_) - np.outer(y_, y_)
    return P


@dataclass
class CIQuantities:
    energy1: float
    gradient1: np.ndarray
    energy2: float
    gradient2: np.ndarray
    energy_diff: float
    gradient_diff: np.ndarray
    gradient_mean: np.ndarray
    P: np.ndarray
    # Difference gradient vector
    x: np.ndarray
    # Derivative coupling vector
    y: np.ndarray
    energy: float
    forces: np.ndarray


class ConicalIntersection(Calculator):
    """Calculator for conical intersection optimization.

    Based on [1].
    """

    def __init__(self, calculator1, calculator2, **kwargs):
        super().__init__(**kwargs)

        self.calculator1 = calculator1
        self.calculator2 = calculator2

        self.x_prev = None
        self.y_prev = None
        self.ciq_store = dict()

    def get_energy(self, atoms, coords, **prepare_kwargs):
        """Energy of calculator 1."""
        return self.calculator1.get_energy(atoms, coords, **prepare_kwargs)

    def get_ci_quantities(self, atoms, coords, **prepare_kwargs):
        """Relavent quantities including branching plane and projector P."""
        hash_ = hash_args(atoms, coords, *prepare_kwargs.values(), precision=6)
        if hash_ in self.ciq_store:
            ciq = self.ciq_store[hash_]
            self.log(f"Returning cached CI quantities for hash='{hash_}'.")
            return ciq
        else:
            self.log(f"Calculating CI quantities for hash='{hash_}'.")

        res1 = self.calculator1.get_forces(atoms, coords, **prepare_kwargs)
        energy1 = res1["energy"]
        gradient1 = -res1["forces"]

        res2 = self.calculator2.get_forces(atoms, coords, **prepare_kwargs)
        energy2 = res2["energy"]
        gradient2 = -res2["forces"]

        # Energy difference
        energy_diff = energy1 - energy2
        # Difference Gradient vector DGV
        gradient_diff = gradient1 - gradient2
        # Gradient mean
        gradient_mean = (gradient1 + gradient2) / 2

        # x vector; unit vector parallel to DGV
        x = gradient_diff / np.linalg.norm(gradient_diff)

        # Initialize y when it is not yet set ...
        if self.y_prev is None:
            y = gradient_mean / np.linalg.norm(gradient_mean)
        # ... or update it in later calculations
        else:
            y = update_y(x, self.x_prev, self.y_prev)

        # Orthogonalize against x
        y -= y.dot(x) * x
        # And renormalize
        y /= np.linalg.norm(y)

        # Projector to remove components along x and y (in BP) from gradient.
        # Eq. (7) in [1].
        P = get_P(x, y)

        self.log(f"ΔE={energy_diff:.6f} au ({energy_diff*AU2KJPERMOL:.2f} kJ mol⁻¹)")
        # Store vectors for next calculation
        self.x_prev = x
        self.y_prev = y

        gradient_diff_dash = 2 * energy_diff * x
        # Eq. (6) in [1].
        forces = -(gradient_diff_dash + P.dot(gradient_mean))

        ciq = CIQuantities(
            energy1=energy1,
            gradient1=gradient1,
            energy2=energy2,
            gradient2=gradient2,
            energy_diff=energy_diff,
            gradient_diff=gradient_diff,
            gradient_mean=gradient_mean,
            P=P,
            x=x,
            y=y,
            energy=energy1,
            forces=forces,
        )
        self.ciq_store[hash_] = deepcopy(ciq)
        return ciq

    def get_forces(self, atoms, coords, **prepare_kwargs):
        """Projected gradient for CI optimization."""
        ciq = self.get_ci_quantities(atoms, coords, **prepare_kwargs)

        return {
            "energy": ciq.energy,
            "forces": ciq.forces,
        }

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        """Projected Hessian."""
        ciq = self.get_ci_quantities(atoms, coords, **prepare_kwargs)

        hessian1 = self.calculator1.get_hessian(atoms, coords, **prepare_kwargs)[
            "hessian"
        ]
        hessian2 = self.calculator2.get_hessian(atoms, coords, **prepare_kwargs)[
            "hessian"
        ]
        hessian_mean = (hessian1 + hessian2) / 2
        hessian_proj = ciq.P.dot(hessian_mean).dot(ciq.P)

        return {
            "energy": ciq.energy,
            "forces": ciq.forces,
            "hessian": hessian_proj,
        }

# [1] https://doi.org/10.1063/1.5090303
#     Geodesic interpolation for reaction pathways
#     Zhu, Thompson, Martinez, 2019

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.constants import BOHR2ANG

try:
    from geodesic_interpolate.interpolation import redistribute
    from geodesic_interpolate.geodesic import Geodesic as Geo

    can_geodesic = True
except (ImportError, ModuleNotFoundError):
    can_geodesic = False


class Geodesic(Interpolator):
    def __init__(
        self,
        *args,
        align: bool = False,
        tol: float = 2e-3,
        scaling: float = 1.7,
        dist_cutoff: float = 3.0,
        friction: float = 1e-2,
        maxiter: int = 15,
        microiter: int = 20,
        **kwargs,
    ):
        """Geodesic Interpolation for Reaction Pathways.

        Requires the 'geodesic-interpolate' package found at

            https://github.com/virtualzx-nad/geodesic-interpolate.git

        Parameters
        ----------
        align : bool, optional
            Whether to align geometries, by default False
        tol : float, optional
            Convergence tolerance, by default 2e-3
        scaling : float, optional
            Exponential parameter for morse potential, by default 1.7
        dist_cutoff : float, optional
            Cut-off value for the distance between a pair of atoms
            to be included in the coordinate system, by default 3.0
        friction : float, optional
            Size of friction term used to prevent very large
            change of geometry, by default 1e-2
        maxiter : int, optional
            Maximum number of minimization iterations, by default 15
        microiter : int, optional
            Maximum number of micro iterations for
            sweeping algorithm, by default 20
        """
        super().__init__(*args, align=align, **kwargs)
        self.tol = tol
        self.scaling = scaling
        self.dist_cutoff = dist_cutoff
        self.friction = friction
        self.maxiter = maxiter
        self.microiter = microiter

    def interpolate(self, initial_geom, final_geom, **kwargs):
        if not can_geodesic:
            raise ModuleNotFoundError(
                "Geodesic interpolation requires the geodesic_interpolate package."
            )
        coords3d = np.array((initial_geom.coords3d, final_geom.coords3d))
        coords3d *= BOHR2ANG

        nimages = 2 + self.between
        raw = redistribute(self.atoms, coords3d, nimages, tol=self.tol * 5)

        smoother = Geo(
            self.atoms,
            raw,
            self.scaling,
            threshold=self.dist_cutoff,
            friction=self.friction,
        )
        if len(self.atoms) > 35:
            path = smoother.sweep(
                tol=self.tol, max_iter=self.maxiter, micro_iter=self.microiter
            )
        else:
            path = smoother.smooth(tol=self.tol, max_iter=self.maxiter)
        path /= BOHR2ANG
        interpolated_geoms = [Geometry(self.atoms, coords) for coords in path]
        return interpolated_geoms[1:-1]

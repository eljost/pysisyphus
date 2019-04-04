#!/usr/bin/env python3

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize

from pysisyphus.calculators.LSTPot import LSTPot
from pysisyphus.constants import BOHR2ANG, ANG2BOHR
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import fit_rigid, geom_from_xyz_file
from pysisyphus.interpolate.Interpolator import Interpolator


np.set_printoptions(suppress=True, precision=3)

# https://www.sciencedirect.com/science/article/pii/S0927025603001113


class LST(Interpolator):

    def cost_function(self, wa_c, f, rab, wab):
        wa_c = wa_c.reshape(-1, 3)
        rab_c = pdist(wa_c)

        rab_i = rab(f)
        wa_i = wab(f)

        first_term = np.sum(
            ((rab_i - rab_c)**2) / (rab_i**4)
        )
        second_term = 1e-6 * np.sum(
            (wa_i - wa_c)**2
        )
        return first_term + second_term

    def interpolate(self, initial_geom, final_geom):
        coords3d = np.array((initial_geom.coords3d, final_geom.coords3d))
        # coords3d *= BOHR2ANG

        # Calculate the condensed distances matrices
        pdists = [pdist(c3d) for c3d in coords3d]

        def rab_(f, pdist_r, pdist_p):
            return (1-f)*pdist_r + f*pdist_p
        rab = lambda f: rab_(f, pdists[0], pdists[1])

        def wab_(f, coords_r, coords_p):
            return (1-f)*coords_r + f*coords_p
        wab = lambda f: wab_(f, coords3d[0], coords3d[1])
        G = lambda w_c, f: self.cost_function(w_c, f, rab, wab)

        interpolated_geoms = list()
        x0_flat = wab(0).flatten()
        for f in np.linspace(0, 1, self.between):
            x0_flat = wab(f)
            res = minimize(self.cost_function, x0=x0_flat, args=(f, rab, wab), method="BFGS",
                           options={
                               "gtol": 1e-4,
                            }
            )
            x0_flat = res.x
            print(f"{f:.04f}, success: {res.success}")
            interpolated_geoms.append(
                # Geometry(self.atoms, res.x*ANG2BOHR)
                Geometry(self.atoms, res.x)
            )
        return interpolated_geoms

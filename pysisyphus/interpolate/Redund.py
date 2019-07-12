#!/usr/bin/env python3

import itertools as it

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.intcoords.helpers import get_tangent, form_coordinate_union
from pysisyphus.xyzloader import write_geoms_to_trj


class Redund(Interpolator):

    def __init__(self, geoms, between, align=True):
        super().__init__(geoms, between, align)

        self.geoms = [Geometry(geom.atoms, geom.cart_coords, coord_type="redund")
                      for geom in self.geoms]

    def interpolate(self, initial_geom, final_geom):
        print(f"No. of primitives at initial structure: {initial_geom.coords.size}")
        print(f"No. of primitives at final structure: {final_geom.coords.size}")

        prim_indices = form_coordinate_union(initial_geom, final_geom)
        union_length = len(list(it.chain(*prim_indices)))
        print("Union of primitives: ", union_length)

        geom1 = Geometry(initial_geom.atoms, initial_geom.cart_coords,
                         coord_type="redund", prim_indices=prim_indices
        )
        geom2 = Geometry(final_geom.atoms, final_geom.cart_coords,
                         coord_type="redund", prim_indices=prim_indices
        )

        dihed_start = geom1.internal.dihed_start
        initial_tangent = get_tangent(geom1.coords, geom2.coords, dihed_start)
        initial_diff = np.linalg.norm(initial_tangent)
        approx_stepsize = initial_diff / (self.between+1)
        final_prims = geom2.internal.prim_coords

        geoms = [geom1, ]
        for i in range(self.between):
            print(f"Interpolating {i+1:03d}/{self.between:03d}")
            new_geom = geoms[-1].copy()
            prim_tangent = get_tangent(new_geom.coords, final_prims, dihed_start)
            # Form active set
            B = new_geom.internal.B_prim
            G = B.dot(B.T)
            eigvals, eigvectors = np.linalg.eigh(G)
            U = eigvectors[:, np.abs(eigvals) > 1e-6]
            reduced_tangent = (np.einsum("i,ij->j", prim_tangent, U) * U).sum(axis=1)
            reduced_tangent /= np.linalg.norm(reduced_tangent)
            step = approx_stepsize * reduced_tangent
            try:
                new_coords = new_geom.coords + step
            except ValueError as err:
                fn = "redund_interpol_fail.trj"
                write_geoms_to_trj(geoms, fn)
                print(f"Chosen set of primitives is not valid for step {i} "
                      f"of {self.between}. Wrote interpolation progress to "
                       "'{fn}'."
                )
                raise err

            new_geom.coords = new_coords
            geoms.append(new_geom)
        return geoms[1:]

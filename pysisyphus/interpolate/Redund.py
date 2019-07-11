#!/usr/bin/env python3

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.xyzloader import write_geoms_to_trj


class Redund(Interpolator):

    def __init__(self, geoms, between, align=True):
        super().__init__(geoms, between, align)

        self.geoms = [Geometry(geom.atoms, geom.coords, coord_type="redund")
                      for geom in self.geoms]

    def to_set(self, iterable):
        return set([tuple(_) for _ in iterable])

    def get_ind_sets(self, geom):
        bonds, bends, dihedrals = geom.internal.prim_indices
        return self.to_set(bonds), self.to_set(bends), self.to_set(dihedrals)

    def merge_coordinate_definitions(self, geom1, geom2):
        bonds1, bends1, dihedrals1 = self.get_ind_sets(geom1)
        bonds2, bends2, dihedrals2 = self.get_ind_sets(geom2)
        # Form new superset of coordinate definitions that contain
        # all definitions from geom1 and geom2.
        all_bonds = bonds1 | bonds2
        all_bends = bends1 | bends2
        all_dihedrals = dihedrals1 | dihedrals2
        all_prim_indices = (all_bonds, all_bends, all_dihedrals)
        # Check if internal coordinates that are only present in one
        # of the two geometries are valid in the other. If not we omit
        # this coordinate definition in the end.
        redundant = RedundantCoords(geom1.atoms, geom1.cart_coords,
                                    prim_indices=all_prim_indices)
        bonds, bends, dihedrals = redundant.prim_indices
        return self.to_set(bonds), self.to_set(bends), self.to_set(dihedrals)

    def interpolate(self, initial_geom, final_geom):
        print(f"No. of primitives at initial structure: {initial_geom.coords.size}")
        print(f"No. of primitives at final structure: {final_geom.coords.size}")

        bonds1, bends1, dihedrals1 = self.merge_coordinate_definitions(initial_geom, final_geom)
        bonds2, bends2, dihedrals2 = self.merge_coordinate_definitions(final_geom, initial_geom)

        # Only use primitive coordinate definitions that are valid for both
        valid_bonds = bonds1 & bonds2
        valid_bends = bends1 & bends2
        valid_dihedrals = dihedrals1 & dihedrals2
        prim_indices = (valid_bonds, valid_bends, valid_dihedrals)
        print("Union of primitives",
              len(valid_bonds) + len(valid_bends) + len(valid_dihedrals))

        geom1 = Geometry(initial_geom.atoms, initial_geom.cart_coords,
                         coord_type="redund", prim_indices=prim_indices
        )
        geom2 = Geometry(final_geom.atoms, final_geom.cart_coords,
                         coord_type="redund", prim_indices=prim_indices
        )

        def get_tangent(prims1, prims2, bonds_bends):
            diff = prims2 - prims1
            diheds = diff[bonds_bends:].copy()
            diheds_plus = diheds.copy() + 2*np.pi
            diheds_minus = diheds.copy() - 2*np.pi
            bigger = np.abs(diheds) > np.abs(diheds_plus)
            diheds[bigger] = diheds_plus[bigger]
            bigger = np.abs(diheds) > np.abs(diheds_minus)
            diheds[bigger] = diheds_minus[bigger]
            diff[bonds_bends:] = diheds
            return diff

        bonds_bends = len(valid_bonds) + len(valid_bends)
        initial_tangent = get_tangent(geom1.coords, geom2.coords, bonds_bends)
        initial_diff = np.linalg.norm(initial_tangent)
        approx_stepsize = initial_diff / (self.between+1)
        final_prims = geom2.internal.prim_coords

        geoms = [geom1, ]
        for i in range(self.between):
            print(f"Interpolating {i+1:03d}/{self.between:03d}")
            new_geom = geoms[-1].copy()
            prim_tangent = get_tangent(new_geom.coords, final_prims, bonds_bends)
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

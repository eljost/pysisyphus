import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.intcoords.exceptions import DifferentPrimitivesException
from pysisyphus.intcoords.helpers import get_tangent, form_coordinate_union
from pysisyphus.xyzloader import write_geoms_to_trj


class Redund(Interpolator):
    def __init__(self, geoms, between, align=True):
        super().__init__(geoms, between, align)

        self.geoms = [
            Geometry(geom.atoms, geom.cart_coords, coord_type="redund")
            for geom in self.geoms
        ]

    def dump_progress(self, geoms, out_fn="redund_interpol_fail.trj"):
        write_geoms_to_trj(geoms, out_fn)
        print(f"Dumped interpolation progress to '{out_fn}'.")

    def interpolate(self, initial_geom, final_geom, typed_prims=None):
        print(f"No. of primitives at initial structure: {initial_geom.coords.size}")
        print(f"No. of primitives at final structure: {final_geom.coords.size}")

        if typed_prims is None:
            typed_prims = form_coordinate_union(initial_geom, final_geom)
            print(f"Union of primitives: {len(typed_prims)}")
        else:
            print(f"Using supplied primitive internals ({len(typed_prims)}).")

        geom1 = Geometry(
            initial_geom.atoms,
            initial_geom.cart_coords,
            coord_type="redund",
            coord_kwargs={
                "typed_prims": typed_prims,
            },
        )
        geom2 = Geometry(
            final_geom.atoms,
            final_geom.cart_coords,
            coord_type="redund",
            coord_kwargs={
                "typed_prims": typed_prims,
            },
        )

        dihedral_indices = geom1.internal.dihedral_indices
        initial_tangent = get_tangent(geom1.coords, geom2.coords, dihedral_indices)
        initial_diff = np.linalg.norm(initial_tangent)
        approx_stepsize = initial_diff / (self.between + 1)
        final_prims = geom2.internal.prim_coords

        geoms = [
            geom1,
        ]
        for i in range(self.between):
            print(f"Interpolating {i+1:03d}/{self.between:03d}")
            new_geom = geoms[-1].copy()
            try:
                prim_tangent = get_tangent(new_geom.coords, final_prims, dihedral_indices)
            # This will result in a different number of internals at the two outer and
            # the inner geometries.
            except DifferentPrimitivesException:
                new_typed_prims = new_geom.internal.typed_prims
                new_num = len(new_typed_prims)
                print(
                    f"Encountered different number of primitive internals ({new_num} vs. "
                    f"{len(typed_prims)}) at geometry {i+1}!\n"
                    "Restarting interpolation with reduced number of internals "
                    "at all geometries."
                )
                self.dump_progress(geoms, out_fn=f"redund_interpol_fail_at_{i+1:03d}.trj")
                print()
                # Recursive call with reduced set of primitive internals
                return self.interpolate(
                    initial_geom, final_geom, typed_prims=new_typed_prims
                )

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
                self.dump_progress(geoms)
                raise err

            new_geom.coords = new_coords
            geoms.append(new_geom)
        print()
        return geoms[1:]

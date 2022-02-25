import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.intcoords.exceptions import (
    DifferentPrimitivesException,
    RebuiltInternalsException,
)
from pysisyphus.intcoords.helpers import get_tangent, form_coordinate_union
from pysisyphus.xyzloader import write_geoms_to_trj


class Redund(Interpolator):
    def __init__(self, *args, align=True, **kwargs):
        super().__init__(*args, align=align, **kwargs)

        self.geoms = [
            Geometry(geom.atoms, geom.cart_coords, coord_type="redund")
            for geom in self.geoms
        ]

    def dump_progress(self, geoms, out_fn="redund_interpol_fail.trj"):
        write_geoms_to_trj(geoms, out_fn)
        print(f"Dumped interpolation progress to '{out_fn}'.")

    def interpolate(
        self, initial_geom, final_geom, interpolate_only=0, extrapolate=False,
        typed_prims=None,
    ):
        print(f"No. of primitives at initial structure: {initial_geom.coords.size}")
        print(f"No. of primitives at final structure: {final_geom.coords.size}")

        if typed_prims is None:
            typed_prims = form_coordinate_union(initial_geom, final_geom)
            print(f"Union of primitives: {len(typed_prims)}")
        else:
            print(f"Using supplied primitive internals ({len(typed_prims)}).")

        # Recreate geometries with consistent set of internal coordinates

        geom1 = Geometry(
            initial_geom.atoms,
            initial_geom.cart_coords,
            coord_type="redund",
            coord_kwargs={
                "typed_prims": typed_prims,
                "recalc_B": True,
            },
        )
        geom2 = Geometry(
            final_geom.atoms,
            final_geom.cart_coords,
            coord_type="redund",
            coord_kwargs={
                "typed_prims": typed_prims,
                "recalc_B": True,
            },
        )

        dihedral_indices = geom1.internal.dihedral_indices
        initial_tangent = get_tangent(geom1.coords, geom2.coords, dihedral_indices)
        initial_diff = np.linalg.norm(initial_tangent)
        approx_step_size = initial_diff / (self.between + 1)
        final_prims = geom2.internal.prim_coords

        geoms = [
            geom1,
        ]

        def restart(new_geom):
            return self.restart_interpolate(initial_geom, final_geom, geoms, new_geom)

        interpolations = interpolate_only if interpolate_only else self.between
        for i in range(interpolations):
            print(f"Interpolating {i+1:03d}/{self.between:03d}")
            new_geom = geoms[-1].copy()
            # Try to use the target primtive internals (final_prims) to calculate
            # a tangent at the current, new geometry. As some primitives may be
            # invalid at 'new_geom', DifferentPrimitivesException may be raised.
            try:
                prim_tangent = get_tangent(
                    new_geom.coords, final_prims, dihedral_indices
                )
            # If this is raised we update our target primitives and try to continue
            # with them.
            except DifferentPrimitivesException:
                # return self.restart_interpolate(initial_geom, final_geom, new_geom)
                return restart(new_geom)
            if extrapolate:
                prim_tangent *= -1
                approx_step_size *= self.extrapolate_damp

            step = self.step_along_tangent(new_geom, prim_tangent, approx_step_size)
            try:
                new_coords = new_geom.coords + step
            # Will be raised when the sizes of both vectors differ
            except ValueError:
                # return self.restart_interpolate(initial_geom, final_geom, new_geom)
                return restart(new_geom)

            # The current primitives may also break down when a step along the
            # tangent is taken.
            try:
                new_geom.coords = new_coords
            except RebuiltInternalsException:
                # return self.restart_interpolate(initial_geom, final_geom, new_geom)
                return restart(new_geom)

            geoms.append(new_geom)

        print()
        return geoms[1:]

    def restart_interpolate(self, initial_geom, final_geom, geoms, new_geom):
        new_typed_prims = new_geom.internal.typed_prims
        print(
            f"Encountered breakdown of current primitive internals.\n"
            "Restarting interpolation with reduced number of internals."
        )
        self.dump_progress(geoms, out_fn=f"redund_interpol_fail.trj")
        print()
        # Recursive call with reduced set of primitive internals
        return self.interpolate(initial_geom, final_geom, typed_prims=new_typed_prims)

    def step_along_tangent(self, geom, prim_tangent, step_size):
        # Form active set
        B = geom.internal.B_prim
        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)
        U = eigvectors[:, np.abs(eigvals) > 1e-6]
        reduced_tangent = (np.einsum("i,ij->j", prim_tangent, U) * U).sum(axis=1)
        reduced_tangent /= np.linalg.norm(reduced_tangent)
        step = step_size * reduced_tangent
        return step

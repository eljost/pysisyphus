import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import align_geoms
from pysisyphus.xyzloader import write_geoms_to_trj


class Interpolator:
    def __init__(
        self,
        geoms,
        between,
        extrapolate=0,
        extrapolate_before=0,
        extrapolate_after=0,
        extrapolate_damp=1.0,
        align=False,
    ):
        self.geoms = geoms
        self.between = between
        self.extrapolate = extrapolate
        self.extrapolate_before = (
            extrapolate_before if extrapolate_before else self.extrapolate
        )
        self.extrapolate_after = (
            extrapolate_after if extrapolate_after else self.extrapolate
        )
        self.extrapolate_damp = extrapolate_damp
        one_atom_geom = any([len(geom.atoms) == 1 for geom in geoms])
        # Don't try to align one atom species
        self.align = align and not one_atom_geom

        assert len(geoms) >= 2, "Need at least two geometries to interpolate!"

        # Check for consistent atom ordering
        for i, geom in enumerate(self.geoms[:-1]):
            next_geom = self.geoms[i + 1]
            atoms = geom.atoms
            next_atoms = next_geom.atoms
            assert len(atoms) == len(
                next_atoms
            ), f"Geometries {i} and {i+1} have a different number of atoms!"
            assert (
                atoms == next_atoms
            ), f"Different atom ordering in geometries {i} and {i+1}!"

        self.atoms = self.geoms[0].atoms
        if self.align:
            align_geoms(geoms)
        self.all_geoms = None

    def interpolate_all(self):
        all_geoms = list()

        if self.extrapolate_before:
            geoms_extrapol_before = self.interpolate(
                self.geoms[0],
                self.geoms[1],
                interpolate_only=self.extrapolate_before,
                extrapolate=True,
            )
            all_geoms += geoms_extrapol_before[::-1]
        # Interpolate between all pairs of geometries
        for i, initial_geom in enumerate(self.geoms[:-1]):
            final_geom = self.geoms[i + 1]
            interpolated_geoms = self.interpolate(initial_geom, final_geom)
            assert len(interpolated_geoms) == self.between, (
                "Something is wrong with the number of interpolated "
                "geometries. Maybe you accidentally also return the "
                "initial and final geometry?"
            )
            all_geoms.append(self.geoms[i])
            all_geoms.extend(interpolated_geoms)
        # As we only added the i-th geometry and the new interpolated
        # geometries of every geometry pair we also have to add the last
        # ((i+1)-th) geometry at the end.
        all_geoms.append(self.geoms[-1])

        if self.extrapolate_after:
            geoms_extrapol_after = self.interpolate(
                self.geoms[-1],
                self.geoms[-2],
                interpolate_only=self.extrapolate_after,
                extrapolate=True,
            )
            all_geoms += geoms_extrapol_after

        self.all_geoms = all_geoms
        if self.align:
            align_geoms(self.all_geoms)
        return all_geoms

    def interpolate(
        self, initial_geom, final_geom, interpolate_only=0, extrapolate=False
    ):
        initial_coords = initial_geom.coords

        # Linear interpolation
        step = (final_geom.coords - initial_geom.coords) / (self.between + 1)
        if extrapolate:
            step *= -1
            step *= self.extrapolate_damp
        # When we extrapolate we probably want a number of geometries that is
        # different from self.between
        interpolations = interpolate_only if interpolate_only else self.between
        # initial + i*step
        i_array = np.arange(1, interpolations + 1)
        new_coords = initial_coords + i_array[:, None] * step
        return [Geometry(self.atoms, nc) for nc in new_coords]

    def all_geoms_to_trj(self, trj_fn):
        write_geoms_to_trj(self.all_geoms, trj_fn)

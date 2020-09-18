import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import align_geoms
from pysisyphus.xyzloader import write_geoms_to_trj


class Interpolator:

    def __init__(self, geoms, between, align=False):
        self.geoms = geoms
        self.between = between
        self.align = align

        assert len(geoms) >= 2, "Need at least two geometries to interpolate!"

        # Check for consistent atom ordering
        for i, geom in enumerate(self.geoms[:-1]):
            next_geom = self.geoms[i+1]
            atoms = geom.atoms
            next_atoms = next_geom.atoms
            assert len(atoms) == len(next_atoms), \
                f"Geometries {i} and {i+1} have a different number of atoms!"
            assert atoms == next_atoms, \
                f"Different atom ordering in geometries {i} and {i+1}!"

        self.atoms = self.geoms[0].atoms
        if self.align:
            align_geoms(geoms)
        self.all_geoms = None

    def interpolate_all(self):
        all_geoms = list()
        # Interpolate between all pairs of geometries
        for i, initial_geom in enumerate(self.geoms[:-1]):
            final_geom = self.geoms[i+1]
            interpolated_geoms = self.interpolate(initial_geom, final_geom)
            assert len(interpolated_geoms) == self.between, \
                "Something is wrong with the number of interpolated " \
                "geometries. Maybe you accidentally also return the " \
                "initial and final geometry?"
            all_geoms.append(self.geoms[i])
            all_geoms.extend(interpolated_geoms)
        # As we only added the i-th geometry and the new interpolated
        # geometries of every geometry pair we also have to add the last
        # ((i+1)-th) geometry at the end.
        all_geoms.append(self.geoms[-1])
        self.all_geoms = all_geoms
        if self.align:
            align_geoms(self.all_geoms)
        return all_geoms

    def interpolate(self, initial_geom, final_geom):
        initial_coords = initial_geom.coords
        final_coords = final_geom.coords

        # Linear interpolation
        step = (final_coords-initial_coords) / (self.between+1)
        # initial + i*step
        i_array = np.arange(1, self.between+1)
        new_coords = initial_coords + i_array[:, None]*step
        return [Geometry(self.atoms, nc) for nc in new_coords]
    
    def all_geoms_to_trj(self, trj_fn):
        write_geoms_to_trj(self.all_geoms, trj_fn)

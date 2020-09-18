import numpy as np
import rmsd

from pysisyphus.stocastic.Pipeline import Pipeline


class Kick(Pipeline):

    def __init__(self, geom, radius=0.5, **kwargs):
        super().__init__(geom, **kwargs)

        self.radius = radius

    def get_kick(self):
        # Interval [0, 1)
        kick = np.random.random(4*self.coords_size)
        # Stretch [0, 1) to  [-r ... r)
        kick = self.radius * (2*kick - 1)
        kick = kick.reshape(-1, 3)
        # Filter for kicks within a sphere with radius self.radius
        kick_lengths = np.linalg.norm(kick, axis=1)
        valid_kicks = kick[kick_lengths <= self.radius].flatten()
        # Do a recursion when not enough valid kicks were found
        if valid_kicks.size < self.coords_size:
            valid_kicks = self.get_kick()
        # Don't return more than we need
        return valid_kicks[:self.coords_size]

    def get_input_geom(self, geom):
        kick = self.get_kick()
        new_geom = geom.copy()
        new_coords = new_geom.coords + kick
        # Rotate the newly generated coordinates on the initial
        # coordinates.
        # TODO: this may be not needed as we align+match later on ...
        new_coords = rmsd.kabsch_rotate(new_coords.reshape(-1, 3),
                                        self.initial_coords3d
        ).flatten()
        new_geom.coords = new_coords
        return new_geom

    def run_cycle(self, geom):
        print(f"##### Cycle {self.cur_cycle:03d}, "
              f"Micro Cycle {self.cur_micro_cycle:03d} #####")
        opt_geoms = [self.run_kicked_geom(geom) for i in range(self.cycle_size)]
        # Filter out None
        opt_geoms = [geom for geom in opt_geoms if geom]
        opt_num = len(opt_geoms)
        print(f"{opt_num}/{self.cycle_size} optimizations converged.")

        # Comparing to the initial geometry is only useful when the initial
        # geometry is optimized. Otherwise all (small) kicks will converge
        # to the same optimized structure, that is still very different from
        # the initial one.
        """
        initial_rmsds = [
            rmsd.kabsch_rmsd(self.initial_coords3d, ogeom.coords3d)
            for ogeom in opt_geoms
        ]
        print("Initial RMDS")
        print(initial_rmsds)
        """
        kept_geoms = self.get_unique_geometries(opt_geoms)

        # cycle_str = f"{self.cur_cycle:03d}_{self.cur_micro_cycle:03d}"
        # fn_base = f"cycle_{cycle_str}"
        # trj_fn = f"{fn_base}.trj"
        # with open(trj_fn, "w") as handle:
            # handle.write(make_trj_str_from_geoms(opt_geoms))

        # trj_filtered_fn = f"{fn_base}_filtered.trj"
        # with open(trj_filtered_fn, "w") as handle:
            # handle.write(make_trj_str_from_geoms(kept_geoms))

        # for i, ogeom in enumerate(kept_geoms):
            # fn = f"geom_{i:02d}_{cycle_str}.xyz"
            # with open(fn, "w") as handle:
                # handle.write(ogeom.as_xyz())
        self.cur_micro_cycle += 1

        return kept_geoms

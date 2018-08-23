#!/usr/bin/env python3

import time

import numpy as np
from scipy.spatial.distance import pdist

from pysisyphus.InternalCoordinates import get_cov_radii_sum_array


class Pipeline:

    def __init__(self, geom, seed=None, cycles=5, cycle_size=15,
                 rmsd_thresh=.1, energy_thresh=1e-3):
        self.initial_geom = geom
        self.cycles = cycles
        self.cycle_size = cycle_size
        if seed is None:
            self.seed = int(time.time())
        else:
            self.seed = seed
        self.rmsd_thresh = rmsd_thresh
        self.energy_thresh = energy_thresh

        np.random.seed(self.seed)
        print(f"Seed: {self.seed}")
        self.coords_size = self.initial_geom.coords.size

        self.calc_counter = 0
        self.cur_cycle = 0
        self.cur_micro_cycle = 0

        self.new_geoms = []
        self.new_energies = []

        self.initial_coords3d = self.initial_geom.copy()
        self.atoms = self.initial_geom.atoms

    def __str__(self):
        return f"Pipeline(seed={self.seed})"

    def atoms_are_too_close(self, geom, factor=.7):
        """Determine if atoms are too close."""
        dist_mat = pdist(geom.coords3d)
        cov_rad_mat = get_cov_radii_sum_array(geom.atoms, geom.coords)
        too_close = dist_mat < factor*cov_rad_mat
        return any(too_close)

    def geometry_is_new(self, geom):
        """Determine if geometry is not already known."""
        return True
    
    def get_input_geom(self):
        raise Exception("Implement me!")

    """
    def run_new_geom(self, geom):
        input_geom = self.get_input_geom(geom)
        input_coords = input_geom.coords.copy()

        # Check if the geometry is similar to an already known starting
        # geometry.
        overlaps = new_coords.dot(np.array(self.starting_coords).T)/self.coords_size
        # print("overlaps with already known starting coordinates")
        # print(overlaps)
        max_overlap_ind = overlaps.argmax()
        max_overlap = overlaps[max_overlap_ind]
        similar_fn = f"similar_{self.similar_ind:03d}.trj"
        # print(f"max overlap is {max_overlap:.1f}, {similar_fn}, index {max_overlap_ind}")
        max_coords = self.starting_coords[max_overlap_ind]
        self.similar_ind += 1
        rmsds = list()
        for sc in self.starting_coords:
            sc3d = sc.reshape(-1, 3)
            rm = rmsd.kabsch_rmsd(new_coords.reshape(-1,3), sc3d)
            rmsds.append(rm)
        rmsds = np.array(rmsds)
        self.starting_coords.append(new_coords)
    """

    def get_unique_geometries(self, geoms):
        geom_num = len(geoms)
        rmsds = np.full((geom_num, geom_num), np.inf)
        for i, j in it.combinations(range(geom_num), 2):
            coords1 = geoms[i].coords.reshape(-1, 3)
            coords2 = geoms[j].coords.reshape(-1, 3)
            rmsds[i, j] = rmsd.kabsch_rmsd(coords1, coords2)
        is_, js = np.where(rmsds < self.rmsd_thresh)
        similar_inds = np.unique(js)
        all_inds = np.arange(geom_num)
        unique_geoms = [geoms[i] for i in range(geom_num) if not i in similar_inds]
        unique_num = len(unique_geoms)
        return unique_geoms

    def run_geom_opt(self, geom):
        calc = XTB(calc_number=self.calc_counter)
        self.calc_counter += 1
        opt_geom = calc.run_opt(geom.atoms, geom.coords, keep=False)
        return opt_geom

    def run(self):
        while self.cur_cycle < self.cycles:
            print(f"Starting cycle {self.cur_cycle} with " \
                  f"{len(self.geoms_to_kick)} geometries.")
            new_geoms = list(
                it.chain(
                    *[self.run_cycle(geom) for geom in self.geoms_to_kick]
                )
            )
            new_num = len(new_geoms)
            print(f"Kicks in cycle {self.cur_cycle} produced "
                  f"{new_num} new geometries.")
            kept_geoms = self.get_unique_geometries(new_geoms)
            kept_num = len(kept_geoms)
            self.geoms_to_kick = kept_geoms
            self.cur_cycle += 1
            self.cur_micro_cycle = 0
            print()
            geoms_sorted = sorted(kept_geoms, key=lambda g: g._energy)
            trj_filtered_fn = f"cycle_{self.cur_cycle:03d}.trj"
            with open(trj_filtered_fn, "w") as handle:
                handle.write(make_trj_str_from_geoms(geoms_sorted))
        # keine optimierungen von bereits bekannten startgeometrien starten
        #                                           endgeometrien starten
        # energie einbeziehen

    def geoms_to_trj(self, geoms, fn):
        with open(fn, "w") as handle:
            handle.write(make_trj_str_from_geoms(geoms))

    def run(self):
        while self.cur_cycle < self.cycles:
            print(f"Cycle {self.cur_cycle}")
            input_geoms = [self.get_input_geom() for _ in range(self.cycle_size)]
            # Write input geometries to disk
            self.geoms_to_trj(input_geoms, f"cycle_{self.cur_cycle:03d}_input.trj")
            opt_geoms = [self.run_geom_opt(geom) for geom in input_geoms]

            kept_geoms = list()
            rejected_geoms = 0
            for geom in opt_geoms:
                # Filter out None and reject geoms where atoms are too close
                # or geometries already known
                if (geom is None
                    or self.atoms_are_too_close(geom)
                    or not self.geom_is_new(geom)):
                    continue

                energy = geom.energy
                i = bisect.bisect_left(self.new_energies, energy)
                energies_arr = np.array(self.new_energies)
                diffs = np.abs(energies_arr - energy)
                if len(self.new_energies) > 0 and diffs.min() < self.energy_thresh:
                    #print("min(diffs) is", diffs.min())
                    rejected_geoms += 1
                    continue
                self.new_energies.insert(i, energy)
                self.new_geoms.insert(i, geom)
                kept_geoms.append(geom)

            kept_num = len(kept_geoms)
            print(f"Kicks in cycle {self.cur_cycle} produced "
                  f"{kept_num} new geometries.")
            if rejected_geoms:
                print(f"{rejected_geoms} geometries were rejected by energy criteria.")

            trj_filtered_fn = f"cycle_{self.cur_cycle:03d}.trj"
            # Sort by energy
            kept_geoms = sorted(kept_geoms, key=lambda g: g.energy)
            self.geoms_to_trj(kept_geoms, trj_filtered_fn)
            print()

            # diffs = np.diff(self.new_energies)
            # print(f"min(diffs) {diffs.min():.4f}")

            self.cur_cycle += 1
        fn = "final.trj"
        self.geoms_to_trj(self.new_geoms, fn)


if __name__ == "__main__":
    from pysisyphus.helpers import geom_from_library
    geom = geom_from_library("benzene.xyz")
    p = Pipeline(geom)
    print(p)

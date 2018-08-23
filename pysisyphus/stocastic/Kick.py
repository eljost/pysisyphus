#!/usr/bin/env python3

import itertools as it
import time

import numpy as np
import rmsd
from scipy.spatial.distance import pdist, squareform

from pysisyphus.calculators.XTB import XTB
from pysisyphus.Geometry import Geometry
from pysisyphus.InternalCoordinates import get_cov_radii_sum_array
from pysisyphus.xyzloader import make_trj_str_from_geoms


np.set_printoptions(suppress=True, precision=2)

class Kick:

    def __init__(self, geom, radius=0.5, cycles=5, cycle_size=15,
                 seed=None, rmsd_thresh=1):
        self.initial_geom = geom
        self.radius = radius
        self.cycles = cycles
        self.cycle_size = cycle_size
        if seed is None:
            self.seed = int(time.time())
        else:
            self.seed = seed
        self.rmsd_thresh = rmsd_thresh

        np.random.seed(self.seed)
        print("seed", self.seed)
        self.coords_size = self.initial_geom.coords.size
        self.calc_counter = 0
        self.cur_cycle = 0
        self.cur_micro_cycle = 0

        self.initial_coords3d = geom.coords3d
        self.atoms = self.initial_geom.atoms
        self.geoms_to_kick = [geom]
        self.kept_geoms = [geom]
        self.starting_coords = [self.initial_geom.coords]

        self.similar_ind = 0

    def get_kick(self):
        # Interval [0, 1)
        kick = np.random.random(4*self.coords_size)
        # Stretch [0, 1) to  [-r ... r)
        kick = self.radius * (2*kick - 1)
        kick = kick.reshape(-1, 3)
        # Filter for kicks within a sphere with radius self.radius
        kick_lengths = np.linalg.norm(kick, axis=1)
        valid_kicks = kick[kick_lengths <= self.radius].flatten()
        if valid_kicks.size < self.coords_size:
            print("not enough kicks, recursing")
            valid_kicks = self.get_kick()
        # Don't return more than we need
        return valid_kicks[:self.coords_size]

    def get_kicked_geom(self, geom):
        kick = self.get_kick()
        new_geom = Geometry(geom.atoms, geom.coords)
        new_coords = new_geom.coords + kick
        # Rotate the newly generated coordinates on the initial
        # coordinates.
        new_coords = rmsd.kabsch_rotate(new_coords.reshape(-1, 3),
                                        self.initial_coords3d
        ).flatten()
        new_geom.coords = new_coords
        return new_geom

    def run_kicked_geom(self, geom):
        new_geom = self.get_kicked_geom(geom)
        new_coords = new_geom.coords.copy()
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
        # with open(similar_fn, "w") as handle:
            # ovlp_geom = Geometry(self.atoms, max_coords)
            # handle.write(make_trj_str_from_geoms((new_geom, ovlp_geom)))
        self.similar_ind += 1
        rmsds = list()
        for sc in self.starting_coords:
            sc3d = sc.reshape(-1, 3)
            rm = rmsd.kabsch_rmsd(new_coords.reshape(-1,3), sc3d)
            rmsds.append(rm)
        rmsds = np.array(rmsds)
        # print(f"RMSDs, min_rmsds with {rmsds.argmin()}")
        # print(rmsds)
        # print("Quotient")
        # print(rmsds/overlaps)
        # import pdb; pdb.set_trace()
        self.starting_coords.append(new_coords)
        calc = XTB(calc_number=self.calc_counter)
        self.calc_counter += 1
        opt_geom = calc.run_opt(new_geom.atoms, new_geom.coords, keep=False)
        # print()
        return opt_geom

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
        kept_geoms = [geoms[i] for i in range(geom_num) if not i in similar_inds]
        kept_num = len(kept_geoms)
        # print(rmsds)
        # print("similar inds", similar_inds)
        # print(f"Keeping {kept_num}/{geom_num} geometries.")
        return kept_geoms

    def atoms_are_too_close(self, geom, factor=.7):
        """Determine if atoms are too close."""
        dist_mat = pdist(geom.coords3d)
        cov_rad_mat = get_cov_radii_sum_array(geom.atoms, geom.coords)
        to_reject = dist_mat < factor*cov_rad_mat
        return any(to_reject)

    def run_cycle(self, geom):
        print(f"##### Cycle {self.cur_cycle:03d}, "
              f"Micro Cycle {self.cur_micro_cycle:03d} #####")
        opt_geoms = [self.run_kicked_geom(geom) for i in range(self.cycle_size)]
        # Filter out None
        opt_geoms = [geom for geom in opt_geoms if geom]
        opt_num = len(opt_geoms)
        energies = np.array([geom.energy for geom in opt_geoms])
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

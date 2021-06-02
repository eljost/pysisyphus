import bisect
import itertools as it
import logging
import time

import numpy as np
import rmsd
from scipy.spatial.distance import pdist

from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import check_for_stop_sign
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.intcoords.setup import get_pair_covalent_radii
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.stocastic.align import matched_rmsd
from pysisyphus.xyzloader import make_trj_str_from_geoms


class Pipeline:

    def __init__(self, geom, calc_getter=None, seed=None, max_cycles=5, cycle_size=15,
                 rmsd_thresh=.1, energy_thresh=1e-3, energy_range=.125,
                 compare_num=25, break_after=2,
                 calc_kwargs=None):
        self.initial_geom = geom
        self.calc_getter = calc_getter
        self.max_cycles = max_cycles
        self.cycle_size = cycle_size
        if seed is None:
            self.seed = int(time.time())
        else:
            self.seed = seed
        self.rmsd_thresh = rmsd_thresh
        self.energy_thresh = energy_thresh
        self.energy_range = energy_range
        self.compare_num = compare_num
        self.break_after = break_after
        self.calc_kwargs = {
                "charge": 0,
                "mult": 1,
        }
        if calc_kwargs is not None:
            self.calc_kwargs.update(calc_kwargs)

        np.random.seed(self.seed)
        self.logger = logging.getLogger("stocastic")

        self.is_analytical2d = len(geom.atoms) == 1
        self.log(f"Seed: {self.seed}")
        self.coords_size = self.initial_geom.coords.size

        self.calc_counter = 0
        self.cur_cycle = 0
        self.cur_micro_cycle = 0
        # Indicates if the next cycle is the last one
        self.break_in = self.break_after

        self.new_geoms = []
        self.new_energies = []

        self.initial_coords3d = self.initial_geom.coords3d
        self.atoms = self.initial_geom.atoms

    def __str__(self):
        return f"{self.__class__.__name__}(seed={self.seed})"

    def log(self, message):
        """Write a log message.

        Wraps the logger variable.

        Parameters
        ----------
        message : str
            Message to be logged.
        """

        self.logger.debug(f"{message}")

    def get_valid_index_set(self, to_intersect):
        return set(range(len(self.new_energies))) & set(to_intersect)

    def atoms_are_too_close(self, geom, factor=.7):
        """Determine if atoms are too close."""
        dist_mat = pdist(geom.coords3d)
        cov_rad_mat = get_pair_covalent_radii(geom.atoms)
        too_close = dist_mat < factor*cov_rad_mat
        return any(too_close)

    def geom_is_close_in_energy(self, geom):
        energy = geom.energy
        i = bisect.bisect_left(self.new_energies, energy)
        # Determine if there are neighbours that are close in energy
        # as we insert left/before the most similary energy the indices
        # of the (to be) neighbours in the current self.new_energies list
        # are i-1 and i.
        valid_inds = self.get_valid_index_set((i-1, i))
        diffs = [abs(self.new_energies[j] - energy) for j in valid_inds]
        return len(diffs) > 0 and min(diffs) < self.energy_thresh

    def geom_is_new(self, geom):
        """Determine if geometry is not already known."""
        if len(self.new_geoms) == 0:
            self.log("Found first geometry!")
            return True

        i = bisect.bisect_left(self.new_energies, geom.energy)
        to_intersect = range(i-self.compare_num, i+self.compare_num)
        valid_inds = np.array(list(self.get_valid_index_set(to_intersect)))
        new_energies = np.array(self.new_energies)[valid_inds]
        # Restrict geometries for RMSD comparison to an energy range
        # around the energy of the geometry to check.
        in_range = np.abs(new_energies - geom.energy) < self.energy_range
        valid_inds = valid_inds[in_range]
        # If this evalutes to True the energy of the current in geometry is
        # quite different and we add the geometry.
        if valid_inds.size == 0:
            # print("Energy of geometry is very different from the remaining "
                  # "ones. Adding geometry!")
            reason = "different energy."
            is_new = True
        # Otherwise check the RMSD values for the remaining geometries that
        # are close in energy.
        else:
            rmsds = [matched_rmsd(geom, self.new_geoms[i])[0] for i in valid_inds]
            rmsds = np.array(rmsds)
            is_new = rmsds.min() > self.rmsd_thresh
            reason = f"different RMSD (min(RMSD) = {rmsds.min():.3f})"

        if is_new:
            self.log(f"Found new geometry based on {reason}")
        return is_new

    def geom_is_valid(self, geom):
        """Filter out geometries that are None, or were the atoms are too close
        or when they are already known."""

        valid = (geom is not None
                 and not self.geom_is_close_in_energy(geom)
        )

        if not self.is_analytical2d:
            valid = (valid
                     and not self.atoms_are_too_close(geom)
                     and self.geom_is_new(geom)
            )
        return valid
    
    def get_input_geom(self, geom):
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
        unique_geoms = [geoms[i] for i in range(geom_num) if i not in similar_inds]
        return unique_geoms

    def run_geom_opt(self, geom):
        if self.calc_getter is not None:
            calc = self.calc_getter(calc_number=self.calc_counter)
            geom.set_calculator(calc)
            opt_kwargs = {
                "gdiis": False,
                "thresh": "gau_loose",
                "overachieve_factor": 2,
                "max_cycles": 75,
            }
            opt = RFOptimizer(geom, **opt_kwargs)
            opt.run()
            opt_geom = geom if opt.is_converged else None
        else:
            calc = XTB(calc_number=self.calc_counter, **self.calc_kwargs)
            opt_result = calc.run_opt(geom.atoms, geom.coords, keep=False)
            opt_geom = opt_result.opt_geom

        self.calc_counter += 1
        return opt_geom

    def write_geoms_to_trj(self, geoms, fn, comments=None):
        with open(fn, "w") as handle:
            handle.write(
                make_trj_str_from_geoms(geoms, comments, energy_comments=True)
            )

    def run(self):
        for self.cur_cycle in range(self.max_cycles):
            cycle_start = time.time()
            self.log(highlight_text(f"Cycle {self.cur_cycle}"))
            input_geoms = [self.get_input_geom(self.initial_geom)
                           for _ in range(self.cycle_size)]
            # Write input geometries to disk
            self.write_geoms_to_trj(input_geoms, f"cycle_{self.cur_cycle:03d}_input.trj")
            # Run optimizations on input geometries
            calc_start = time.time()
            opt_geoms = list()
            for i, geom in enumerate(input_geoms, 1):
                print(f"Optimizing geometry {i:03d}/{self.cycle_size:03d}", end="\r")
                opt_geoms.append(self.run_geom_opt(geom))
            print()
            calc_end = time.time()
            calc_duration = calc_end - calc_start
            self.log(f"Optimizations took {calc_duration:.0f} s.")

            kept_geoms = list()
            for geom in opt_geoms:
                # Do all the filtering and reject all invalid geometries
                if not self.geom_is_valid(geom):
                    continue

                energy = geom.energy
                i = bisect.bisect_left(self.new_energies, energy)
                self.new_energies.insert(i, energy)
                self.new_geoms.insert(i, geom)
                kept_geoms.append(geom)
                if i == 0 and len(self.new_energies) > 1:
                    last_minimum = self.new_energies[1]
                    diff = abs(energy - last_minimum)
                    self.log(f"It is a new global minimum at {energy:.4f} au! "
                             f"Last one was {diff:.4f} au higher.")

            kept_num = len(kept_geoms)

            trj_filtered_fn = f"cycle_{self.cur_cycle:03d}.trj"
            # Sort by energy
            kept_geoms = sorted(kept_geoms, key=lambda g: g.energy)
            if kept_geoms:
                self.write_geoms_to_trj(kept_geoms, trj_filtered_fn)
                self.log(f"Kicks in cycle {self.cur_cycle} produced "
                         f"{kept_num} new geometries.")
                self.break_in = self.break_after
            elif self.break_in == 0:
                self.log("Didn't find any new geometries in the last "
                      f"{self.break_after} cycles. Exiting!")
                break
            else:
                self.log(f"Cycle {self.cur_cycle} produced no new geometries.")
                self.break_in -= 1

            cycle_end = time.time()
            cycle_duration = cycle_end - cycle_start
            self.log(f"Cycle {i} took {cycle_duration:.0f} s.")
            self.log("")
            if check_for_stop_sign():
                break

        self.log(f"Run produced {len(self.new_energies)} geometries!")
        # Return empty list of nothing was found
        if not self.new_energies:
            return []

        fn = "final.trj"
        self.write_geoms_to_trj(self.new_geoms, fn)
        # self.new_energies = np.array(new_energies)
        np.savetxt("energies.dat", self.new_energies)
        first_geom = self.new_geoms[0]
        first_geom.standard_orientation()
        first_geom.energy = self.new_energies[0]

        if self.is_analytical2d:
            return self.new_geoms

        matched_geoms = [first_geom, ]
        for geom, energy in zip(self.new_geoms[1:], self.new_energies):
            rmsd, (_, matched_geom) = matched_rmsd(first_geom, geom)
            matched_geom.energy = energy
            matched_geoms.append(matched_geom)
        fn_matched = "final_matched.trj"
        self.write_geoms_to_trj(matched_geoms, fn_matched)
        return matched_geoms

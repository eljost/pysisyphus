#!/usr/bin/env python3

import abc
import logging
import os
from pathlib import Path
import sys
import textwrap
import time

import numpy as np
import yaml

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.helpers import check_for_stop_sign, highlight_text, get_coords_diffs


class Optimizer(metaclass=abc.ABCMeta):
    CONV_THRESHS = {
        #              max_force, rms_force, max_step, rms_step
        "gau_loose":  (2.5e-3,    1.7e-3,    1.0e-2,   6.7e-3),
        "gau":        (4.5e-4,    3.0e-4,    1.8e-3,   1.2e-3),
        "gau_tight":  (1.5e-5,    1.0e-5,    6.0e-5,   4.0e-5),
        "gau_vtight": (2.0e-6,    1.0e-6,    6.0e-6,   4.0e-6),
        "baker":      (3.0e-4,    2.0e-4,    3.0e-4,   2.0e-4),
    }

    def __init__(self, geometry, thresh="gau_loose", max_step=0.04,
                 rms_force=None, align=False, dump=False, dump_restart=None,
                 prefix="", reparam_thresh=1e-3, overachieve_factor=0.,
                 restart_info=None, check_coord_diffs=True, coord_diff_thresh=0.01,
                 **kwargs):
        self.geometry = geometry

        self.is_cos = issubclass(type(self.geometry), ChainOfStates)

        assert thresh in self.CONV_THRESHS.keys()
        self.thresh = thresh
        self.convergence = self.make_conv_dict(thresh, rms_force)
        self.align = align
        self.dump = dump
        self.dump_restart = dump_restart
        self.prefix = prefix
        self.reparam_thresh = reparam_thresh
        self.overachieve_factor = float(overachieve_factor)
        self.check_coord_diffs = check_coord_diffs
        self.coord_diff_thresh = float(coord_diff_thresh)

        for key, value in self.convergence.items():
            setattr(self, key, value)

        # Setting some default values
        self.resetted = False
        self.max_cycles = 50
        self.max_step = max_step
        self.rel_step_thresh = 1e-3
        self.out_dir = os.getcwd()

        assert(self.max_step > self.rel_step_thresh)

        if self.is_cos:
            image_num = len(self.geometry.moving_indices)
            print(f"Path with {image_num} moving images.")

        # Overwrite default values if they are supplied as kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.out_dir = Path(self.out_dir)
        if not self.out_dir.exists():
            os.mkdir(self.out_dir)

        current_fn = "current_geometries.trj" if self.is_cos else "current_geometry.xyz"
        self.current_fn = self.get_path_for_fn(current_fn)
        final_fn = "final_geometries.trj" if self.is_cos else "final_geometry.xyz"
        self.final_fn = self.get_path_for_fn(final_fn)

        self.logger = logging.getLogger("optimizer")

        # Setting some empty lists as default
        self.list_attrs = "cart_coords coords energies forces steps " \
                          "max_forces rms_forces max_steps rms_steps " \
                          "cycle_times tangents modified_forces".split()
        for la in self.list_attrs:
            setattr(self, la, list())

        self.opt_results_fn = "optimizer_results.yaml"
        self.image_results_fn = "image_results.yaml"
        self.image_results = list()

        if self.dump:
            out_trj_fn = self.get_path_for_fn("optimization.trj")
            self.out_trj_handle= open(out_trj_fn, "w")
        if self.prefix:
            self.log(f"Created optimizer with prefix {self.prefix}")

        self.restarted = False
        self.last_cycle = 0
        if restart_info is not None:
            if isinstance(restart_info, str):
                restart_info = yaml.load(restart_info, Loader=yaml.SafeLoader)
            self.set_restart_info(restart_info)
            self.restarted = True

    def get_path_for_fn(self, fn):
        return self.out_dir / (self.prefix + fn)

    def make_conv_dict(self, key, rms_force=None):
        if not rms_force:
            threshs = self.CONV_THRESHS[key]
        else:
            print( "Deriving convergence threshold from supplied "
                  f"rms_force={rms_force}.")
            threshs = (1.5*rms_force,
                       rms_force,
                       6*rms_force,
                       4*rms_force,
            )
        keys = ["max_force_thresh",
                "rms_force_thresh",
        ]
        # Only used gradient information for CoS optimizations
        if not self.is_cos:
            keys += ["max_step_thresh", "rms_step_thresh"]
        conv_dict = {
            k: v for k, v in zip(keys, threshs)
        }
        return conv_dict

    def log(self, message):
        # self.logger.debug(f"Cycle {self.cur_cycle:03d}, {message}")
        self.logger.debug(message)

    def check_convergence(self, step=None, multiple=1.0, overachieve_factor=None,
                          energy_thresh=1e-6):
        """Check if the current convergence of the optimization
        is equal to or below the required thresholds, or a multiple
        thereof. The latter may be used in initiating the climbing image.
        """

        if step is None:
            step = self.steps[-1]
        if overachieve_factor is None:
            overachieve_factor = self.overachieve_factor

        # When using a ChainOfStates method we are only interested
        # in optimizing the forces perpendicular to the MEP.
        if self.is_cos:
            forces = self.geometry.perpendicular_forces
        elif len(self.modified_forces) == len(self.forces):
            self.log("Using modified forces to determine convergence!")
            forces = self.modified_forces[-1]
        else:
            forces = self.forces[-1]

        # The forces of fixed images may be zero and this may distort the RMS
        # values. So we take into account the number of moving images with
        # non-zero forces vectors.
        if self.is_cos:
            non_zero_elements = (len(self.geometry.moving_indices)
                                 * self.geometry.coords_length
            )
            rms_force = np.sqrt(np.sum(np.square(forces))/non_zero_elements)
            rms_step = np.sqrt(np.sum(np.square(step))/non_zero_elements)
        else:
            rms_force = np.sqrt(np.mean(np.square(forces)))
            rms_step = np.sqrt(np.mean(np.square(step)))

        max_force = np.abs(forces).max()
        max_step = np.abs(step).max()

        self.max_forces.append(max_force)
        self.rms_forces.append(rms_force)
        self.max_steps.append(max_step)
        self.rms_steps.append(rms_step)

        this_cycle = {
            "max_force_thresh": max_force,
            "rms_force_thresh": rms_force,
            "max_step_thresh": max_step,
            "rms_step_thresh": rms_step
        }

        # Check if force convergence is overachieved
        overachieved = False
        if overachieve_factor > 0:
            max_thresh = self.convergence["max_force_thresh"] / overachieve_factor
            rms_thresh = self.convergence["rms_force_thresh"] / overachieve_factor
            max_ = max_force < max_thresh
            rms_ = rms_force < rms_thresh
            overachieved = max_ and rms_
            if max_:
                self.log("max(force) is overachieved")
            if rms_:
                self.log("rms(force) is overachieved")
            if max_ and rms_:
                self.log("Force convergence overachieved!")

        normal_convergence = all(
            [this_cycle[key] <= getattr(self, key)*multiple
             for key in self.convergence.keys()]
        )

        if self.thresh == "baker":
            energy_converged = False
            if self.cur_cycle > 0:
                cur_energy = self.energies[-1]
                prev_energy = self.energies[-2]
                energy_converged = abs(cur_energy - prev_energy) < 1e-6
            converged = (max_force < 3e-4) and (energy_converged or (max_step < 3e-4))
            return converged
        return any((normal_convergence, overachieved))

    def print_header(self):
        hs = "max(force) rms(force) max(step) rms(step) s/cycle".split()
        header = "cycle" + " ".join([h.rjust(13) for h in hs])
        print(header)

    def print_opt_progress(self):
        int_fmt = "{:>5d}"
        float_fmt = "{:>12.6f}"
        conv_str = int_fmt + " " + (float_fmt + " ") * 4 + "{:>12.1f}"
        print(conv_str.format(
            self.cur_cycle, self.max_forces[-1], self.rms_forces[-1],
            self.max_steps[-1], self.rms_steps[-1], self.cycle_times[-1])
        )
        try:
            # Geometries/ChainOfStates objects can also do some printing.
            add_info = self.geometry.get_additional_print()
            print(add_info)
        except AttributeError:
            pass

    def scale_by_max_step(self, steps):
        steps_max = np.abs(steps).max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def prepare_opt(self):
        pass

    @abc.abstractmethod
    def optimize(self):
        raise Exception("Not implemented!")

    def write_to_out_dir(self, out_fn, content, mode="w"):
        out_path = self.out_dir / out_fn
        with open(out_path, mode) as handle:
            handle.write(content)

    def write_image_trjs(self):
        base_name = "image_{:03d}.trj"
        for i, image in enumerate(self.geometry.images):
            image_fn = base_name.format(i)
            comment = f"cycle {self.cur_cycle}"
            as_xyz = image.as_xyz(comment)
            self.write_to_out_dir(image_fn, as_xyz+"\n", "a")

    def write_results(self):
        # Save results from the Geometry.
        results = self.geometry.results
        # Results will be a list for COS geometries, instead of a
        # dictionary.
        if not self.is_cos:
            results["cart_coords"] = self.cart_coords[-1]
            results["atoms"] = self.geometry.atoms
        self.image_results.append(results)
        self.write_to_out_dir(self.image_results_fn,
                              yaml.dump(self.image_results))

        # Save results from the Optimizer
        opt_results = {la: getattr(self, la) for la in self.list_attrs}
        self.write_to_out_dir(self.opt_results_fn,
                              yaml.dump(opt_results))

    def write_cycle_to_file(self):
        as_xyz_str = self.geometry.as_xyz()

        if self.is_cos:
            out_fn = "cycle_{:03d}.trj".format(self.cur_cycle)
            self.write_to_out_dir(out_fn, as_xyz_str)
            # Also write separate .trj files for every image in the cos
            self.write_image_trjs()
        else:
            # Append to .trj file
            self.out_trj_handle.write(as_xyz_str+"\n")
            self.out_trj_handle.flush()
        self.write_results()

    def final_summary(self):
        # If the optimization was stopped _forces may not be set, so
        # then we force a calculation if it was not already set.
        _ = self.geometry.forces
        cart_forces = self.geometry._forces
        max_cart_forces = np.abs(cart_forces).max()
        rms_cart_forces = np.sqrt(np.mean(cart_forces**2))
        int_str = ""
        if self.geometry.coord_type != "cart":
            int_forces = self.geometry.forces
            max_int_forces = np.abs(int_forces).max()
            rms_int_forces = np.sqrt(np.mean(int_forces**2))
            int_str = f"""
            \tmax(forces,internal): {max_int_forces:.6f} hartree/(bohr,rad)
            \trms(forces,internal): {rms_int_forces:.6f} hartree/(bohr,rad)"""
        energy = self.geometry.energy
        final_summary = f"""
        Final summary:{int_str}
        \tmax(forces,cartesian): {max_cart_forces:.6f} hartree/bohr
        \trms(forces,cartesian): {rms_cart_forces:.6f} hartree/bohr
        \tenergy: {energy:.8f} hartree
        """
        return textwrap.dedent(final_summary.strip())

    def run(self):
        if not self.restarted:
            prep_start_time = time.time()
            self.prepare_opt()
            prep_end_time = time.time()
            prep_time = prep_end_time - prep_start_time
            print(f"Spent {prep_time:.1f} s preparing the first cycle.")

        self.print_header()
        self.stopped = False
        # Actual optimization loop
        for self.cur_cycle in range(self.last_cycle, self.max_cycles):
            start_time = time.time()
            self.log(highlight_text(f"Cycle {self.cur_cycle:03d}"))

            if self.is_cos and self.check_coord_diffs:
                image_coords = [image.cart_coords for image in self.geometry.images]
                align = len(image_coords[0]) > 3
                cds = get_coords_diffs(image_coords, align=align)
                cds_str = np.array2string(cds, precision=4)
                self.log(f"Coordinate differences: {cds_str}")
                # Differences of coordinate differences ;)
                cds_diffs = np.diff(cds)
                min_ind = cds_diffs.argmin()
                if cds_diffs[min_ind] < self.coord_diff_thresh:
                    similar_inds = min_ind, min_ind+1
                    msg = f"Cartesian coordinates of images {similar_inds} are " \
                           "too similar. Stopping optimization!"
                    # I should improve my logging :)
                    print(msg)
                    self.log(msg)
                    break

            # Check if something considerably changed in the optimization,
            # e.g. new images were added/interpolated. Then the optimizer
            # should be reset.
            reset_flag = False
            if self.cur_cycle > 0 and self.is_cos:
                reset_flag = self.geometry.prepare_opt_cycle(self.coords[-1],
                                                             self.energies[-1],
                                                             self.forces[-1])
            self.coords.append(self.geometry.coords.copy())
            self.cart_coords.append(self.geometry.cart_coords.copy())
            if reset_flag:
                self.reset()

            step = self.optimize()

            if step is None:
                # Remove the previously added coords
                self.coords.pop(-1)
                self.cart_coords.pop(-1)
                continue

            if self.is_cos:
                self.tangents.append(self.geometry.get_tangents())

            self.steps.append(step)

            # Convergence check
            self.is_converged = self.check_convergence()

            end_time = time.time()
            elapsed_seconds = end_time - start_time
            self.cycle_times.append(elapsed_seconds)

            if self.dump:
                self.write_cycle_to_file()
                with open(self.current_fn, "w") as handle:
                    handle.write(self.geometry.as_xyz())

            if self.dump and self.dump_restart \
               and (self.cur_cycle % self.dump_restart) == 0:
                self.dump_restart_info()

            self.print_opt_progress()
            if self.is_converged:
                print("Converged!")
                print()
                break

            # Update coordinates
            new_coords = self.geometry.coords.copy() + step
            self.geometry.coords = new_coords

            if hasattr(self.geometry, "reparametrize"):
                reparametrized = self.geometry.reparametrize()
                cur_coords = self.geometry.coords
                prev_coords = self.coords[-1]

                if reparametrized and (cur_coords.size == prev_coords.size):
                    self.log("Did reparametrization")

                    rms = np.sqrt(np.mean((prev_coords - cur_coords)**2))
                    self.log(f"rms of coordinates after reparametrization={rms:.6f}")
                    self.is_converged = rms < self.reparam_thresh
                    if self.is_converged:
                        print("Insignificant coordinate change after "
                              "reparametrization. Signalling convergence!"
                        )
                        print()
                        break

            sys.stdout.flush()
            if check_for_stop_sign():
                self.stopped = True
                break

            self.log("")
        else:
            print("Number of cycles exceeded!")

        # Outside loop
        if self.dump:
            self.out_trj_handle.close()

        if (not self.is_cos) and (not self.stopped):
            print(self.final_summary())
            # Remove 'current_geometry.xyz' file
            try:
                os.remove(self.current_fn)
            except FileNotFoundError:
                self.log(f"Tried to delete '{self.current_fn}'. Couldn't find it.")
        with open(self.final_fn, "w") as handle:
            handle.write(self.geometry.as_xyz())
        print(f"Wrote final, hopefully optimized, geometry to '{self.final_fn.name}'")
        sys.stdout.flush()

    def _get_opt_restart_info(self):
        """To be re-implemented in the derived classes."""
        return dict()

    def _set_opt_restart_info(self, opt_restart_info):
        """To be re-implemented in the derived classes."""
        return

    def get_restart_info(self):
        restart_info = {
            "geom_info": self.geometry.get_restart_info(),
            "last_cycle": self.cur_cycle,
            "max_cycles": self.max_cycles,
            "energies": self.energies,
            "coords": self.coords,
            "forces": [forces.tolist() for forces in self.forces],
            "steps": [step.tolist() for step in self.steps],
        }
        restart_info.update(self._get_opt_restart_info())
        return restart_info

    def set_restart_info(self, restart_info):
        # Set restart information general to all optimizers
        self.last_cycle = restart_info["last_cycle"] + 1

        if self.last_cycle >= self.max_cycles:
            self.max_cycles += restart_info["max_cycles"]

        self.coords = [np.array(coords) for coords in restart_info["coords"]]
        self.energies = restart_info["energies"]
        self.forces = [np.array(forces) for forces in restart_info["forces"]]
        self.steps = [np.array(step) for step in restart_info["steps"]]

        # Set subclass specific information
        self._set_opt_restart_info(restart_info)

        # Propagate restart information downwards to the geometry
        self.geometry.set_restart_info(restart_info["geom_info"])

    def dump_restart_info(self):
        restart_info = self.get_restart_info()

        restart_fn = f"restart_{self.cur_cycle:03d}.yaml"
        restart_yaml = yaml.dump(restart_info)
        self.write_to_out_dir(restart_fn, restart_yaml)

#!/usr/bin/env python3

# https://verahill.blogspot.de/2013/06/439-calculate-frequencies-from-hessian.html
# https://chemistry.stackexchange.com/questions/74639

import logging
import pathlib
import sys

import h5py
import numpy as np

from pysisyphus.xyzloader import make_trj_str, make_xyz_str
from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers import check_for_stop_sign, highlight_text, eigval_to_wavenumber, rms
from pysisyphus.TablePrinter import TablePrinter


class IRC:

    def __init__(self, geometry, step_length=0.1, max_cycles=150,
                 downhill=False, forward=True, backward=True,
                 mode=0, hessian_init=None,
                 displ="energy", displ_energy=5e-4, displ_length=0.1,
                 rms_grad_thresh=5e-4, dump_fn="irc_data.h5", dump_every=5):
        assert(step_length > 0), "step_length must be positive"
        assert(max_cycles > 0), "max_cycles must be positive"

        self.logger = logging.getLogger("irc")

        self.geometry = geometry
        assert self.geometry.coord_type == "cart"

        self.step_length = step_length
        self.max_cycles = max_cycles
        self.downhill = downhill
        # Disable forward/backward when downhill is set
        self.forward = not self.downhill and forward
        self.backward = not self.downhill and backward
        self.mode = mode
        # Load initial (not massweighted) cartesian hessian if provided
        self.hessian_init = hessian_init
        if self.hessian_init is not None:
            self.hessian_init = np.loadtxt(hessian_init)
        self.displ = displ
        assert self.displ in ("energy", "length"), \
            "displ must be either 'energy' or 'length'"
        self.displ_energy = float(displ_energy)
        self.displ_length = float(displ_length)
        self.rms_grad_thresh = float(rms_grad_thresh)
        self.dump_fn = dump_fn
        self.dump_every = int(dump_every)

        self.all_energies = list()
        self.all_coords = list()
        self.all_gradients = list()
        self.all_mw_coords = list()
        self.all_mw_gradients = list()

        # step length dE max(|grad|) rms(grad)
        col_fmts = "int float float float float".split()
        header = ("Step", "IRC length", "dE / au", "max(|grad|)", "rms(grad)")
        self.table = TablePrinter(header, col_fmts)

        self.cur_cycle = 0
        self.converged = False

    @property
    def coords(self):
        return self.geometry.coords

    @coords.setter
    def coords(self, coords):
        self.geometry.coords = coords

    @property
    def mw_coords(self):
        return self.geometry.mw_coords

    @mw_coords.setter
    def mw_coords(self, mw_coords):
        self.geometry.mw_coords = mw_coords

    @property
    def energy(self):
        return self.geometry.energy

    @property
    def gradient(self):
        return self.geometry.gradient

    @property
    def mw_gradient(self):
        return self.geometry.mw_gradient

    @property
    def mw_hessian(self):
        # TODO: This can be removed when the mw_hessian property is updated
        #       in Geometry.py.
        self.geometry.hessian
        return self.geometry.mw_hessian

    def log(self, msg):
        # self.logger.debug(f"step {self.cur_cycle:03d}, {msg}")
        self.logger.debug(msg)

    # def un_massweight(self, vec):
        # return vec * np.sqrt(self.geometry.masses_rep)

    def prepare(self, direction):
        self.cur_cycle = 0
        self.converged = False

        self.irc_energies = list()
        # Not mass-weighted
        self.irc_coords = list()
        self.irc_gradients = list()
        # Mass-weighted
        self.irc_mw_coords = list()
        self.irc_mw_gradients = list()

        # We don't need an initiald displacement when going downhill
        if self.downhill:
            return

        # Over the course of the IRC the hessian may get updated.
        # Copying the TS hessian here ensures a clean start in combined
        # forward and backward runs. Otherwise we would accidently use
        # the updated hessian from the end of the first run for the second
        # run.
        self.hessian = self.ts_hessian.copy()

        # Do inital displacement from the TS
        init_factor = 1 if (direction == "forward") else -1
        initial_step = init_factor*self.init_displ
        self.coords = self.ts_coords + initial_step
        initial_step_length = np.linalg.norm(initial_step)
        self.logger.info(f"Did inital step of length {initial_step_length:.4f} "
                          "from the TS.")

    def initial_displacement(self):
        """Returns a non-mass-weighted step in angstrom for an initial
        displacement from the TS along the transition vector.

        See https://aip.scitation.org/doi/pdf/10.1063/1.454172?class=pdf
        """
        mm_sqr_inv = self.geometry.mm_sqrt_inv
        mw_hessian = self.geometry.mw_hessian
        try:
            if not self.geometry.calculator.analytical_2d:
                mw_hessian = self.geometry.eckart_projection(mw_hessian)
        except AttributeError:
            pass
        eigvals, eigvecs = np.linalg.eigh(mw_hessian)
        neg_inds = eigvals < -1e-8
        assert sum(neg_inds) > 0, "The hessian does not have any negative eigenvalues!"
        min_eigval = eigvals[self.mode]
        self.log(f"Transition vector is mode {self.mode} with wavenumber "
                 f"{eigval_to_wavenumber(min_eigval):.2f} cm⁻¹.")
        mw_trans_vec = eigvecs[:,self.mode]
        self.mw_transition_vector = mw_trans_vec
        # Un-mass-weight the transition vector
        trans_vec = mm_sqr_inv.dot(mw_trans_vec)
        self.transition_vector = trans_vec / np.linalg.norm(trans_vec)

        if self.downhill:
            step = np.zeros_like(self.transition_vector)
        elif self.displ == "length":
            self.log("Using length-based initial displacement from the TS.")
            step = self.displ_length * self.transition_vector
        else:
            # Calculate the length of the initial step away from the TS to initiate
            # the IRC/MEP. We assume a quadratic potential and calculate the
            # displacement for a given energy lowering.
            # dE = (k*dq**2)/2 (dE = energy lowering, k = eigenvalue corresponding
            # to the transition vector/imaginary mode, dq = step length)
            # dq = sqrt(dE*2/k)
            # See 10.1021/ja00295a002 and 10.1063/1.462674
            # 10.1002/jcc.540080808 proposes 3 kcal/mol as initial energy lowering
            self.log("Using energy-based initial displacement from the TS.")
            step_length = np.sqrt(self.displ_energy*2
                                  / np.abs(min_eigval)
            )
            # This calculation is derived from the mass-weighted hessian, so we
            # probably have to multiply this step length with the mass-weighted
            # mode and un-weigh it.
            mw_step = step_length * mw_trans_vec
            step = mw_step / np.sqrt(self.geometry.masses_rep)
        print(f"Norm of initial displacement step: {np.linalg.norm(step):.4f}")
        return step

    def irc(self, direction):
        self.log(highlight_text(f"IRC {direction}"))
        self.prepare(direction)
        # Calculate gradient
        self.gradient
        self.irc_energies.append(self.energy)
        # Non mass-weighted
        self.irc_coords.append(self.coords)
        self.irc_gradients.append(self.gradient)
        # Mass-weighted
        self.irc_mw_coords.append(self.mw_coords)
        self.irc_mw_gradients.append(self.mw_gradient)

        self.table.print_header()
        while True:
            self.log(highlight_text(f"IRC step {self.cur_cycle:03d}") + "\n")
            if self.cur_cycle == self.max_cycles:
                print("IRC steps exceeded. Stopping.")
                print()
                break

            # Do macroiteration/IRC step to update the geometry
            self.step()

            # Calculate energy and gradient on the new geometry
            self.irc_energies.append(self.energy)
            # Non mass-weighted
            self.irc_coords.append(self.coords)
            self.irc_gradients.append(self.gradient)
            # Mass-weighted
            self.irc_mw_coords.append(self.mw_coords)
            self.irc_mw_gradients.append(self.mw_gradient)

            rms_grad = rms(self.gradient)

            irc_length = np.linalg.norm(self.irc_mw_coords[0] - self.irc_mw_coords[-1])
            dE = self.irc_energies[-1] - self.irc_energies[-2]
            max_grad = np.abs(self.gradient).max()

            row_args = (self.cur_cycle, irc_length, dE, max_grad, rms_grad)
            self.table.print_row(row_args)
            try:
                # The derived IRC classes may want to do some printing
                add_info = self.get_additional_print()
                self.table.print(add_info)
            except AttributeError:
                pass
            last_energy = self.irc_energies[-2]
            this_energy = self.irc_energies[-1]

            break_msg = ""
            if self.converged:
                break_msg = "Integrator indicated convergence!"
            elif rms_grad <= self.rms_grad_thresh:
                break_msg = "RMS of gradient converged!"
                self.converged = True
            # TODO: Allow some threshold?
            elif this_energy > last_energy:
                break_msg = "Energy increased!"
            elif abs(last_energy - this_energy) <= 1e-6:
                break_msg = "Energy converged!"
                self.converged = True

            dumped = (self.cur_cycle % self.dump_every) == 0
            if dumped:
                dump_fn = f"{direction}_{self.dump_fn}"
                self.dump_data(dump_fn)

            if break_msg:
                self.table.print(break_msg)
                break

            self.cur_cycle += 1
            if check_for_stop_sign():
                break
            self.log("")
            sys.stdout.flush()

        if direction == "forward":
            self.irc_energies.reverse()
            self.irc_coords.reverse()
            self.irc_gradients.reverse()
            self.irc_mw_coords.reverse()
            self.irc_mw_gradients.reverse()

        if not dumped:
            self.dump_data

    def set_data(self, prefix):
        energies_name = f"{prefix}_energies"
        coords_name = f"{prefix}_coords"
        grad_name = f"{prefix}_gradients"
        mw_coords_name = f"{prefix}_mw_coords"
        mw_grad_name = f"{prefix}_mw_gradients"

        setattr(self, coords_name, self.irc_coords)
        setattr(self, grad_name, self.irc_gradients)
        setattr(self, mw_coords_name, self.irc_mw_coords)
        setattr(self, mw_grad_name, self.irc_mw_gradients)
        setattr(self, energies_name, self.irc_energies)

        self.all_energies.extend(getattr(self, energies_name))
        self.all_coords.extend(getattr(self, coords_name))
        self.all_gradients.extend(getattr(self, grad_name))
        self.all_mw_coords.extend(getattr(self, mw_coords_name))
        self.all_mw_gradients.extend(getattr(self, mw_grad_name))

        setattr(self, f"{prefix}_step", self.cur_cycle)
        self.write_trj(".", prefix, getattr(self, mw_coords_name))

    def run(self):
        # Calculate data at TS and create backup
        self.ts_coords = self.coords.copy()
        self.ts_mw_coords = self.mw_coords.copy()
        self.ts_gradient = self.gradient.copy()
        self.ts_mw_gradient = self.mw_gradient.copy()
        self.ts_energy = self.energy

        ts_grad_norm = np.linalg.norm(self.ts_gradient)
        ts_grad_max = np.abs(self.ts_gradient).max()
        ts_grad_rms = rms(self.ts_gradient)

        self.log( "Transition state (TS):\n"
                 f"\tnorm(grad)={ts_grad_norm:.8f}\n"
                 f"\t max(grad)={ts_grad_max:.8f}\n"
                 f"\t rms(grad)={ts_grad_rms:.8f}"
        )

        print("IRC length in mw. coords, max(|grad|) and rms(grad) in non-"
              "mass-weighted coords.")

        # For forward/backward runs we need an intial displacement
        # and for this we need a hessian, that we calculate now.
        # For downhill runs we probably dont need a hessian.
        if not self.downhill:
            if self.hessian_init is not None:
                self.ts_hessian = self.hessian_init.copy()
            else:
                self.ts_hessian = self.geometry.hessian.copy()
            self.init_displ = self.initial_displacement()

        if self.forward:
            print("\n" + highlight_text("IRC - Forward") + "\n")
            self.irc("forward")
            self.set_data("forward")

        # Add TS/starting data
        self.all_energies.append(self.ts_energy)
        self.all_coords.append(self.ts_coords)
        self.all_gradients.append(self.ts_gradient)
        self.all_mw_coords.append(self.ts_mw_coords)
        self.all_mw_gradients.append(self.ts_mw_gradient)
        self.ts_index = len(self.all_energies) - 1

        if self.backward:
            print("\n" + highlight_text("IRC - Backward") + "\n")
            self.irc("backward")
            self.set_data("backward")

        if self.downhill:
            print("\n" + highlight_text("IRC - Downhill") + "\n")
            self.irc("downhill")
            self.set_data("downhill")

        self.all_mw_coords = np.array(self.all_mw_coords)
        self.all_energies = np.array(self.all_energies)
        self.postprocess()
        if not self.downhill:
            self.write_trj(".", "finished")

            # Dump the whole IRC to HDF5
            dump_fn = "finished_" + self.dump_fn
            self.dump_data(dump_fn, full=True)

        # Convert to arrays
        [setattr(self, name, np.array(getattr(self, name)))
         for name in "all_energies all_coords all_gradients "
                     "all_mw_coords all_mw_gradients".split()
        ]

        # Right now self.all_mw_coords is still in mass-weighted coordinates.
        # Convert them to un-mass-weighted coordinates.
        self.all_mw_coords_umw = self.all_mw_coords / self.geometry.masses_rep**0.5

    def postprocess(self):
        pass

    def write_trj(self, path, prefix, coords=None):
        path = pathlib.Path(path)
        atoms = self.geometry.atoms
        if coords is None:
            coords = self.all_mw_coords
        coords = coords.copy()
        coords /= self.geometry.masses_rep**0.5
        coords = coords.reshape(-1, len(atoms), 3) * BOHR2ANG
        # all_mw_coords = self.all_mw_coords.flatten()
        trj_string = make_trj_str(atoms, coords, comments=self.all_energies)
        trj_fn = f"{prefix}_irc.trj"
        with open(path / trj_fn, "w") as handle:
            handle.write(trj_string)

        first_coords = coords[0]
        first_fn = f"{prefix}_first.xyz"
        with open(path / first_fn, "w") as handle:
            handle.write(make_xyz_str(atoms, first_coords))

        last_coords = coords[-1]
        first_fn = f"{prefix}_last.xyz"
        with open(path / first_fn, "w") as handle:
            handle.write(make_xyz_str(atoms, last_coords))

    def get_irc_data(self):
        data_dict = {
            "energies": np.array(self.irc_energies, dtype=float),
            "coords": np.array(self.irc_coords, dtype=float),
            "gradients": np.array(self.irc_gradients, dtype=float),
            "mw_coords": np.array(self.irc_mw_coords, dtype=float),
            "mw_gradients": np.array(self.irc_mw_gradients, dtype=float),
        }
        return data_dict

    def get_full_irc_data(self):
        data_dict = {
            "energies": np.array(self.all_energies, dtype=float),
            "coords": np.array(self.all_coords, dtype=float),
            "gradients": np.array(self.all_gradients, dtype=float),
            "mw_coords": np.array(self.all_mw_coords, dtype=float),
            "mw_gradients": np.array(self.all_mw_gradients, dtype=float),
            "ts_index": np.array(self.ts_index, dtype=int),
        }
        return data_dict

    def dump_data(self, dump_fn=None, full=False):
        get_data = self.get_full_irc_data if full else self.get_irc_data
        data_dict = get_data()

        data_dict.update({
                "atoms": np.array(self.geometry.atoms, dtype="S"),
                "rms_grad_thresh": np.array(self.rms_grad_thresh),
        })

        if dump_fn is None:
            dump_fn = self.dump_fn

        with h5py.File(dump_fn, "w") as handle:
            for key, val in data_dict.items():
                handle.create_dataset(name=key, dtype=val.dtype, data=val)

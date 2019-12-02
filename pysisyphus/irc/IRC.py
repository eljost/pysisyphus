#!/usr/bin/env python3

# https://verahill.blogspot.de/2013/06/439-calculate-frequencies-from-hessian.html
# https://chemistry.stackexchange.com/questions/74639

import copy
import logging
import pathlib
import sys

import numpy as np

from pysisyphus.xyzloader import make_trj_str, make_xyz_str
from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers import check_for_stop_sign, highlight_text
from pysisyphus.TablePrinter import TablePrinter


class IRC:

    def __init__(self, geometry, step_length=0.1, max_cycles=150,
                 downhill=False, forward=True, backward=True, mode=0,
                 displ="energy", displ_energy=5e-4, displ_length=0.1,
                 rms_grad_thresh=1e-4):
        assert(step_length > 0), "step_length must be positive"
        assert(max_cycles > 0), "max_cycles must be positive"

        self.logger = logging.getLogger("irc")

        self.geometry = geometry
        self.step_length = step_length
        self.max_cycles = max_cycles
        self.downhill = downhill
        # Disable forward/backward when downhill is set
        self.forward = not self.downhill and forward
        self.backward = not self.downhill and backward
        self.mode = mode
        self.displ = displ
        assert self.displ in ("energy", "length"), \
            "displ must be either 'energy' or 'length'"
        self.displ_energy = float(displ_energy)
        # assert self.displ_energy > 0, \
            # "displ_energy must be positive"
        self.displ_length = float(displ_length)
        # assert self.displ_length > 0, \
            # "displ_displ must be positive"
        self.rms_grad_thresh = float(rms_grad_thresh)

        self.all_coords = list()
        self.all_energies = list()

        # Backup TS data
        self.ts_coords = self.coords.copy()
        self.ts_mw_coords = self.mw_coords.copy()
        # self.ts_gradient = self.geometry.gradient.copy()
        self.ts_mw_gradient = self.mw_gradient.copy()
        self.ts_energy = self.energy.copy()
        self.ts_hessian = self.geometry.hessian.copy()

        self.cur_step = 0
        self.converged = False
        # With downhill=True we shouldn't need any initial displacement.
        # We still call the method because here the initial hessian is
        # calculated and some sanity checks are made. The returned init_displ
        # will be the zero vector though.
        self.init_displ = self.initial_displacement()
        # step length dE max(|grad|) rms(grad)
        col_fmts = "int float float float float".split()
        header = ("Step", "IRC length", "dE", "max(|grad|)", "rms(grad)")
        self.table = TablePrinter(header, col_fmts)

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
        self.logger.debug(f"step {self.cur_step:03d}, {msg}")

    # def un_massweight(self, vec):
        # return vec * np.sqrt(self.geometry.masses_rep)

    def prepare(self, direction):
        self.cur_step = 0
        self.converged = False
        # Over the course of the IRC the hessian may get updated.
        # Copying the TS hessian here ensures a clean start in combined
        # forward and backward runs. Otherwise we would accidently use
        # the updated hessian from the end of the first run for the second
        # run.
        self.hessian = self.ts_hessian

        self.irc_mw_coords = list()
        self.irc_energies = list()
        self.irc_gradients = list()
        self.irc_mw_gradients = list()

        # We don't need an initiald displacement when going downhill
        if self.downhill:
            return

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
        neg_inds = eigvals < -1e-10
        assert sum(neg_inds) > 0, "The hessian does not have any negative eigenvalues!"
        min_eigval = eigvals[self.mode]
        mw_trans_vec = eigvecs[:,self.mode]
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
        self.logger.info(f"IRC {direction}")
        self.prepare(direction)
        gradient = self.gradient
        self.irc_mw_coords.append(self.mw_coords)
        self.irc_energies.append(self.energy)
        self.irc_gradients.append(self.gradient)
        self.irc_mw_gradients.append(self.mw_gradient)

        self.table.print_header()
        while True:
            if self.cur_step == self.max_cycles:
                print("IRC steps exceeded. Stopping.")
                print()
                break

            # Do macroiteration/IRC step to update the geometry
            self.step()
            # Calculate energy and gradient on the new geometry
            self.irc_mw_coords.append(self.mw_coords)
            self.irc_gradients.append(self.gradient)
            self.irc_mw_gradients.append(self.mw_gradient)
            self.irc_energies.append(self.energy)
            rms_grad = np.sqrt(np.mean(np.square(self.gradient)))

            irc_length = np.linalg.norm(self.irc_mw_coords[0] - self.irc_mw_coords[-1])
            dE = self.irc_energies[-1] - self.irc_energies[-2]
            max_grad = np.abs(self.gradient).max()

            if self.converged:
                break_msg = "Optimizer indicated convergence!"
                self.table.print(break_msg)
                break

            row_args = (self.cur_step, irc_length, dE, max_grad, rms_grad)
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
            if rms_grad <= self.rms_grad_thresh:
                break_msg = "RMS of gradient converged!"
            elif this_energy > last_energy:
                break_msg = "Energy increased!"
            elif abs(last_energy - this_energy) <= 1e-6:
                break_msg = "Energy converged!"

            if break_msg:
                self.table.print(break_msg)
                break

            self.cur_step += 1
            if check_for_stop_sign():
                break
            self.log("")
            sys.stdout.flush()

        if direction == "forward":
            self.irc_mw_coords.reverse()
            self.irc_energies.reverse()

    def run(self):
        if self.forward:
            print(highlight_text("Forward"))
            # try:
                # self.irc("forward")
            # except Exception as error:
                # logging.error(error)
            self.irc("forward")
            self.forward_coords = self.irc_mw_coords
            self.forward_energies = self.irc_energies
            self.all_coords.extend(self.forward_coords)
            self.all_energies.extend(self.forward_energies)
            self.forward_step = self.cur_step
            self.write_trj(".", "forward", self.forward_coords)

        # Add TS/starting data
        self.all_coords.append(self.ts_mw_coords)
        self.all_energies.append(self.ts_energy)

        if self.backward:
            print(highlight_text("Backward"))
            # try:
                # self.irc("backward")
            # except Exception as error:
                # logging.error(error)
            self.irc("backward")
            self.backward_coords = self.irc_mw_coords
            self.backward_energies = self.irc_energies
            self.all_coords.extend(self.backward_coords)
            self.all_energies.extend(self.backward_energies)
            self.backward_step = self.cur_step
            self.write_trj(".", "backward", self.backward_coords)

        if self.downhill:
            print(highlight_text("Downhill"))
            self.irc("downhill")
            self.downhill_coords = self.irc_mw_coords
            self.downhill_energies = self.irc_energies
            self.all_coords.extend(self.downhill_coords)
            self.all_energies.extend(self.downhill_energies)
            self.downhill_step = self.cur_step
            self.write_trj(".", "downhill", self.downhill_coords)

        self.all_coords = np.array(self.all_coords)
        self.all_energies = np.array(self.all_energies)
        self.postprocess()
        if not self.downhill:
            self.write_trj(".", "finished")

        # Right now self.all_coords is still in mass-weighted coordinates.
        # Convert them to un-mass-weighted coordinates.
        self.all_coords_umw = self.all_coords / self.geometry.masses_rep**0.5

    def postprocess(self):
        pass

    def write_trj(self, path, prefix, coords=None):
        path = pathlib.Path(path)
        atoms = self.geometry.atoms
        if coords is None:
            coords = self.all_coords
        coords = coords.copy()
        coords /= self.geometry.masses_rep**0.5
        coords = coords.reshape(-1, len(atoms), 3) * BOHR2ANG
        # all_coords = self.all_coords.flatten()
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

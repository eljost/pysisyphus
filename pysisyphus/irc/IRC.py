# https://verahill.blogspot.de/2013/06/439-calculate-frequencies-from-hessian.html
# https://chemistry.stackexchange.com/questions/74639

import logging
from math import ceil, log
import os
from pathlib import Path
import sys

import h5py
import numpy as np

from pysisyphus.constants import BOHR2ANG, AU2KJPERMOL
from pysisyphus.helpers import check_for_end_sign, rms
from pysisyphus.helpers_pure import (
    highlight_text,
    eigval_to_wavenumber,
    report_isotopes,
)
from pysisyphus.irc.initial_displ import cubic_displ, third_deriv_fd, cubic_displ_for_h5
from pysisyphus.irc.Instanton import T_crossover_from_eigval
from pysisyphus.io import save_third_deriv
from pysisyphus.optimizers.guess_hessians import get_guess_hessian
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.xyzloader import make_trj_str, make_xyz_str


class IRC:
    valid_displs = ("energy", "length", "energy_cubic")

    def __init__(
        self,
        geometry,
        step_length=0.1,
        max_cycles=125,
        downhill=False,
        forward=True,
        backward=True,
        root=0,
        hessian_init=None,
        displ="energy",
        displ_energy=1e-3,
        displ_length=0.1,
        displ_third_h5=None,
        rms_grad_thresh=1e-3,
        hard_rms_grad_thresh=None,
        energy_thresh=1e-6,
        imag_below=0.0,
        force_inflection=True,
        check_bonds=True,
        out_dir=".",
        prefix="",
        dump_fn="irc_data.h5",
        dump_every=5,
    ):
        """Base class for IRC calculations.

        Parameters
        ----------
        geometry : Geometry
            Transtion state geometry, or initial geometry for downhill run.
        step_length : float, optional
            Step length in unweighted coordinates.
        max_cycles : int, optional
            Positive integer, controlloing the maximum number of IRC steps
            taken in a direction (forward/backward/downhill).
        downhill : bool, default=False
            Downhill run from a non-stationary point with non-vanishing
            gradient. Disables forward and backward runs.
        forward : bool, default=True
            Integrate IRC in positive s direction.
        backward : bool, default=True
            Integrate IRC in negative s direction.
        root : int, default=0
            Use n-th root for initial displacement from TS.
        hessian_init : str, default=None
            Path to Hessian HDF5 file, e.g., from a previous TS calculation.
        displ: str, one of ("energy", "length", "energy_cubic")
            Controlls initial displacement from the TS. 'energy' assumes a
            quadratic model, from which a step length for a given energy
            lowering (see 'displ_energy') is determined. 'length' corresponds
            to a displacement along the transition vector. 'energy_cubic' considers
            3rd derivatives of the energy along the transition vector.
        displ_energy : float, default=1e-3
            Required energy lowering from the TS in au (Hartree). Used with
            'displ: energy|energy_cubic'.
        displ_length : float, default=0.1
            Step length along the transition vector. Used only with
            'displ: length'.
        displ_third_h5: str, optional
            Path to HDF5 file containing 3rd derivative information. Used with
            'displ: energy_cubic'.
        rms_grad_thresh : float, default=1e-3,
            Convergence is signalled when to root mean square of the unweighted
            gradient is less than or equal to this value
        energy_thresh : float, default=1e-6,
            Signal convergence when the energy difference between two points
            is equal to or less than 'energy_thresh'.
        imag_below : float, default=0.0
            Require the wavenumber of the imaginary mode to be below the
            given threshold.
        force_inflection : bool, optional
            Don't indicate convergence before passing an inflection point.
        check_bonds : bool, optional, default=True
            Report whether bonds are formed/broken along the IRC, w.r.t the TS.
        out_dir : str, optional
            Dump everything into 'out_dir' directory instead of the CWD.
        prefix : str, optional
            Short string that is prepended to all files that are created
            by this class, e.g., trajectories and HDF5 dumps.
        dump_fn : str, optional
            Base name for the HDF5 files.
        dump_every : int, optional
            Dump to HDF5 every n-th cycle.
        """
        assert step_length > 0, "step_length must be positive"
        assert max_cycles > 0, "max_cycles must be positive"

        self.logger = logging.getLogger("irc")

        self.geometry = geometry
        self.atoms = self.geometry.atoms
        assert self.geometry.coord_type == "cart"

        report_isotopes(self.geometry, "the IRC")

        self.step_length = step_length
        self.max_cycles = max_cycles
        self.downhill = downhill
        # Disable forward/backward when downhill is set
        self.forward = not self.downhill and forward
        self.backward = not self.downhill and backward
        self.root = root
        if hessian_init is None:
            hessian_init = "calc" if not self.downhill else "unit"
        self.hessian_init = hessian_init
        self.displ = displ
        assert (
            self.displ in self.valid_displs
        ), f"'displ: {self.displ}' not in {self.valid_displs}"
        self.displ_energy = float(displ_energy)
        self.displ_length = float(displ_length)
        self.displ_third_h5 = displ_third_h5
        self.rms_grad_thresh = float(rms_grad_thresh)
        self.hard_rms_grad_thresh = hard_rms_grad_thresh
        self.energy_thresh = float(energy_thresh)
        assert imag_below <= 0.0
        self.imag_below = imag_below
        self.force_inflection = force_inflection
        self.check_bonds = check_bonds
        self.out_dir = out_dir
        self.out_dir = Path(self.out_dir)
        if not self.out_dir.exists():
            os.mkdir(self.out_dir)
        self.prefix = f"{prefix}_" if prefix else prefix
        self.dump_fn = dump_fn
        self.dump_every = int(dump_every)

        # Determine bonds at TS
        self.ts_bond_sets = self.geometry.bond_sets
        self.ref_bond_sets = {}

        self._m_sqrt = np.sqrt(self.geometry.masses_rep)
        self.all_energies = list()
        self.all_coords = list()
        self.all_gradients = list()
        self.all_mw_coords = list()
        self.all_mw_gradients = list()

        # step length dE max(|grad|) rms(grad)
        col_fmts = "int float float float float".split()
        header = ("Step", "IRC length", "dE / au", "max(|grad|)", "rms(grad)")
        self.table = TablePrinter(header, col_fmts)

        self.cycle_places = ceil(log(self.max_cycles, 10))

    def get_path_for_fn(self, fn):
        return self.out_dir / f"{self.prefix}{fn}"

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

    # @property
    # def mw_hessian(self):
    # # TODO: This can be removed when the mw_hessian property is updated
    # #       in Geometry.py.
    # return self.geometry.mw_hessian

    def log(self, msg):
        # self.logger.debug(f"step {self.cur_cycle:03d}, {msg}")
        self.logger.debug(msg)

    @property
    def m_sqrt(self):
        return self._m_sqrt

    def unweight_vec(self, vec):
        return self.m_sqrt * vec

    def mass_weigh_hessian(self, hessian):
        return self.geometry.mass_weigh_hessian(hessian)

    def prepare(self, direction):
        self.direction = direction
        self.converged = False
        self.energy_increased = False
        self.energy_converged = False
        self.past_inflection = not self.force_inflection

        self.irc_energies = list()
        # Not mass-weighted
        self.irc_coords = list()
        self.irc_gradients = list()
        # Mass-weighted
        self.irc_mw_coords = list()
        self.irc_mw_gradients = list()

        self.ref_bond_sets = self.ts_bond_sets.copy()

        # Over the course of the IRC the hessian may get updated.
        # Copying the initial hessian here ensures a clean start in combined
        # forward and backward runs. Otherwise we may accidentally use
        # the updated hessian from the end of the first run for the second
        # run.
        self.mw_hessian = self.mass_weigh_hessian(self.init_hessian)

        trj_fn = self.get_path_for_fn(f"{direction}_irc.trj")
        self.trj_handle = open(trj_fn, "w")

        # We don't need an initial displacement when going downhill
        if self.downhill:
            return

        # Do inital displacement from the TS
        if direction == "forward":
            initial_step = self.init_displ_plus
        elif direction == "backward":
            initial_step = self.init_displ_minus
        else:
            raise Exception("Invalid direction='{direction}'!")
        self.coords = self.ts_coords + initial_step

        if self.displ in ("energy", "energy_cubic"):
            actual_energy = self.energy
            actual_lowering = self.ts_energy - actual_energy
            diff = self.displ_energy - actual_lowering

            def en_str(en):
                return f"{en: .4f} au ({en*AU2KJPERMOL: .2f} kJ mol⁻¹)"

            print(
                f"Requested energy lowering: {en_str(self.displ_energy)}\n"
                f"   Actual energy lowering: {en_str(actual_lowering)}\n"
                f"                        Δ: {en_str(diff)}"
            )
            if actual_lowering < 0.0:
                print("Displaced geometry is higher in energy compared to TS!")
            print("\n")
            sys.stdout.flush()
        initial_step_length = np.linalg.norm(initial_step)
        self.logger.info(
            f"Did inital step of length {initial_step_length:.4f} " "from the TS."
        )

    def initial_displacement(self):
        """Returns non-mass-weighted steps in +s and -s direction
        for initial displacement from the TS. Earlier version only
        returned one step, that was later multiplied by either 1 or -1,
        depending on the desired IRC direction (forward/backward).
        The current implementation directly returns two steps for forward
        and backward direction. Whereas for plus and minus steps for
        displ 'length' and displ 'energy'
            step_plus = -step_minus
        is valid, it is not valid for dipsl 'energy_cubic' anymore. The
        latter step is formed as
            x(ds) = ds * v0 + ds**2 * v1
        so
            x(ds) != -x(ds)
        as
            ds * v0 + ds**2 * v1 != -ds * v0 - ds**2 * v1 .

        So, all required step are formed directly and later used as appropriate.

        See
            https://aip.scitation.org/doi/pdf/10.1063/1.454172
            https://pubs.acs.org/doi/10.1021/j100338a027
            https://aip.scitation.org/doi/pdf/10.1063/1.459634
        """

        mw_hessian = self.geometry.mass_weigh_hessian(self.init_hessian)
        if self.coords.size > 3:
            proj_hessian, P = self.geometry.eckart_projection(mw_hessian, return_P=True)
        # Don't project single atom species and analytical potentials
        else:
            proj_hessian = mw_hessian
            P = np.eye(self.coords.size)
        eigvals, eigvecs = np.linalg.eigh(proj_hessian)
        mw_cart_displs = P.T.dot(eigvecs)
        cart_displs = self.geometry.mm_sqrt_inv.dot(mw_cart_displs)
        nus = eigval_to_wavenumber(eigvals)
        nu_root = nus[self.root]
        assert nu_root <= self.imag_below, (
            f"Wavenumber {nu_root:.2f} cm⁻¹ of imaginary mode {self.root} is above "
            f"the threshold of {self.imag_below:.2f} cm⁻¹."
        )
        neg_inds = eigvals < -1e-8
        assert sum(neg_inds) > 0, "The hessian does not have any negative eigenvalues!"

        min_eigval = eigvals[self.root]
        min_nu = nus[self.root]
        T_c = T_crossover_from_eigval(min_eigval)
        min_msg = (
            f"Transition vector is mode {self.root} with wavenumber {min_nu:.2f} cm⁻¹.\n"
            f"Crossover temperature T_c: {T_c:.2f} K"
        )
        # Doing it this way hurts ... I'll have to improve my logging game...
        self.log(min_msg)
        print(min_msg)

        # Mass-weighted
        mw_trans_vec = mw_cart_displs[:, self.root]
        self.mw_transition_vector = mw_trans_vec
        # Not mass-weighted
        trans_vec = cart_displs[:, self.root]
        self.transition_vector = trans_vec / np.linalg.norm(trans_vec)

        if self.downhill:
            mw_step_plus = mw_step_minus = np.zeros_like(self.transition_vector)
            msg = "Downhill run. No initial displacement from the TS."
        elif self.displ == "length":
            msg = "Using length-based initial displacement from the TS."
            mw_step_plus = self.displ_length * mw_trans_vec
            mw_step_minus = -mw_step_plus
        elif self.displ == "energy":
            # Calculate the length of the initial step away from the TS to initiate
            # the IRC/MEP. We assume a quadratic potential and calculate the
            # displacement for a given energy lowering.
            # dE = (k*dq**2)/2 (dE = energy lowering, k = eigenvalue corresponding
            # to the transition vector/imaginary mode, dq = step length)
            # dq = sqrt(dE*2/k)
            # See 10.1021/ja00295a002 and 10.1063/1.462674
            # 10.1002/jcc.540080808 proposes 3 kcal/mol as initial energy lowering
            msg = (
                f"Energy-based (ΔE={self.displ_energy} au) initial displacement from "
                "the TS using 2rd derivatives."
            )
            step_length = np.sqrt(self.displ_energy * 2 / np.abs(min_eigval))
            # This calculation is derived from the mass-weighted hessian, so we
            # have to multiply this step length with the mass-weighted
            # mode and un-weigh it.
            mw_step_plus = step_length * mw_trans_vec
            mw_step_minus = -mw_step_plus
        elif self.displ == "energy_cubic":
            if self.displ_third_h5:
                self.log(f"Loaded 3rd derivative data from '{self.displ_third_h5}'")
                mw_step_plus, mw_step_minus = cubic_displ_for_h5(
                    self.displ_third_h5, -self.displ_energy
                )
            else:
                Gv, third_deriv_res = third_deriv_fd(self.geometry, mw_trans_vec)
                h5_fn = self.get_path_for_fn("third_deriv.h5")
                save_third_deriv(h5_fn, self.geometry, third_deriv_res, H_mw=mw_hessian)
                mw_step_plus, mw_step_minus = cubic_displ(
                    proj_hessian,
                    eigvecs[:, self.root],
                    min_eigval,
                    Gv,
                    -self.displ_energy,
                )
                mw_step_plus = P.T.dot(mw_step_plus)
                mw_step_minus = P.T.dot(mw_step_minus)
            msg = (
                f"Energy-based (ΔE={self.displ_energy} au) initial displacement from "
                "the TS using 3rd derivatives."
            )
        else:
            raise Exception(f"self.displ={self.displ} is invalid!")

        step_plus = mw_step_plus / self.m_sqrt
        step_minus = mw_step_minus / self.m_sqrt
        self.log(msg)
        print(msg)
        print(
            "Initial step lengths (not mass-weighted):\n"
            f"\t Forward: {np.linalg.norm(step_plus):.4f} au\n"
            f"\tBackward: {np.linalg.norm(step_minus):.4f} au"
        )
        return step_plus, step_minus

    def get_conv_fact(self, mw_grad, min_fact=2.0):
        # Numerical integration of differential equations requires a step length and/or
        # we have to terminate the integration at some point, e.g. when the desired
        # step length is reached. IRCs are integrated in mass-weighted coordinates,
        # but self.step_length is given in unweighted coordinates. Unweighting a step
        # in mass-weighted coordinates will reduce its norm as we divide by sqrt(m).
        #
        # If we want to do an Euler-integration we have to decide on a step size
        # when a desired integration length is to be reached in a given number of steps.
        # [3] proposes using Δs/250 with a maximum of 500 steps, so something like
        # Δs/(max_steps / 2). It seems we can't use this because (at
        # least for the systems I tested) this will lead to a step length that is too
        # small, so the predictor Euler-integration will fail to converge in the
        # prescribed number of cycles. It fails because simply dividing the desired
        # step length in unweighted coordinates does not take into account the mass
        # dependence. Such a step size is appropriate for integrations in unweighted
        # coordinates, but not when using mass-weighted coordinates.
        #
        # We determine a conversion factor from comparing the magnitudes (norms) of
        # the mass-weighted and un-mass-weighted gradients. This takes into account
        # which atoms are actually moving, so it should be a good guess.
        norm_mw_grad = np.linalg.norm(mw_grad)
        norm_grad = np.linalg.norm(self.unweight_vec(mw_grad))
        conv_fact = norm_grad / norm_mw_grad
        conv_fact = max(min_fact, conv_fact)
        self.log(f"Un-weighted / mass-weighted conversion factor {conv_fact:.4f}")
        return conv_fact

    def report_bonds(self, prefix, bonds):
        if len(bonds) == 0:
            return

        plural = "s" if len(bonds) > 1 else ""
        bond_strs = list()
        for from_, to_ in bonds:
            from_atom = self.atoms[from_]
            to_atom = self.atoms[to_]
            bond_strs.append(f"[{from_atom}{from_}-{to_atom}{to_}]")
        bonds_str = ", ".join(bond_strs)
        self.table.print(f"Bond{plural} {prefix}: {bonds_str}")

    def irc(self, direction):
        self.log(highlight_text(f"IRC {direction}", level=1))
        self.cur_direction = direction
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
        for self.cur_cycle in range(self.max_cycles):
            self.log(highlight_text(f"IRC step {self.cur_cycle:03d}") + "\n")

            # Dump current coordinates to trj
            comment = f"{direction} IRC, step {self.cur_cycle}"
            coords_str = make_xyz_str(
                self.atoms, BOHR2ANG * self.coords.reshape((-1, 3)), comment
            )
            self.trj_handle.write(coords_str + "\n")
            self.trj_handle.flush()

            self.log(f"Current energy: {self.energy:.6f} au")
            #
            # Take IRC step.
            #
            self.step()

            # Calculate gradient and energy on the new geometry
            # Non mass-weighted
            self.log("Calculating energy and gradient at new geometry.")
            self.irc_coords.append(self.coords)
            self.irc_gradients.append(self.gradient)
            self.irc_energies.append(self.energy)
            # Mass-weighted
            self.irc_mw_coords.append(self.mw_coords)
            self.irc_mw_gradients.append(self.mw_gradient)

            rms_grad = rms(self.gradient)

            # Only update once
            if not self.past_inflection:
                self.past_inflection = rms_grad >= self.rms_grad_thresh
                _ = "" if self.past_inflection else "not yet"
                self.log(f"(rms(grad) > threshold) {_} fullfilled!")

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

            if self.check_bonds:
                cur_bond_sets = self.geometry.bond_sets
                formed = cur_bond_sets - self.ref_bond_sets
                broken = self.ref_bond_sets - cur_bond_sets
                self.report_bonds("formed", formed)
                self.report_bonds("broken", broken)
                # Update bond sets to avoid repeated reporting of bond topology changes
                self.ref_bond_sets -= broken
                self.ref_bond_sets |= formed  # union

            last_energy = self.irc_energies[-2]
            this_energy = self.irc_energies[-1]

            break_msg = ""
            self.energy_increased = this_energy > last_energy
            self.energy_converged = abs(last_energy - this_energy) <= self.energy_thresh
            if self.converged:
                break_msg = "Integrator indicated convergence!"
            elif self.past_inflection and (rms_grad <= self.rms_grad_thresh):
                break_msg = "rms(grad) converged!"
                self.converged = True
            elif (
                self.hard_rms_grad_thresh
                and (not self.past_inflection)
                and (rms_grad <= self.hard_rms_grad_thresh)
            ):
                break_msg = "rms(grad) below hard threshold."
            # TODO: Allow some threshold?
            elif self.energy_increased:
                break_msg = "Energy increased!"
            elif self.energy_converged:
                break_msg = "Energy converged!"
                self.converged = True

            dumped = (self.cur_cycle % self.dump_every) == 0
            if dumped:
                dump_fn = self.get_path_for_fn(f"{direction}_{self.dump_fn}")
                self.dump_data(dump_fn)

            if break_msg:
                self.table.print(break_msg)
                break

            if check_for_end_sign():
                break
            self.log("")
            sys.stdout.flush()
        else:
            print("IRC steps exceeded. Stopping.")
            print()

        if direction == "forward":
            self.irc_energies.reverse()
            self.irc_coords.reverse()
            self.irc_gradients.reverse()
            self.irc_mw_coords.reverse()
            self.irc_mw_gradients.reverse()

        if not dumped:
            self.dump_data(dump_fn)

        self.cur_direction = None
        self.trj_handle.close()

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

        setattr(self, f"{prefix}_is_converged", self.converged)
        setattr(self, f"{prefix}_energy_increased", self.energy_increased)
        setattr(self, f"{prefix}_energy_converged", self.energy_converged)
        setattr(self, f"{prefix}_cycle", self.cur_cycle)
        self.dump_ends(".", prefix, getattr(self, mw_coords_name))

    def report_conv_thresholds(self):
        threshs = [
            f"\t     rms(|gradient|) <= {self.rms_grad_thresh:.6f} E_h a_0⁻¹",
            f"\t             Δenergy <= {self.energy_thresh:.6f} E_h",
        ]
        # Drop hard rms grad item
        if self.hard_rms_grad_thresh is not None:
            threshs.insert(
                1,
                f"\thard rms(|gradient|) <= {self.hard_rms_grad_thresh:.6f} E_h a_0⁻¹",
            )
        print(
            "Convergence thresholds (non mass-weighted gradient):\n"
            + "\n".join(threshs)
            + "\n"
        )

    def run(self):
        self.report_conv_thresholds()
        # Calculate data at TS and create backup
        self.ts_coords = self.coords.copy()
        self.ts_mw_coords = self.mw_coords.copy()
        print("Calculating energy and gradient at TS.")
        self.ts_gradient = self.gradient.copy()
        self.ts_mw_gradient = self.mw_gradient.copy()
        self.ts_energy = self.energy

        ts_grad_norm = np.linalg.norm(self.ts_gradient)
        ts_grad_max = np.abs(self.ts_gradient).max()
        ts_grad_rms = rms(self.ts_gradient)

        self.log(
            "Transition state (TS):\n"
            f"\t    energy={self.ts_energy:.6f} au\n"
            f"\tnorm(grad)={ts_grad_norm:.6f}\n"
            f"\t max(grad)={ts_grad_max:.6f}\n"
            f"\t rms(grad)={ts_grad_rms:.6f}"
        )

        self.init_hessian, hess_str = get_guess_hessian(
            self.geometry,
            self.hessian_init,
            cart_gradient=self.ts_gradient,
            h5_fn=self.get_path_for_fn(f"hess_init_irc.h5"),
        )
        self.log(f"Initial hessian: {hess_str}")

        # For forward/backward runs from a TS we need an intial displacement,
        # calculated from the transition vector (imaginary mode) of the TS
        # hessian. If we need/want a Hessian for a downhill run from a
        # non-stationary point (with non-vanishing gradient) depends on the
        # actual IRC integrator (e.g. EulerPC and LQA need a Hessian).
        if not self.downhill:
            self.init_displ_plus, self.init_displ_minus = self.initial_displacement()

        print(
            "IRC length in mw. coords, max(|grad|) and rms(grad) in "
            "unweighted coordinates."
        )

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
            self.dump_ends(".", "finished", trj=True)

            # Dump the whole IRC to HDF5
            dump_fn = self.get_path_for_fn("finished_" + self.dump_fn)
            self.dump_data(dump_fn, full=True)

        # Convert to arrays
        [
            setattr(self, name, np.array(getattr(self, name)))
            for name in "all_energies all_coords all_gradients "
            "all_mw_coords all_mw_gradients".split()
        ]

        # Right now self.all_mw_coords is still in mass-weighted coordinates.
        # Convert them to un-mass-weighted coordinates.
        self.all_mw_coords_umw = self.all_mw_coords / self.m_sqrt

    def postprocess(self):
        pass

    def dump_ends(self, path, prefix, coords=None, trj=False):
        if coords is None:
            coords = self.all_mw_coords
        coords = coords.copy()
        coords /= self.m_sqrt
        coords = coords.reshape(-1, len(self.atoms), 3) * BOHR2ANG
        if trj:
            trj_string = make_trj_str(self.atoms, coords, comments=self.all_energies)
            trj_fn = self.get_path_for_fn(f"{prefix}_irc.trj")
            with open(trj_fn, "w") as handle:
                handle.write(trj_string)

        first_coords = coords[0]
        first_fn = self.get_path_for_fn(f"{prefix}_first.xyz")
        with open(first_fn, "w") as handle:
            handle.write(make_xyz_str(self.atoms, first_coords))

        last_coords = coords[-1]
        first_fn = self.get_path_for_fn(f"{prefix}_last.xyz")
        with open(first_fn, "w") as handle:
            handle.write(make_xyz_str(self.atoms, last_coords))

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

        data_dict.update(
            {
                "atoms": np.array(self.atoms, dtype="S"),
                "rms_grad_thresh": np.array(self.rms_grad_thresh),
            }
        )

        if dump_fn is None:
            dump_fn = self.get_path_for_fn(self.dump_fn)

        with h5py.File(dump_fn, "w") as handle:
            for key, val in data_dict.items():
                handle.create_dataset(name=key, dtype=val.dtype, data=val)

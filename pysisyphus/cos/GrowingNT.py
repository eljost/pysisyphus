# [1] https://aip.scitation.org/doi/abs/10.1063/1.1885467
#     Quapp, 2005

import logging
import os
from pathlib import Path

import numpy as np

from pysisyphus.helpers import rms
from pysisyphus.intcoords.helpers import get_weighted_bond_mode


class GrowingNT:
    logger = logging.getLogger("cos")

    def __init__(
        self,
        geom,
        step_len=0.5,
        rms_thresh=1.7e-3,
        r=None,
        final_geom=None,
        between=None,
        bonds=None,
        r_update=True,
        r_update_thresh=1.0,
        stop_after_ts=False,
        require_imag_freq=0.0,
        hessian_at_ts=False,
        out_dir=".",
        dump=True,
    ):
        assert geom.coord_type == "cart"

        self.geom = geom
        self.step_len = step_len
        self.rms_thresh = rms_thresh
        self.final_geom = final_geom
        self.between = between
        self.bonds = bonds
        self.r_update = r_update
        self.r_update_thresh = r_update_thresh
        self.stop_after_ts = stop_after_ts
        self.require_imag_freq = require_imag_freq
        self.hessian_at_ts = hessian_at_ts
        self.out_dir = Path(out_dir)
        self.dump = dump

        if not self.out_dir.exists():
            os.mkdir(self.out_dir)

        self.coord_type = self.geom.coord_type
        if self.final_geom:
            self.converge_to_geom = self.final_geom

        # Determine search direction
        self.r = self.get_r(self.geom, self.final_geom, self.bonds, r)
        self.r_org = self.r.copy()
        # Determine appropriate step_len from between
        if final_geom and self.between:
            self.step_len = np.linalg.norm(final_geom.coords - geom.coords) / (
                self.between + 1
            )

        self._initialized = False
        self.images = [self.geom.copy()]
        self.all_energies = list()
        self.all_real_forces = list()
        self.sp_images = [self.geom.copy()]  # Stationary points
        self.ts_images = list()
        self.min_images = list()

        self.ts_imag_freqs = list()

        # Right now this leads to a gradient calculation in the momement,
        # this object is constructed, which is bad.
        self.initialize()

        if self.dump:
            self.trj_fn = self.get_path("newton_trajectory.trj")

    def get_path(self, fn):
        return self.out_dir / fn

    @staticmethod
    def get_r(geom, final_geom, bonds, r):
        if final_geom:
            r = final_geom - geom
            # self.converge_to_geom = self.final_geom
        elif bonds is not None:
            r = get_weighted_bond_mode(bonds, geom.coords3d)
        # Use 'r' as it is
        elif r is not None:
            pass
        else:
            raise Exception("Please supply either 'r' or 'final_geom'!")

        r = r / np.linalg.norm(r)
        return r

    def log(self, message):
        self.logger.debug(message)

    @property
    def r(self):
        """Parallel/search direction."""
        return self._r

    @property
    def P(self):
        """Projector that keeps perpendicular component."""
        return self._P

    @r.setter
    def r(self, r):
        """Update r and calculate new projector."""
        self._r = r
        self._P = np.eye(self.coords.size) - np.outer(self.r, self.r)

    @property
    def atoms(self):
        return self.geom.atoms

    @property
    def coords(self):
        return self.geom.coords

    @coords.setter
    def coords(self, coords):
        self.geom.coords = coords

    @property
    def cart_coords(self):
        return self.geom.cart_coords

    def grow_image(self):
        self.images[-1] = self.geom.copy()
        # Update coordinates of newly grown image. Try to use Eq. (6) in [1].
        if self._initialized and self.final_geom and self.between:
            m = self.between + 2
            k = len(self.images) - 1
            lambda_ = (m - k) / (m + 1 - k)
            # Adapted from [6] to produce a step instead of new coords
            step = self.coords * (lambda_ - 1) + (1 - lambda_) * self.final_geom.coords
        # If no final image is given we just displace along r
        else:
            step = self.step_len * self.r
        self.coords = self.coords + step

        # Calculate energy and forces at newly grown geometry and append new frontier
        # image.
        real_forces = self.geom.forces
        energy = self.geom.energy
        self.all_energies.append(energy)
        self.all_real_forces.append(real_forces)
        self.images.append(self.geom)

    def initialize(self):
        assert not self._initialized, "GrowingNT.initialize() can only be called once!"
        # Calculation at initial geometry
        init_results = self.geom.get_energy_and_forces_at(self.images[0].coords)
        self.all_energies.append(init_results["energy"])
        self.all_real_forces.append(init_results["forces"])
        # Do initial displacement
        self.grow_image()
        # Indicate the GrowingNT was properly initialized
        self._initialized = True

    def calc_hessian_for(self, other_geom):
        res = self.geom.get_energy_and_cart_hessian_at(other_geom.cart_coords)
        cart_hessian = res["hessian"]
        return cart_hessian

    @property
    def energy(self):
        return self.geom.energy

    @property
    def forces(self):
        forces = self.geom.forces
        perp_forces = self.P.dot(forces)
        return perp_forces

    @property
    def cart_forces(self):
        return self.geom.cart_forces

    def get_energy_at(self, coords):
        return self.geom.get_energy_at(coords)

    def get_energy_and_forces_at(self, coords):
        return self.geom.get_energy_and_forces_at(coords)

    def as_xyz(self):
        return self.geom.as_xyz()

    def clear_passed(self):
        self.passed_min = False
        self.passed_ts = False

    def reparametrize(self):
        """Check if GNT can be grown."""

        # Real, unprojected, forces of the underlying geometry
        real_forces = self.geom.forces
        energy = self.energy
        # Update the last element in all_real_forces and energies with the
        # current values.
        self.all_real_forces[-1] = real_forces
        self.all_energies[-1] = energy

        # See if we can grow the NT, by checking the convergence of the frontier
        # image using projected forces.
        forces = self.forces
        can_grow = rms(forces) <= self.rms_thresh

        if can_grow:
            if self.dump:
                with open(self.trj_fn, "w") as handle:
                    handle.write("\n".join([geom.as_xyz() for geom in self.images]))

            r"""
            Check if we passed a stationary point (SP).
            ^ Energy
            |
            | -3     -1     -2
            |   \    /     /  \
            |    \  /     /    \
            |     -2     -3    -1
            | Minimum       TS
            """
            ae = self.all_energies  # Shortcut
            self.passed_min = len(ae) >= 3 and ae[-3] > ae[-2] < ae[-1]
            self.passed_ts = len(ae) >= 3 and ae[-3] < ae[-2] > ae[-1]
            passed_sp = self.passed_min or self.passed_ts
            if passed_sp:
                sp_image = self.images[-2].copy()
                sp_kind = "Minimum" if self.passed_min else "TS"
                self.sp_images.append(sp_image)
                self.log(
                    f"Passed stationary point! It seems to be a {sp_kind}."
                    f"\n{sp_image.as_xyz()}"
                )
                if self.passed_ts:
                    self.ts_images.append(sp_image)
                    if self.hessian_at_ts:
                        sp_hessian = self.calc_hessian_for(sp_image)
                        nus, *_ = sp_image.get_normal_modes(sp_hessian)
                        self.log(f"First 5 frequencies: {nus[:5]}")
                    if self.require_imag_freq < 0.0:
                        try:
                            sp_hessian
                        except NameError:
                            sp_hessian = self.calc_hessian_for(sp_image)
                        self.ts_imag_freqs.append(sp_image.get_imag_frequencies(sp_hessian))
                elif self.passed_min:
                    self.min_images.append(sp_image)

            # Update direction 'r', if requested
            r_new = self.get_r(self.geom, self.final_geom, self.bonds, self.r)
            r_dot = r_new.dot(self.r)
            r_org_dot = r_new.dot(self.r_org)
            self.log(f"r.dot(r')={r_dot:.6f} r_org.dot(r')={r_org_dot:.6f}")
            if (
                self.r_update
                and (r_org_dot <= self.r_update_thresh)
                and self.passed_min
            ):
                self.r = r_new
                self.log("Updated r")

            # Grow new image
            self.grow_image()
            assert (
                len(self.images) == len(self.all_energies) == len(self.all_real_forces)
            )

        self.did_reparametrization = can_grow
        return can_grow

    def check_convergence(self, *args, **kwargs):
        if len(self.ts_images) == 0:
            return False

        converged = self.stop_after_ts
        if self.require_imag_freq:
            converged = converged and self.ts_imag_freqs[-1][0] <= self.require_imag_freq
        return converged

    def get_additional_print(self):
        if self.did_reparametrization:
            img_num = len(self.images)
            str_ = f"Grew Newton trajectory to {img_num} images."
            if self.passed_min:
                str_ += f" Passed minimum geometry at image {img_num-1}."
            elif self.passed_ts:
                str_ += f" Passed transition state geometry at image {img_num-1}."
        else:
            str_ = None

        self.did_reparametrization = False
        self.clear_passed()

        return str_

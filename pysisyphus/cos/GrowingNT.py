# [1] https://aip.scitation.org/doi/abs/10.1063/1.1885467
#     Quapp, 2005

import logging

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
        update_r=False,
    ):
        assert geom.coord_type == "cart"

        self.geom = geom
        self.step_len = step_len
        self.rms_thresh = rms_thresh
        self.final_geom = final_geom
        self.between = between
        self.bonds = bonds
        self.update_r = update_r

        self.coord_type = self.geom.coord_type
        if self.final_geom:
            self.converge_to_geom = self.final_geom

        # Determine search direction
        self.r = self.get_r(self.geom, self.final_geom, self.bonds, r)
        # Determine appropriate step_len from between
        if final_geom and self.between:
            self.step_len = np.linalg.norm(final_geom.coords - geom.coords) / (
                self.between + 1
            )

        self.images = [self.geom.copy()]
        self.sp_images = [self.geom.copy()]  # Stationary points
        self.all_energies = list()
        self.all_forces = list()

        # Do initial displacement towards final_geom along r
        step = self.step_len * self.r
        self.geom.coords += step

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

    def as_xyz(self):
        return self.geom.as_xyz()

    def reparametrize(self):
        # TODO: use already calculated quantities to decide this
        real_forces = self.geom.forces  # Unprojected forces
        energy = self.energy
        if rms(real_forces) <= self.rms_thresh:
            self.sp_images.append(self.geom.copy())

        forces = self.forces  # Projected forces
        reparametrized = rms(forces) <= self.rms_thresh
        # Check if we passed a stationary point by comparing energies
        if reparametrized and len(self.all_energies) > 2:
            prev_prev_energy, prev_energy = self.all_energies[-3:-1]
            self.prev_energy_increased = prev_energy > prev_prev_energy
            sp_passed = energy < prev_energy
            if sp_passed:
                self.log(f"Passed stationary point!\n{self.geom.as_xyz()}")
                self.sp_images.append(self.geom.copy())
            if sp_passed and self.update_r:
                self.r = self.get_r(self.geom, self.final_geom, self.bonds, self.r)
            if len(self.images) == 5:
                self.r = self.get_r(self.geom, self.final_geom, self.bonds, self.r)

        if reparametrized:
            self.images.append(self.geom.copy())
            self.all_forces.append(real_forces)
            self.all_energies.append(energy)
            step = self.step_len * self.r
            self.coords = self.coords + step

        self.did_reparametrization = reparametrized
        return reparametrized

    def get_additional_print(self):
        if self.did_reparametrization:
            str_ = "Grew Newton trajectory."
        else:
            str_ = None
        self.did_reparametrization = False

        return str_

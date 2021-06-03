import numpy as np

from pysisyphus.helpers import rms


class GrowingNT:
    def __init__(
        self,
        geom,
        step_len=0.5,
        rms_thresh=1.7e-3,
        r=None,
        final_geom=None,
        between=None,
        update_r=False,
    ):
        assert geom.coord_type == "cart"

        self.geom = geom
        self.step_len = step_len
        self.rms_thresh = rms_thresh

        self.coord_type = self.geom.coord_type
        if final_geom:
            r = final_geom - geom
            r = r / np.linalg.norm(r)
            self.converge_to_geom = final_geom
        elif r is not None:
            pass
        else:
            raise Exception("Please supply either 'r' or 'final_geom'!")

        if final_geom and between:
            # Determine appropriate step_len from between
            self.step_len = np.linalg.norm(final_geom.coords - geom.coords) / (
                between + 1
            )

        self.r = r
        self.P = np.eye(self.coords.size) - np.outer(r, r)

        self.images = [self.geom.copy()]
        self.sp_images = [self.geom.copy()]  # Stationary points
        self.all_energies = list()
        self.all_forces = list()
        # Do initial displacement towards final_geom along r
        step = self.step_len * self.r
        self.geom.coords += step

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
        real_forces = self.geom.forces  # Unprojected forces
        energy = self.energy
        if rms(real_forces) <= self.rms_thresh:
            self.sp_images.append(self.geom.copy())

        forces = self.forces  # Projected forces
        reparametrized = rms(forces) <= self.rms_thresh
        if reparametrized:
            self.images.append(self.geom.copy())
            self.all_forces.append(forces)
            self.all_energies.append(energy)
            step = self.step_len * self.r
            self.coords = self.coords + step
        return reparametrized

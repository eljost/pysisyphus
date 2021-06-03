import numpy as np


class GrowingNT:
    def __init__(
        self, geom, step_len=0.5, rms_thresh=1e-3, r=None, final_geom=None, between=10
    ):
        assert geom.coord_type == "cart"

        self.geom = geom
        self.step_len = step_len
        self.rms_thresh = rms_thresh

        self.coord_type = self.geom.coord_type
        if final_geom is not None:
            r = final_geom - geom
            r = r / np.linalg.norm(r)
            self.converge_to_geom = final_geom
        elif r is not None:
            pass
        else:
            raise Exception("Please supply either 'r' or 'final_geom'!")

        self.r = r
        self.P = np.eye(self.coords.size) - np.outer(r, r)

        self.images = [self.geom.copy()]
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
        norm = np.linalg.norm(self.forces)
        reparametrized = norm <= self.rms_thresh
        if reparametrized:
            self.images.append(self.geom.copy())
            step = self.step_len * self.r
            new_coords = self.coords + step  # self.step * self.r
            self.coords = new_coords
            print("Reparametrized")
        return reparametrized

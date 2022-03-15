import warnings

import numpy as np
from numpy.typing import NDArray


class CartesianCoords:
    def __init__(
        self,
        atoms,
        coords3d,
        masses,
        freeze_atoms=None,
        *,
        mass_weighted=False,
        **kwargs,
    ):
        for key, val in kwargs.items():
            # f-string don't seem to work when pytest reports the warnings after a test
            # run, so we construct to warning before it is issued.
            msg = (
                f"Keyword '{key}': '{val}' is not supported by this coordinate system!"
            )
            warnings.warn(msg)
        self.atoms = atoms
        self._coords3d = coords3d
        self.atom_num = len(self.atoms)
        self.masses = masses
        if freeze_atoms is None:
            freeze_atoms = list()
        self.freeze_atoms = freeze_atoms
        self.mass_weighted = mass_weighted

        self.move_mask = np.full(self.atom_num, True, dtype=bool)
        self.move_mask[self.freeze_atoms] = False
        self.move_mask_rep = np.repeat(self.move_mask, 3)
        self.zero_vec = np.zeros_like(coords3d)

    @property
    def masses(self) -> NDArray:
        return self._masses

    @property
    def masses_sqrt(self) -> NDArray:
        return self._masses_sqrt

    @property
    def inv_masses_rep_sqrt(self) -> NDArray:
        return self._inv_masses_rep_sqrt

    @masses.setter
    def masses(self, masses):
        assert len(masses) == self.atom_num
        masses = np.array(masses, dtype=float)
        self._masses = masses
        # Also precalculate masses' square roots for mass-weighting
        self._masses_sqrt = np.sqrt(self.masses)
        # and its inverse
        self._inv_masses_rep_sqrt = 1 / np.repeat(self.masses_sqrt, 3)

    @property
    def coords(self):
        coords3d = self.coords3d.copy()
        if self.mass_weighted:
            coords3d *= self.masses_sqrt[:, None]
        coords3d = coords3d[self.move_mask]
        return coords3d.flatten()

    def transform_forces(self, cart_forces):
        forces = cart_forces.reshape(-1, 3)
        if self.mass_weighted:
            forces /= self.masses_sqrt[:, None]
        forces = forces[self.move_mask]
        return forces.flatten()

    def transform_hessian(self, cart_hessian, gradient=None):
        if self.mass_weighted:
            M_inv = np.diag(self.inv_masses_rep_sqrt)
            cart_hessian = M_inv @ cart_hessian @ M_inv
        cart_hessian = cart_hessian[self.move_mask_rep][:, self.move_mask_rep]
        return cart_hessian

    def transform_int_step(self, step, update_constraints=False, pure=False):
        if update_constraints:
            raise Exception("update_constraints is currently ignored!")
        full_step = self.zero_vec.copy()
        full_step[self.move_mask] = step.reshape(-1, 3)
        if self.mass_weighted:
            full_step /= self.masses_sqrt[:, None]
        if not pure:
            self.coords3d += full_step
        return full_step.flatten()

    def project_hessian(self, hessian):
        return hessian

    @property
    def coords3d(self):
        return self._coords3d

    @coords3d.setter
    def coords3d(self, coords3d):
        self._coords3d = coords3d.reshape(-1, 3)

    @property
    def typed_prims(self):
        return list()


class MWCartesianCoords(CartesianCoords):
    def __init__(self, *args, **kwargs):
        kwargs["mass_weighted"] = True
        super().__init__(*args, **kwargs)

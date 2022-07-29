from typing import List, Optional
import tempfile
import warnings

from joblib import Memory
import numpy as np
from numpy.typing import NDArray

from pysisyphus.intcoords.Coords import CoordSys
from pysisyphus.intcoords.Stretch import Stretch


memory = Memory(location=str(tempfile.TemporaryDirectory()), verbose=False)


def rattle(
    bonds, constrain_bonds, coords3d, full_step, tol=1.0e-2, max_cycles=50, lengths=None
):
    tol_sq = tol ** 2

    # Original bond lengths
    if lengths is None:
        lengths = np.array([bond.calculate(coords3d) for bond in bonds])
    lengths_sq = lengths ** 2

    bond_vecs = np.diff(coords3d[constrain_bonds], axis=1).reshape(-1, 3)

    # Unconstrained position update
    coords3d_updated = coords3d + full_step

    # Yields updated positions for t+dt and half-updated velocities for t+dt/2
    for i in range(max_cycles):
        corrected = False
        for j, (atom_1, atom_2) in enumerate(constrain_bonds):
            bond_vec_updated = coords3d_updated[atom_1] - coords3d_updated[atom_2]
            bond_length_sq = bond_vec_updated.dot(bond_vec_updated)
            ref_length_sq = lengths_sq[j]
            diff_sq = ref_length_sq - bond_length_sq
            print(f"{i=:03d}, {j=:>03d}, {diff_sq: >12.6f}")

            # We have to correct coordinates and velocities if the deviations is
            # above the tolerance.
            if abs(diff_sq) <= (ref_length_sq * tol_sq):
                continue

            corrected = True
            # Otherwise try to satisfy the constraint by calculating the
            # approximate lagrange multiplier 'g'.
            dot = bond_vec_updated.dot(bond_vecs[j])
            # TODO: test for constraint failure
            g = diff_sq / (2 * dot)

            # Update positions to satify constraint
            atom_1_factor = g * bond_vecs[j]
            atom_2_factor = g * bond_vecs[j]
            coords3d_updated[atom_1] += atom_1_factor
            coords3d_updated[atom_2] -= atom_2_factor

        # Stop macro-iterations when no correction was done
        if not corrected:
            print(f"RATTLE_1 finished after {i+1} macro cycle(s)!")
            break
        print()
    else:
        raise Exception("First part of RATTLE did not converge!")

    return coords3d_updated


def rattle(
    bonds, constrain_bonds, coords3d, full_step, tol=1e-10, max_cycles=500, lengths=None
):
    old = coords3d  # atoms.positions
    new = coords3d + full_step  #

    for i in range(100):#self.maxiter):
        converged = True
        for j, ab in enumerate(constrain_bonds):  # self.pairs):
            a = ab[0]
            b = ab[1]
            cd = lengths[j]
            d0 = old[a] - old[b]
            d1 = new[a] - new[b]# - d0# + d0
            x = 0.5 * (cd ** 2 - np.dot(d1, d1)) / np.dot(d0, d1)
            # print(f"{i=}, {x=}")
            if abs(x) > tol:
                new[a] += x * d0
                new[b] -= x * d0
                converged = False
        if converged:
            break

    return new  # coords3d_updated


@memory.cache
def get_orth_projector(bonds, coords3d):
    I = np.eye(coords3d.size)
    _, grads = zip(*[bond.calculate(coords3d, gradient=True) for bond in bonds])
    wilson_vecs = np.array(grads).T
    wilson_vecs, _ = np.linalg.qr(wilson_vecs)
    P = wilson_vecs @ wilson_vecs.T
    D = I - P
    return D


class CartesianCoords(CoordSys):
    def __init__(
        self,
        atoms,
        coords3d: NDArray,
        masses,
        freeze_atoms=None,
        constrain_bonds=None,
        *,
        mass_weighted=False,
        **kwargs,
    ):
        for key, val in kwargs.items():
            # f-string don't seem to work when pytest reports the warnings after a test
            # run, so we construct to warning before it is issued.
            msg = (
                # Also printing the value is too much
                # f"Keyword '{key}': '{val}' is not supported by this coordinate system!"
                f"Keyword '{key}' is not supported by this coordinate system!"
            )
            warnings.warn(msg)
        self.atoms = atoms
        self._coords3d = coords3d
        self.atom_num = len(self.atoms)
        self.masses = masses
        if freeze_atoms is None:
            freeze_atoms = list()
        self.freeze_atoms = freeze_atoms
        if constrain_bonds is None:
            constrain_bonds = list()
        self.constrain_bonds = np.array(constrain_bonds, dtype=int)
        self.mass_weighted = mass_weighted

        self.move_mask = np.full(self.atom_num, True, dtype=bool)
        self.move_mask[self.freeze_atoms] = False
        self.move_mask_rep = np.repeat(self.move_mask, 3)
        self.zero_vec = np.zeros_like(coords3d)

        self.bonds = [Stretch(inds) for inds in self.constrain_bonds]
        self.target_bond_lens = self.bond_lengths()

    def orth_projector(self):
        return get_orth_projector(self.bonds, self.coords3d)

    def bond_lengths(self, for_coords=None):
        if for_coords is None:
            for_coords = self.coords3d
        return np.array([bond.calculate(for_coords) for bond in self.bonds])

    @property
    def masses(self) -> NDArray:
        return self._masses

    @masses.setter
    def masses(self, masses: NDArray):
        assert len(masses) == self.atom_num
        masses = np.array(masses, dtype=float)
        self._masses = masses
        # Also precalculate masses' square roots for mass-weighting
        self._masses_sqrt = np.sqrt(self.masses)
        # and its inverse
        self._inv_masses_rep_sqrt = 1 / np.repeat(self.masses_sqrt, 3)

    @property
    def masses_sqrt(self) -> NDArray:
        return self._masses_sqrt

    @property
    def inv_masses_rep_sqrt(self) -> NDArray:
        return self._inv_masses_rep_sqrt

    @property
    def coords(self) -> NDArray:
        coords3d = self.coords3d.copy()
        if self.mass_weighted:
            coords3d *= self.masses_sqrt[:, None]
        coords3d = coords3d[self.move_mask]
        return coords3d.flatten()

    def transform_forces(self, cart_forces: NDArray) -> NDArray:
        if self.constrain_bonds.size > 0:
            cart_forces = self.orth_projector() @ cart_forces
        forces = cart_forces.reshape(-1, 3)
        if self.mass_weighted:
            forces /= self.masses_sqrt[:, None]
        forces = forces[self.move_mask]
        return forces.flatten()

    def transform_hessian(
        self, cart_hessian: NDArray, int_gradient: Optional[NDArray] = None
    ):
        if len(self.constrain_bonds) > 0:
            D = self.orth_projector()
            cart_hessian = D @ cart_hessian @ D
        if self.mass_weighted:
            M_inv = np.diag(self.inv_masses_rep_sqrt)
            cart_hessian = M_inv @ cart_hessian @ M_inv
        cart_hessian = cart_hessian[self.move_mask_rep][:, self.move_mask_rep]
        return cart_hessian

    def transform_int_step(
        self,
        step: NDArray,
        update_constraints: Optional[bool] = False,
        pure: Optional[bool] = False,
    ) -> NDArray:
        if update_constraints:
            raise Exception("update_constraints is currently ignored!")
        full_step = self.zero_vec.copy()
        full_step[self.move_mask] = step.reshape(-1, 3)
        if self.mass_weighted:
            full_step /= self.masses_sqrt[:, None]

        # c3d_updated = rattle(self.bonds, self.constrain_bonds, self.coords3d, full_step)
        c3d_updated = rattle(
            self.bonds,
            self.constrain_bonds,
            self.coords3d,
            full_step,
            lengths=self.target_bond_lens,
        )
        full_step = c3d_updated - self.coords3d
        # if self.constrain_bonds.size > 0:
        # full_step = (self.orth_projector() @ full_step.flatten()).reshape(-1, 3)
        # full_step = full
        # if len(self.constrain_bonds) > 0:
        # # if False:
        # # cur_bond_lens = self.bond_lengths()
        # new_bond_lens = self.bond_lengths(self.coords3d + full_step)
        # diff = self.target_bond_lens - new_bond_lens
        # bond_vecs = np.diff(self.coords3d[self.constrain_bonds], axis=1).reshape(-1, 3)
        # feasibility_corr = diff[:, None] * bond_vecs / new_bond_lens[:, None]
        # full_step -= feasibility_corr
        # new_coords3d = self.coords3d + full_step
        # for (from_, to_), target in zip(self.constrain_bonds, self.target_bond_lens):
        # new_bond_len = np.linalg.norm(new_coords3d[from_] - new_coords3d[to_])
        # diff = target - new_bond_len

        if not pure:
            self.coords3d += full_step
        return full_step.flatten()

    def project_hessian(self, hessian):
        if len(self.constrain_bonds) > 0:
            D = self.orth_projector()
            # When atoms are frozen hessian will have a smaller shape, compared to
            # cart_hessian.
            D = D[self.move_mask_rep][:, self.move_mask_rep]
            hessian = D @ hessian @ D
        return hessian

    @property
    def coords3d(self) -> NDArray:
        return self._coords3d

    @coords3d.setter
    def coords3d(self, coords3d: NDArray):
        self._coords3d = coords3d.reshape(-1, 3)

    @property
    def typed_prims(self) -> List:
        return list()


class MWCartesianCoords(CartesianCoords):
    def __init__(self, *args, **kwargs):
        kwargs["mass_weighted"] = True
        super().__init__(*args, **kwargs)

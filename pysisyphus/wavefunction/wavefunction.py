from typing import Optional

import numpy as np
from numpy.typing import NDArray

from pysisyphus.elem_data import nuc_charges_for_atoms
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.shells import Shells
from pysisyphus.wavefunction.helpers import BFType


class Wavefunction:
    def __init__(
        self,
        atoms,
        coords,
        charge: int,
        mult: int,
        unrestricted: bool,
        occ: int,
        C: NDArray,
        bf_type: BFType,
        shells: Optional[Shells] = None,
    ):
        self.atoms = atoms
        self.coords = np.array(coords).flatten()
        assert 3 * len(self.atoms) == len(self.coords)

        self.charge = charge
        self.mult = mult
        self.unrestricted = unrestricted
        self.occ = occ
        if not self.unrestricted:
            assert self.occ[0] == self.occ[1]
        self.C = np.array(C)
        # Always keep α and β coefficients separate. If we get only one set
        # of MO coefficients (one square matrix with ndim == 2) we assume a
        # restricted calculation and use the same set for alpha and beta electrons.
        if C.ndim == 2:
            self.C = np.array((self.C, self.C.copy()))
        self.bf_type = bf_type
        self.shells = shells

    @property
    def atom_num(self):
        return len(self.atoms)

    @property
    def nuc_charges(self):
        return nuc_charges_for_atoms(self.atoms)

    @property
    def masses(self):
        return np.ones(len(self.atoms))

    @staticmethod
    @file_or_str(".json")
    def from_orca_json(text):
        from pysisyphus.io.orca import wavefunction_from_json

        wf = wavefunction_from_json(text)
        return wf

    @property
    def C_occ(self):
        occ_a, occ_b = self.occ
        C_a, C_b = self.C
        return np.array((C_a[:, :occ_a], C_b[:, :occ_b]))

    @property
    def P(self):
        P = np.array([C_occ @ C_occ.T for C_occ in self.C_occ])
        # Restricted density matrix: P = 2 * C_occ @ C_occ.T
        return P

    @property
    def S(self):
        # Check what type of basis functions we are using.
        return {
            BFType.CARTESIAN: self.shells.S_cart,
            BFType.PURE_SPHERICAL: self.shells.S_sph,
        }[self.bf_type]

    def S_with(self, other):
        return {
            BFType.CARTESIAN: self.shells.get_S_cart(other),
            BFType.PURE_SPHERICAL: self.shells.get_S_sph(other),
        }[self.bf_type]

    @property
    def ao_centers(self):
        ao_centers = list(
            {
                BFType.CARTESIAN: self.shells.cart_ao_centers,
                BFType.PURE_SPHERICAL: self.shells.sph_ao_centers,
            }[self.bf_type]
        )
        return ao_centers

    def get_origin(self, kind="coc"):
        kind = kind.lower()
        assert kind in ("com", "coc")

        if kind == "com":
            raise Exception("Masses are not yet implemented!")

        coords3d = self.coords.reshape(-1, 3).copy()
        _factors = {
            "coc": self.nuc_charges,  # Center of charge
            "com": self.masses,  # Center of mass
        }
        factors = _factors[kind]
        tot_factors = factors.sum()

        return np.sum(coords3d * factors[:, None], axis=0) / tot_factors

    def dipole_ints(
        self, origin: Optional[NDArray] = None, kind="coc"
    ) -> NDArray[float]:
        if origin is None:
            origin = self.get_origin(kind=kind)

        dipole_ints = {
            BFType.CARTESIAN: self.shells.get_dipole_ints_cart,
            BFType.PURE_SPHERICAL: self.shells.get_dipole_ints_sph,
        }[self.bf_type](origin)
        return dipole_ints

    def dipole_moment(
        self, origin: Optional[NDArray] = None, kind="coc"
    ) -> NDArray[float]:
        # Make a copy, as the coordinates must be give w.r.t. the origin, which
        # may be != (0., 0., 0.).
        coords3d = self.coords.reshape(-1, 3).copy()
        nuc_charges = self.nuc_charges

        if origin is None:
            origin = self.get_origin(kind=kind)

        origin = np.array(origin, dtype=float)
        coords3d -= origin[None, :]
        dipole_ints = self.dipole_ints(origin)
        electronic = np.sum(
            [np.einsum("xij,ji->x", dipole_ints, P) for P in self.P], axis=0
        )
        nuclear = np.einsum("i,ix->x", nuc_charges, coords3d)
        molecular = nuclear - electronic

        return molecular

    def transition_dipole_moment(self, trans_dens):
        origin = self.get_origin(kind="coc")
        dipole_ints = self.dipole_ints(origin)
        """
        This assert could be lifted/removed to handle an array of multiple transition
        densities at once, e.g., to avoid repeated recalculation of the dipole integrals.
        """
        assert trans_dens.ndim == 2
        occ, _ = trans_dens.shape
        Ca, _ = self.C
        trans_dens_ao = Ca[:occ].T @ trans_dens @ Ca[occ:]
        tdm = np.einsum("xij,ji->x", dipole_ints, trans_dens_ao)
        return tdm

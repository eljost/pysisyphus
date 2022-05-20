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
        self.C = C
        self.bf_type = bf_type
        self.shells = shells

    @property
    def atom_num(self):
        return len(self.atoms)

    @property
    def nuc_charges(self):
        return nuc_charges_for_atoms(self.atoms)

    @staticmethod
    @file_or_str(".json")
    def from_orca_json(text):
        from pysisyphus.io.orca import wavefunction_from_json

        wf = wavefunction_from_json(text)
        return wf

    @property
    def C_occ(self):
        C = self.C
        occ = self.occ
        if self.unrestricted:
            C_occ = [C[i, :, :occ] for i, occ in enumerate(occ)]
        else:
            C_occ = self.C[:, : occ[0]]
        return C_occ

    @property
    def P(self):
        if self.unrestricted:
            P = np.array([C_occ @ C_occ.T for C_occ in self.C_occ])
        else:
            C_occ = self.C_occ
            P = 2 * C_occ @ C_occ.T
        return P

    @property
    def S(self):
        # Check what type of integrals are
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

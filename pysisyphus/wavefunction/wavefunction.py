# [1] https://doi.org/10.1063/1.4937410
#     The consequences of improperly describing oscillator strengths
#     beyond the electric dipole approximation
#     Lestrange, Egidi, Li, 2015
# [2] https://doi.org/10.1016/0009-2614(95)01036-9
#     Ab initio calculation and display of the rotary strength tensor in
#     the random phase approximation. Method and model studies.
#     Pedersen, Hansen, 1995

from typing import Optional, Tuple

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
        occ: Tuple[int],
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
        # MOs are stored in columns
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

    def transition_dipole_moment(
        self, trans_dens: NDArray[float], renorm: bool = True
    ) -> NDArray[float]:
        """
        Eq. (33) in [1].

        X+Y for hermitian operators, X-Y for anti-hermitian operators.
        """
        origin = self.get_origin(kind="coc")
        dipole_ints = self.dipole_ints(origin)
        trans_dens = np.atleast_3d(trans_dens)

        # The 'tdm' formula in get_tdm assumes a norm of 0.5 for the transition
        # density matrices.
        if renorm:
            norms = np.linalg.norm(trans_dens, axis=(1, 2))
            factors = 0.5 / norms
        else:
            factors = np.ones(trans_dens.shape[0])

        def get_tdm(C: NDArray[float], occ: NDArray[int]) -> NDArray[float]:
            # Transform transition density matrix from MO to AO representation.
            # shape in MO basis: (occ, virt), shape in AO basis (nbf, nbf).
            # trans_dens_ao = C[None, :, :occ] @ trans_dens #@ C[None, :, occ:].T
            trans_dens_ao = C[None, :, :occ] @ trans_dens @ C[:, occ:].T
            # Then contract with dipole integrals. See Eq. (18) in [2].
            # 'factors' is used to take the normalization (or lack of) of the
            # transition density matrix into account.
            #
            # x   : Cartesian direction
            # n   : State index
            # i, j: AO indices
            tdm = 2 ** 0.5 * np.einsum(
                "xij,n,nji->nx", dipole_ints, factors, trans_dens_ao
            )
            return tdm

        # Transitions between α -> α and β -> β
        tdms = [get_tdm(C, occ) for C, occ in zip(self.C, self.occ)]
        tdm = np.sum(tdms, axis=0)
        return tdm

    def oscillator_strength(
        self, exc_ens: NDArray[float], trans_moms: NDArray[float]
    ) -> NDArray[float]:
        exc_ens = np.atleast_1d(exc_ens)
        trans_moms = np.atleast_2d(trans_moms)
        fosc = 2 / 3 * (exc_ens[:, None] * trans_moms ** 2).sum(axis=0)
        return fosc

    def transition_dipole_moment(
        self,
        trans_dens: NDArray[float],
        renorm: bool = True,
        restricted=False,
    ) -> NDArray[float]:
        """
        Eq. (33) in [1].

        X+Y for hermitian operators, X-Y for anti-hermitian operators.
        """
        origin = self.get_origin(kind="coc")
        dipole_ints = self.dipole_ints(origin)

        # Bring into '(nstates, a/b, act, virt)' shape.
        if restricted:
            *nstates, act, virt = trans_dens.shape
            assert len(nstates) in (0, 1)
            nstates = nstates[0]
            if nstates == 0:
                nstates = 1
            # Repeat last two axes
            _ = np.empty((nstates, 2, act, virt))
            _[:, 0] = trans_dens
            _[:, 1] = trans_dens
            # trans_dens is now 4d
            trans_dens = _

        # Deal with unrestricted transition density matrices
        if trans_dens.ndim == 3:
            trans_dens = trans_dens[None, :]

        assert trans_dens.ndim == 4

        def get_tdm(
            C: NDArray[float], occ: NDArray[int], trans_dens: NDArray[float]
        ) -> NDArray[float]:
            # The 'tdm' formula in get_tdm assumes a norm of 0.5 for the transition
            # density matrices.
            if renorm:
                norms = np.linalg.norm(trans_dens, axis=(1, 2))
                factors = 0.5 / norms
            else:
                factors = np.ones(trans_dens.shape[0])

            # Transform transition density matrix from MO to AO representation.
            # shape in MO basis: (occ, virt), shape in AO basis (nbf, nbf).
            # trans_dens_ao = C[None, :, :occ] @ trans_dens #@ C[None, :, occ:].T
            trans_dens_ao = C[None, :, :occ] @ trans_dens @ C[:, occ:].T

            # Then contract with dipole integrals. See Eq. (18) in [2].
            # 'factors' is used to take the normalization (or lack of) of the
            # transition density matrix into account.
            #
            # x   : Cartesian direction
            # n   : State index
            # i, j: AO indices
            tdm = 2 ** 0.5 * np.einsum(
                "xij,n,nji->nx", dipole_ints, factors, trans_dens_ao
            )
            return tdm

        # Transitions between α -> α and β -> β
        C_a, C_b = self.C
        occ_a, occ_b = self.occ
        tdms_a = get_tdm(C_a, occ_a, trans_dens[:, 0])
        tdms_b = get_tdm(C_b, occ_b, trans_dens[:, 1])
        tdms = tdms_a + tdms_b
        return tdms

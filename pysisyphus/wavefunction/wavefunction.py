# [1] https://doi.org/10.1063/1.4937410
#     The consequences of improperly describing oscillator strengths
#     beyond the electric dipole approximation
#     Lestrange, Egidi, Li, 2015
# [2] https://doi.org/10.1016/0009-2614(95)01036-9
#     Ab initio calculation and display of the rotary strength tensor in
#     the random phase approximation. Method and model studies.
#     Pedersen, Hansen, 1995
# [3] https://pubs.acs.org/doi/pdf/10.1021/j100180a030.
#     Toward a Systematic Molecular Orbital Theory for Excited States
#     Foresman, Head-Gordon, Pople, Frisch, 1991

from typing import Literal, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from pysisyphus.elem_data import nuc_charges_for_atoms, MASS_DICT
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.shells import Shells
from pysisyphus.wavefunction.helpers import BFType


Center = Literal["coc", "com"]


class Wavefunction:
    def __init__(
        self,
        atoms: Tuple[str],
        coords: NDArray[float],
        charge: int,
        mult: int,
        unrestricted: bool,
        occ: Tuple[int],
        C: NDArray[float],
        bf_type: BFType,
        shells: Optional[Shells] = None,
    ):
        self.atoms = atoms
        self.coords = np.array(coords).flatten()
        assert 3 * len(self.atoms) == len(self.coords)

        self.charge = charge
        self.mult = mult
        if self.mult != 1:
            unrestricted = True
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

        self._masses = np.array([MASS_DICT[atom.lower()] for atom in self.atoms])

    @property
    def atom_num(self):
        return len(self.atoms)

    @property
    def masses(self):
        return self._masses

    @property
    def nuc_charges(self):
        return nuc_charges_for_atoms(self.atoms)

    @property
    def masses(self):
        return np.ones(len(self.atoms))

    @property
    def mo_num(self):
        return self.C.shape[1]

    @staticmethod
    @file_or_str(".json")
    def from_orca_json(text):
        from pysisyphus.io.orca import wavefunction_from_json

        wf = wavefunction_from_json(text)
        return wf

    @staticmethod
    @file_or_str(".json")
    def from_orca_molden(text):
        from pysisyphus.io.orca import wavefunction_from_molden

        wf = wavefunction_from_molden(text)
        return wf

    @property
    def C_occ(self):
        occ_a, occ_b = self.occ
        C_a, C_b = self.C
        return C_a[:, :occ_a], C_b[:, :occ_b]

    @property
    def C_virt(self):
        occ_a, occ_b = self.occ
        C_a, C_b = self.C
        return C_a[:, occ_a:], C_b[:, occ_b:]

    @property
    def P(self):
        return [C_occ @ C_occ.T for C_occ in self.C_occ]

    def P_exc(self, trans_dens):
        """
        Eqs. (2.25) and (2.26) in [3].
        """
        trans_dens *= 2 ** 0.5 / 2
        occ_a, occ_b = self.occ
        assert occ_a == occ_b
        occ = occ_a
        occupations = np.zeros(self.mo_num)
        occupations[:occ_a] = 2
        P = np.diag(occupations)
        dP_oo = -np.einsum("ia,ja->ij", trans_dens, trans_dens)
        dP_vv = np.einsum("ia,ic->ac", trans_dens, trans_dens)
        P[:occ, :occ] += 2 * dP_oo
        P[occ:, occ:] += 2 * dP_vv
        C, _ = self.C
        # The density matric is currently still in the MO basis. Transform
        # it to the AO basis and return.
        return C @ P @ C.T

    @property
    def S(self):
        # Check what type of basis functions we are using.
        return {
            BFType.CARTESIAN: lambda: self.shells.S_cart,
            BFType.PURE_SPHERICAL: lambda: self.shells.S_sph,
        }[self.bf_type]()

    def S_with_shells(self, other_shells):
        return {
            BFType.CARTESIAN: lambda: self.shells.get_S_cart(other_shells),
            BFType.PURE_SPHERICAL: lambda: self.shells.get_S_sph(other_shells),
        }[self.bf_type]()

    def S_with(self, other):
        other_shells = other.shells
        return self.S_with_shells(other_shells)

    def S_MO_with(self, other):
        assert (not self.unrestricted) and (
            self.mult == 1
        ), "Currently only Cα is considered!"
        S_AO = self.S_with(other)
        C = self.C[0]
        C_other = other.C[0]
        return C.T @ S_AO @ C_other

    @property
    def ao_centers(self) -> List[int]:
        return list(
            {
                BFType.CARTESIAN: lambda: self.shells.cart_ao_centers,
                BFType.PURE_SPHERICAL: lambda: self.shells.sph_ao_centers,
            }[self.bf_type]()
        )

    @property
    def ao_center_map(self) -> dict[int, List[int]]:
        ao_center_map = dict()
        for i, aoc in enumerate(self.ao_centers):
            ao_center_map.setdefault(aoc, list()).append(i)
        return ao_center_map

    def get_origin(self, kind="coc"):
        kind = kind.lower()
        assert kind in ("com", "coc")

        coords3d = self.coords.reshape(-1, 3).copy()
        _factors = {
            "coc": self.nuc_charges,  # Center of charge
            "com": self.masses,  # Center of mass
        }
        factors = _factors[kind]
        tot_factors = factors.sum()

        return np.sum(coords3d * factors[:, None], axis=0) / tot_factors

    def dipole_ints(
        self, origin: Optional[NDArray] = None, kind: Center = "coc"
    ) -> NDArray[float]:
        if origin is None:
            origin = self.get_origin(kind=kind)

        dipole_ints = {
            BFType.CARTESIAN: lambda *args: self.shells.get_dipole_ints_cart(*args),
            BFType.PURE_SPHERICAL: lambda *args: self.shells.get_dipole_ints_sph(*args),
        }[self.bf_type](origin)
        return dipole_ints

    def dipole_moment(
        self,
        P: NDArray[float] = None,
        origin: Optional[NDArray] = None,
        kind: Center = "coc",
    ) -> NDArray[float]:
        # Make a copy, as the coordinates must be give w.r.t. the origin, which
        # may be != (0., 0., 0.).
        coords3d = self.coords.reshape(-1, 3).copy()
        nuc_charges = self.nuc_charges

        if P is None:
            P = self.P

        if origin is None:
            origin = self.get_origin(kind=kind)

        # origin = (-1.487808,  3.034774,  0.292895)
        origin = np.array(origin, dtype=float)
        coords3d -= origin[None, :]
        dipole_ints = self.dipole_ints(origin)
        electronic = np.sum(
            [np.einsum("xij,ji->x", dipole_ints, P_) for P_ in P], axis=0
        )
        nuclear = np.einsum("i,ix->x", nuc_charges, coords3d)
        molecular = nuclear - electronic

        return molecular

    def transition_dipole_moment(
        self,
        # trans_dens: NDArray[float],
        trans_dens_a: NDArray[float],
        trans_dens_b: NDArray[float] = None,
    ) -> NDArray[float]:
        """
        Transition dipole moments from transition density matrices.

        Eq. (33) in [1].

        Parameters
        ----------
        trans_dens_a
            Numpy array containing transition density matrices of shape
            (nstates, occ, virt) or (occ, virt) in the α-electron space.
        trans_dens_b
            See trans_dens_a. Optional. If not provided, trans_dens_a is used.

        Returns
        -------
        tdms
            Numpy array of shape (nstates, 3) containing the Cartesian transition
            dipole moments in the length gauge.
        """

        origin = self.get_origin(kind="coc")
        dipole_ints = self.dipole_ints(origin)

        # Deal with unrestricted transition density matrices that contain
        # only one state.
        def atleast_3d(tden):
            if tden.ndim == 2:
                tden = tden[None, :]
            return tden

        trans_dens_a = atleast_3d(trans_dens_a)
        if trans_dens_b is None:
            trans_dens_a = 1 / 2 ** 0.5 * trans_dens_a
            trans_dens_b = trans_dens_a
        trans_dens_b = atleast_3d(trans_dens_b)

        # Expected shape: (nstates, occ, virt)
        assert trans_dens_a.ndim == trans_dens_b.ndim == 3

        def get_tdm(
            C: NDArray[float], occ: NDArray[int], trans_dens: NDArray[float]
        ) -> NDArray[float]:
            dipole_ints_mo = np.einsum(
                "jl,ilm,mk->ijk", C[:, :occ].T, dipole_ints, C[:, occ:]
            )

            # Then contract with dipole integrals. See Eq. (18) in [2].
            #
            # j : Cartesian direction
            # k : occ. MO space
            # l : virt. MO space
            # i : Excited state number
            tdm = np.einsum("jkl,ikl->ij", dipole_ints_mo, trans_dens)
            return tdm

        # Transitions between α -> α and β -> β
        C_a, C_b = self.C
        occ_a, occ_b = self.occ
        tdms_a = get_tdm(C_a, occ_a, trans_dens_a)
        tdms_b = get_tdm(C_b, occ_b, trans_dens_b)
        tdms = tdms_a + tdms_b
        return tdms

    def oscillator_strength(
        self, exc_ens: NDArray[float], trans_moms: NDArray[float]
    ) -> NDArray[float]:
        exc_ens = np.atleast_1d(exc_ens)
        if trans_moms.ndim == 1:
            trans_moms = trans_moms[None, :]
        fosc = 2 / 3 * exc_ens * (trans_moms ** 2).sum(axis=1)
        return fosc

    def as_geom(self):
        return Geometry(self.atoms, self.coords)

    def __str__(self):
        is_restricted = "unrestricted" if self.unrestricted else "restricted"
        return (
            f"Wavefunction({self.atom_num} atoms, charge={self.charge}, {is_restricted}, "
            f"mult={self.mult})"
        )

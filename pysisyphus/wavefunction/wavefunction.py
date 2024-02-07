# [1] https://pubs.acs.org/doi/pdf/10.1021/j100180a030.
#     Toward a Systematic Molecular Orbital Theory for Excited States
#     Foresman, Head-Gordon, Pople, Frisch, 1991

import contextlib
from enum import Enum
import operator
from pathlib import Path
from typing import Dict, Literal, List, Optional, Tuple
import warnings

import numpy as np
from numpy.typing import NDArray

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.elem_data import nuc_charges_for_atoms, MASS_DICT
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.helpers import BFType
from pysisyphus.wavefunction.multipole import (
    get_multipole_moment,
    get_transition_multipole_moment,
)
from pysisyphus.wavefunction.shells import Shells


Center = Literal["coc", "com"]
DensityType = Enum("DensityType", ["UNRELAXED", "RELAXED"])


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
        ecp_electrons=None,
        strict=True,
        warn_charge=4,
    ):
        self.atoms = atoms
        self.coords = np.array(coords).flatten()
        assert 3 * len(self.atoms) == len(self.coords)

        self.mult = mult
        if self.mult != 1:
            unrestricted = True
        self.unrestricted = unrestricted
        self.restricted = not self.unrestricted
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
        if ecp_electrons is None:
            ecp_electrons = np.zeros(len(self.atoms))
        elif isinstance(ecp_electrons, dict):
            _ecp_electrons = np.zeros(len(self.atoms))
            for k, v in ecp_electrons.items():
                _ecp_electrons[k] = v
            ecp_electrons = _ecp_electrons
        self.ecp_electrons = np.array(ecp_electrons, dtype=int)
        self.charge = charge
        if abs(self.charge) > warn_charge:
            warnings.warn(
                f"Encountered charge={self.charge} with high absolute value!\n"
                "This may happen in systems with ECPs. In such cases please set "
                "'ecp_electrons' or at least the correct total charge."
            )

        self._masses = np.array([MASS_DICT[atom.lower()] for atom in self.atoms])

        # Wavefunction must be initialized with the GS density
        self._current_density_key = (0, DensityType.RELAXED)
        self._densities = {
            # self._current_density_key: self.P_tot,
            self._current_density_key: self.P,
        }

        if strict:
            self.check_sanity()

    def check_sanity(self):
        assert self.shells is not None
        S = self.S  # Overlap matrix
        diag_S = np.diag(S)
        np.testing.assert_allclose(diag_S, np.ones_like(diag_S), atol=1e-12)
        Pa_mo, Pb_mo = [
            C.T @ S @ P @ S @ C for C, P in zip(self.C, self.P)
        ]  # Density matrices in MO basis
        calc_occ_a = np.trace(Pa_mo)
        calc_occ_b = np.trace(Pb_mo)
        np.testing.assert_allclose((calc_occ_a, calc_occ_b), self.occ)

        assert (
            self.ecp_electrons >= 0
        ).all(), "All entries in 'ecp_electrons' must be >= 0!"
        # We can't use self.nuc_charges, as it already contains ecp_electrons
        electrons_expected = (
            nuc_charges_for_atoms(self.atoms).sum()
            - self.charge  # Negative charge means MORE electrons, not less!
            - self.ecp_electrons.sum()
        )
        electrons_present = sum(self.occ)
        electrons_missing = electrons_expected - electrons_present
        assert electrons_missing == 0, (
            f"{electrons_missing} electrons are missing! Did you forget to specify "
            "'ecp_electrons' and/or the correct 'charge'?"
        )

    def get_permut_matrix(self):
        bf_type = self.bf_type
        if bf_type == BFType.CARTESIAN:
            return self.shells.P_cart
        elif bf_type == BFType.PURE_SPHERICAL:
            return self.shells.P_sph
        else:
            raise Exception(f"Unknown {bf_type=}!")

    @property
    def atom_num(self):
        return len(self.atoms)

    @property
    def coords3d(self):
        return self.coords.reshape(-1, 3)

    @property
    def masses(self):
        return self._masses

    @property
    def nuc_charges(self):
        return nuc_charges_for_atoms(self.atoms) - self.ecp_electrons

    @property
    def mo_num(self):
        return self.C.shape[1]

    @staticmethod
    def from_file(fn, **kwargs):
        if str(fn).startswith("lib:"):
            fn = WF_LIB_DIR / fn[4:]
        path = Path(fn)

        if not path.exists():
            raise FileNotFoundError(path)

        from_funcs = {
            ".bson": Wavefunction.from_orca_bson,
            ".fchk": Wavefunction.from_fchk,
            ".json": Wavefunction.from_orca_json,
            ".molden": Wavefunction.from_molden,
        }
        from_funcs_for_line = (
            # Molden format
            ("[Molden Format]", Wavefunction.from_molden),
            # AOMix, e.g. from Turbomole
            ("[AOMix Format", Wavefunction.from_aomix),
        )
        # If possible I would advise to stay away from .molden files :)
        try:
            from_func = from_funcs[path.suffix.lower().strip()]
        except KeyError:
            # Try to guess wavefunction kind from first line
            with open(fn) as handle:
                first_line = handle.readline().strip()
            for key, func in from_funcs_for_line:
                if first_line.startswith(key):
                    from_func = func
                    break
            else:
                raise Exception("Could not determine file format!")
        return from_func(path, **kwargs)

    @staticmethod
    @file_or_str(".molden", ".molden.input")
    def from_molden(text, **kwargs):
        from pysisyphus.io.molden import wavefunction_from_molden

        wf = wavefunction_from_molden(text, **kwargs)
        return wf

    @staticmethod
    @file_or_str(".json")
    def from_orca_json(text, **kwargs):
        """Create wavefunction from ORCA JSON.

        As of version 5.0.3 ORCA does not create JSON files for systems
        containing an ECP, so this method does not take any additional
        args or kwargs in contrast to from_orca_molden."""
        from pysisyphus.io.orca import wavefunction_from_json

        wf = wavefunction_from_json(text, **kwargs)
        return wf

    @staticmethod
    @file_or_str(".bson", mode="rb")
    def from_orca_bson(text, **kwargs):
        """Create wavefunction from ORCA BSON.

        See from_orca_json for further comments."""
        from pysisyphus.io.orca import wavefunction_from_bson

        return wavefunction_from_bson(text, **kwargs)

    @staticmethod
    @file_or_str(".molden", ".molden.input")
    def from_orca_molden(text, **kwargs):
        """Create wavefunction from ORCA molden file.

        While ORCA refuses to create JSON files for systems containing
        ECPs, this is not the case for 'orca_2mkl [fn] -molden'. So we may
        encounter ECPs here. To handle this 'wavefunction_from_molden' accepts
        an additional charge argument, to specify the correct charge, e.g. as
        stored in an ORCA calculator. If the actual, true charge is not availble
        wavefunction_from_molden will come up with an absurdly high charge.
        """

        from pysisyphus.io.orca import wavefunction_from_orca_molden

        return wavefunction_from_orca_molden(text, **kwargs)

    @staticmethod
    @file_or_str(".molden", ".molden.input")
    def from_fchk(text, **kwargs):
        from pysisyphus.io.fchk import wavefunction_from_fchk

        return wavefunction_from_fchk(text, **kwargs)

    @staticmethod
    @file_or_str(".in")
    def from_aomix(text, **kwargs):
        from pysisyphus.io.aomix import wavefunction_from_aomix

        return wavefunction_from_aomix(text, **kwargs)

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
        return np.array([C_occ @ C_occ.T for C_occ in self.C_occ])

    @property
    def P_tot(self):
        return operator.add(*self.P)

    def P_exc(self, trans_dens):
        """
        Eqs. (2.25) and (2.26) in [1].
        """
        trans_dens *= 2**0.5 / 2
        occ_a, occ_b = self.occ
        assert occ_a == occ_b
        occ = occ_a
        occupations = np.zeros(self.mo_num)
        occupations[:occ_a] = 2
        if self.unrestricted:
            raise Exception("Fix density matrix construction!")
        P = np.diag(occupations)
        dP_oo = np.einsum("ia,ja->ij", trans_dens, trans_dens)
        dP_vv = np.einsum("ia,ic->ac", trans_dens, trans_dens)
        P[:occ, :occ] -= 2 * dP_oo
        P[occ:, occ:] += 2 * dP_vv
        C, _ = self.C
        # The density matric is currently still in the MO basis. Transform
        # it to the AO basis and return.
        return C @ P @ C.T

    """
    def store_density(self, P, name, ao_or_mo="ao"):
        assert P.ndim == 2, "Handling of alpha/beta densities is not yet implemented!"
        if ao_or_mo == "mo":
            C, _ = self.C
            P_ao = C @ P @ C.T
        elif ao_or_mo == "ao":
            P_ao = P.copy()
        self._densities[name] = P_ao

    def get_density(self, name):
        return self._densities[name]
    """

    @property
    def densities(self):
        return self._densities

    def set_current_density(self, key):
        self._current_density_key = key

    def get_current_density(self):
        return self._densities[self._current_density_key]

    @contextlib.contextmanager
    def current_density(self, key):
        key_bak = self._current_density_key
        self.set_current_density(key)
        try:
            yield self
        finally:
            self.set_current_density(key_bak)

    def get_current_total_density(self):
        return self._densities[self._current_density_key].sum(axis=0)

    def get_relaxed_density(self, root: int):
        return self._densities[(root, DensityType.RELAXED)]

    def get_total_relaxed_density(self, root: int):
        return self.get_relaxed_density(root).sum(axis=0)

    def set_relaxed_density(self, root: int, density: np.ndarray):
        key = (root, DensityType.RELAXED)
        self._densities[key] = density
        return key

    def eval_density(self, coords3d, P=None):
        if P is None:
            P = self.P_tot
        spherical = self.bf_type == BFType.PURE_SPHERICAL
        vals = self.shells.eval(coords3d, spherical=spherical)
        rho = np.einsum(
            "uv,iu,iv->i", P, vals, vals, optimize=["einsum_path", (0, 1), (0, 1)]
        )
        return rho

    def eval_esp(self, coords3d, with_nuc=True):
        charges = np.ones(1)  # Assume positive charge of 1 e
        esp = np.zeros(len(coords3d))
        nuc_coords3d = self.coords3d
        Pa, Pb = self.P
        for i, c3d in enumerate(coords3d):
            V = self.get_V(c3d[None, :], charges)  # 1el Coulomb integrals
            el = np.einsum("uv,uv->", Pa, V)
            if self.unrestricted:
                el += np.einsum("uv,uv->", Pb, V)
            else:
                el *= 2.0
            if with_nuc:
                nuc = (
                    self.nuc_charges
                    / np.linalg.norm(nuc_coords3d - c3d[None, :], axis=1)
                ).sum()
            else:
                nuc = 0.0
            esp[i] = el + nuc
        return esp

    @property
    def ao_centers(self) -> List[int]:
        return list(
            {
                BFType.CARTESIAN: lambda: self.shells.cart_ao_centers,
                BFType.PURE_SPHERICAL: lambda: self.shells.sph_ao_centers,
            }[self.bf_type]()
        )

    @property
    def ao_center_map(self) -> Dict[int, List[int]]:
        ao_center_map = dict()
        for i, aoc in enumerate(self.ao_centers):
            ao_center_map.setdefault(aoc, list()).append(i)
        return ao_center_map

    @property
    def is_cartesian(self) -> bool:
        return self.bf_type == BFType.CARTESIAN

    #####################
    # Overlap Integrals #
    #####################

    @property
    def S(self):
        # Check what type of basis functions we are using.
        return {
            BFType.CARTESIAN: lambda: self.shells.S_cart,
            BFType.PURE_SPHERICAL: lambda: self.shells.S_sph,
        }[self.bf_type]()

    @property
    def S_from_C(self):
        """Reconstructed overlap-matrix from 1 = C.T @ S @ C.

        It will have the same BFType as the underlying orbitals!"""
        C, _ = self.C
        C_inv = np.linalg.pinv(C, rcond=1e-8)
        return C_inv.T @ C_inv

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

    #####################
    # Multipole Moments #
    #####################

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

    def quadrupole_ints(
        self, origin: Optional[NDArray] = None, kind: Center = "coc"
    ) -> NDArray[float]:
        if origin is None:
            origin = self.get_origin(kind=kind)

        quadrupole_ints = {
            BFType.CARTESIAN: lambda *args: self.shells.get_quadrupole_ints_cart(*args),
            BFType.PURE_SPHERICAL: lambda *args: self.shells.get_quadrupole_ints_sph(
                *args
            ),
        }[self.bf_type](origin)
        return quadrupole_ints

    def get_dipole_moment(
        self,
        P: Optional[NDArray[float]] = None,
        origin: Optional[NDArray[float]] = None,
        kind: Center = "coc",
    ) -> NDArray[float]:
        if origin is None:
            origin = self.get_origin(kind=kind)
        if P is None:
            P = self.get_current_total_density()
        dipole_ints = self.dipole_ints(origin)
        return get_multipole_moment(
            1, self.coords3d, origin, dipole_ints, self.nuc_charges, P
        )

    @property
    def dipole_moment(self) -> NDArray[float]:
        return self.get_dipole_moment()

    def get_quadrupole_moment(
        self,
        P: Optional[NDArray[float]] = None,
        origin: Optional[NDArray[float]] = None,
        kind: Center = "coc",
    ) -> NDArray[float]:
        if origin is None:
            origin = self.get_origin(kind=kind)
        if P is None:
            P = self.get_current_total_density()
        quadrupole_ints = self.quadrupole_ints(origin)
        return get_multipole_moment(
            2, self.coords3d, origin, quadrupole_ints, self.nuc_charges, P
        )

    @property
    def quadrupole_moment(self) -> NDArray[float]:
        return self.get_quadrupole_moment()

    ################################
    # Transition Multipole Moments #
    ################################

    def get_transition_multipole_moment(
        self,
        order: int,
        T_a: NDArray[float],
        T_b: NDArray[float] = None,
        origin: Optional[NDArray[float]] = None,
        kind: Center = "coc",
        full=False,
    ) -> NDArray[float]:
        if origin is None:
            origin = self.get_origin(kind=kind)
        if order == 1:
            multipole_ints = self.dipole_ints(origin)
        elif order == 2:
            multipole_ints = self.quadrupole_ints(origin)
        else:
            raise Exception("Multipoles of order {order} are not implemented!")
        C_a, C_b = self.C
        occ_a, occ_b = self.occ
        return get_transition_multipole_moment(
            multipole_ints,
            C_a,
            C_b,
            T_a,
            T_b,
            occ_a,
            occ_b,
            full=full,
        )

    def get_transition_dipole_moment(self, *args, **kwargs):
        return self.get_transition_multipole_moment(1, *args, **kwargs)

    def get_transition_quadrupole_moment(self, *args, **kwargs):
        return self.get_transition_multipole_moment(2, *args, **kwargs)

    def oscillator_strength(
        self, exc_ens: NDArray[float], trans_moms: NDArray[float]
    ) -> NDArray[float]:
        """Oscillator strength from TDMs and excitation energies."""
        exc_ens = np.atleast_1d(exc_ens)
        if trans_moms.ndim == 1:
            trans_moms = trans_moms[None, :]
        fosc = 2 / 3 * exc_ens * (trans_moms**2).sum(axis=1)
        return fosc

    #########################
    # 1el Coulomb Integrals #
    #########################

    def get_V(self, coords3d, charges):
        return {
            BFType.CARTESIAN: self.shells.get_V_cart,
            BFType.PURE_SPHERICAL: self.shells.get_V_sph,
        }[self.bf_type](coords3d, charges)

    def __str__(self):
        is_restricted = "unrestricted" if self.unrestricted else "restricted"
        return (
            f"Wavefunction({self.atom_num} atoms, charge={self.charge}, {is_restricted}, "
            f"mult={self.mult})"
        )

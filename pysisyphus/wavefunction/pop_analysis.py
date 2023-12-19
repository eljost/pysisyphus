from dataclasses import dataclass
import functools
from typing import List, Tuple

import numpy as np

import scipy as sp

from pysisyphus.wavefunction.helpers import symmetric_orthogonalization
from pysisyphus.wavefunction import Wavefunction


@dataclass
class PopAnalysis:
    # Also store atoms and coordinates?!
    # atoms: List[str]
    # coords3d: np.ndarray

    # Population of alpha and beta electrons
    pop_a: np.ndarray
    pop_b: np.ndarray
    nuc_charges: np.ndarray

    @property
    def charges(self):
        return self.nuc_charges - self.pop_a - self.pop_b

    @property
    def tot_charge(self):
        return self.charges.sum()

    @property
    def spin_pop(self):
        return np.abs(self.pop_a - self.pop_b)

    @property
    def alpha_spin_pop(self):
        spin_pop = self.pop_a - self.pop_b
        beta_mask = spin_pop < 0.0
        spin_pop[beta_mask] = 0.0
        return spin_pop

    @property
    def beta_spin_pop(self):
        spin_pop = self.pop_a - self.pop_b
        alpha_mask = spin_pop > 0.0
        spin_pop[alpha_mask] = 0.0
        return np.abs(spin_pop)


@functools.singledispatch
def mulliken_charges(
    P: Tuple[np.ndarray, np.ndarray],
    S: np.ndarray,
    nuc_charges: np.ndarray,
    ao_centers: List[int],
) -> PopAnalysis:
    def mulliken_atom_pops(P: np.ndarray, S: np.ndarray) -> np.ndarray:
        ao_populations = np.einsum("ij,ji->i", P, S)
        atom_populations = np.zeros(len(nuc_charges))
        for i, center in enumerate(ao_centers):
            atom_populations[center] += ao_populations[i]
        return atom_populations

    atom_populations_a = mulliken_atom_pops(P[0], S)
    atom_populations_b = mulliken_atom_pops(P[1], S)

    pop_ana = PopAnalysis(
        pop_a=atom_populations_a,
        pop_b=atom_populations_b,
        nuc_charges=nuc_charges,
    )
    return pop_ana


@mulliken_charges.register
def _(wf: Wavefunction) -> PopAnalysis:
    return mulliken_charges(
        wf.P,
        wf.S,
        wf.nuc_charges,
        wf.ao_centers,
    )


def make_iaos(
    C_occ: np.ndarray,
    S_org: np.ndarray,
    S_minao: np.ndarray,
    S_cross: np.ndarray,
) -> np.ndarray:
    """Intrinsic atomic orbitals.

    [1] https://doi.org/10.1021/ct400687b
    """
    P12 = sp.linalg.solve(S_org, S_cross)  # Projector on basis1, S1⁻¹ S12
    P21 = sp.linalg.solve(S_minao, S_cross.T)  # Projector on basis2, S2⁻¹ S21

    # Depolarized MOs in basis2
    C2_depol = P12.dot(P21).dot(C_occ)
    C2_depol = symmetric_orthogonalization(C2_depol, S_org)

    O = C_occ.dot(C_occ.T).dot(S_org)
    O_ = C2_depol.dot(C2_depol.T).dot(S_org)
    iaos = (2 * O.dot(O_) - O_ - O).dot(P12) + P12
    iaos = symmetric_orthogonalization(iaos, S_org)

    return iaos


@functools.singledispatch
def iao_charges(
    P: Tuple[np.ndarray, np.ndarray],
    nuc_charges: np.ndarray,
    ao_centers: List[int],
) -> PopAnalysis:
    return mulliken_charges(
        P,
        # IAOs are are orthonormal, so the overlap matriy is the identity matrix
        np.eye(len(ao_centers)),
        nuc_charges,
        ao_centers,
    )


@iao_charges.register
def _(wf: Wavefunction) -> PopAnalysis:
    """IAO charges.

    Extension to ES is described here:
        https://chemistry.stackexchange.com/a/75913
    """
    S_org = wf.S
    C_occ = wf.C_occ
    minao_shells = wf.shells.from_basis("minao")
    S_minao = getattr(minao_shells, "S_cart" if wf.is_cartesian else "S_sph")
    S_cross = wf.S_with_shells(
        minao_shells
    )  # Overlap between original and MINAO  basis

    def get_iao_P(C_occ: np.ndarray):
        iaos = make_iaos(C_occ, S_org, S_minao, S_cross)
        C_iao = iaos.T @ S_org @ C_occ  # Projection of original MOs onto IAOs
        P_iao = C_iao @ C_iao.T
        return P_iao

    # Always assume separate α C_occ and β C_occ matrices.
    P_iao = np.array([get_iao_P(c_occ) for c_occ in C_occ])

    minao_ao_centers = list(minao_shells.sph_ao_centers)
    return iao_charges(
        P_iao,
        wf.nuc_charges,
        minao_ao_centers,
    )

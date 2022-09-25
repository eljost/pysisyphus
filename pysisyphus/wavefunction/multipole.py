# [1] https://doi.org/10.1063/1.4937410
#     The consequences of improperly describing oscillator strengths
#     beyond the electric dipole approximation
#     Lestrange, Egidi, Li, 2015
# [2] https://doi.org/10.1016/0009-2614(95)01036-9
#     Ab initio calculation and display of the rotary strength tensor in
#     the random phase approximation. Method and model studies.
#     Pedersen, Hansen, 1995

import numpy as np
from numpy.typing import NDArray


def get_multipole_moment(
    order: int,
    coords3d: NDArray[float],
    origin: NDArray[float],
    multipole_ints: NDArray[float],
    nuc_charges: NDArray[int],
    P: NDArray[float],
) -> NDArray[float]:
    """
    order
        Kind of requested multipoles: (1, dipole moment), (2, quadrupole moment).
    coords3d
        Cartesian coordinates of shape (natoms, 3).
    origin
        Origin of the multipole expansion. The supplied integrales must
        be calculated w.r.t. this origin.
    multipole_ints
        Linear or quadratic multipole moment integrals.
    nuc_charges
        Nuclear charges.
    P
        Density matrix in the AO basis.
    """
    origin = np.array(origin, dtype=float)
    assert origin.size == 3
    assert len(coords3d) == len(nuc_charges)

    # Make a copy, as the coordinates must be give w.r.t. the origin, which
    # may be != (0., 0., 0.).
    coords3d = coords3d.copy()
    coords3d -= origin[None, :]

    electronic = np.einsum("...ij,ji->...", multipole_ints, P)
    if order == 1:
        nuclear = np.einsum("i,ix->x", nuc_charges, coords3d)
    elif order == 2:
        nuclear = np.einsum("i,ix,iy->xy", nuc_charges, coords3d, coords3d)
    else:
        raise Exception(f"Multipoles of order={order} are not implemented!")
    moment = nuclear - electronic
    return moment


def get_transition_multipole_moment(
    multipole_ints: NDArray[float],
    C_a: NDArray[float],
    C_b: NDArray[float],
    occ_a: int,
    occ_b: int,
    P_a: NDArray[float],
    P_b: NDArray[float] = None,
) -> NDArray[float]:
    """
    Transition multipole moments from transition density matrices.

    See Eq. (33) in [1].

    Parameters
    ----------
    multipole_integrals
        Multipole integrals in the AO basis of the desired ordered,
        e.g., linear moment for transition dipole moments or quadratic
        moments for quadrupole transition moments.
    C_a
        MO-coefficients for α-electrons.
    C_b
        MO-coefficients for β-electrons.
    occ_a
        Number of α-electrons.
    occ_b
        Number of β-electrons.
    P_a
        Numpy array containing transition density matrices of shape
        (nstates, occ, virt) or (occ, virt) in the α-electron space.
        The transition density must be given in the MO-basis!
    P_b
        See P_a. Optional. If not provided, P_b is derived from P_a.

    Returns
    -------
    trans_moments
        Transition moments of the respective order for all provided
        states.
    """

    # Deal with unrestricted transition density matrices that contain
    # only one state.
    def atleast_3d(tden):
        if tden.ndim == 2:
            tden = tden[None, :]
        return tden

    P_a = atleast_3d(P_a)
    if P_b is None:
        P_a = 1 / 2**0.5 * P_a
        P_b = P_a
    P_b = atleast_3d(P_b)

    # Expected shape: (nstates, occ, virt)
    assert P_a.ndim == P_b.ndim == 3

    def get_trans_moment(
        C: NDArray[float], occ: NDArray[int], trans_dens: NDArray[float]
    ) -> NDArray[float]:
        # Transform AO integrals to MO basis
        multipole_ints_mo = np.einsum(
            "jl,...lm,mk->...jk", C[:, :occ].T, multipole_ints, C[:, occ:]
        )

        """
        Then contract with multipole integrals. For TDMs see Eq. (18) in [2].
        i : Excited state number
        j : occ. MO space
        k : virt. MO space
        """
        trans_moment = np.einsum("...jk,ijk->i...", multipole_ints_mo, trans_dens)
        return trans_moment

    # Transitions between α -> α and β -> β
    trans_moment_a = get_trans_moment(C_a, occ_a, P_a)
    trans_moment_b = get_trans_moment(C_b, occ_b, P_b)
    trans_moment = trans_moment_a + trans_moment_b
    return trans_moment

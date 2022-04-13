import itertools as it
from typing import Literal, get_args

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import dok_matrix

from pysisyphus.helpers_pure import log
from pysisyphus.intcoords.setup import get_pair_covalent_radii
from pysisyphus.intcoords.setup_fast import (
    find_bonds_for_geom,
    find_bonds_bends,
    find_bonds_bends_dihedrals,
)
from pysisyphus.intcoords.derivatives import dq_b, dq_a, dq_d
from pysisyphus.intcoords import RedundantCoords
from pysisyphus.optimizers.guess_hessians import get_lindh_alpha


# [1] https://www.nature.com/articles/s41598-018-32105-x
#     Mones, 2018


def get_lindh_k(atoms, coords3d, bonds=None, angles=None, torsions=None):
    if bonds is None:
        bonds = list()
    if angles is None:
        angles = list()
    if torsions is None:
        torsions = list()

    atoms = [a.lower() for a in atoms]

    alphas = [get_lindh_alpha(a1, a2) for a1, a2 in it.combinations(atoms, 2)]
    pair_cov_radii = get_pair_covalent_radii(atoms)
    cdm = pdist(coords3d)
    rhos = squareform(np.exp(alphas * (pair_cov_radii ** 2 - cdm ** 2)))

    k_dict = {
        2: 0.45,  # Stretches/bonds
        3: 0.15,  # Bends/angles
        4: 0.005,  # Torsions/dihedrals
    }
    ks = list()
    for inds in it.chain(bonds, angles, torsions):
        rho_product = 1
        for i in range(len(inds) - 1):
            i1, i2 = inds[i : i + 2]
            rho_product *= rhos[i1, i2]
        ks.append(k_dict[len(inds)] * rho_product)
    return ks


def get_lindh_precon(
    atoms,
    coords,
    bonds=None,
    bends=None,
    dihedrals=None,
    c_stab=0.0103,
    logger=None,
):
    """c_stab = 0.00103 hartree/bohr² corresponds to 0.1 eV/Å² as
    given in the paper."""

    if bonds is None:
        bonds = list()
    if bends is None:
        bends = list()
    if dihedrals is None:
        dihedrals = list()

    dim = coords.size
    c3d = coords.reshape(-1, 3)

    # Calculate Lindh force constants
    ks = get_lindh_k(atoms, c3d, bonds, bends)

    grad_funcs = {
        2: dq_b,  # Bond
        3: dq_a,  # Bend
        4: dq_d,  # Dihedral
    }

    P = dok_matrix((dim, dim))
    for inds, k in zip(it.chain(bonds, bends, dihedrals), ks):
        # First derivatives of internal coordinates w.r.t cartesian coordinates
        int_grad = grad_funcs[len(inds)](*c3d[inds].flatten())
        # Assign to the correct cartesian indices
        cart_inds = np.array(list(it.chain(*[range(3 * i, 3 * i + 3) for i in inds])))
        P[cart_inds[:, None], cart_inds[None, :]] += abs(k) * np.outer(
            int_grad, int_grad
        )

    # Add stabilization to diagonal
    P[np.diag_indices(dim)] += c_stab
    P = P.tocsc()
    filled = P.size / dim ** 2
    log(logger, f"Preconditioner P has {P.size} entries ({filled:.2%} filled)")

    return P


PreconKind = Literal["full", "full_fast", "bonds", "bonds_bends"]


def precon_getter(geom, c_stab=0.0103, kind="full", logger=None):
    valid_kinds = get_args(PreconKind)
    assert kind in valid_kinds, f"Invalid kind='{kind}'! Valid kinds are: {valid_kinds}"

    atoms = geom.moving_atoms
    if len(geom.freeze_atoms) > 0:
        assert (
            kind == "full"
        ), "Preconditioning with frozen atoms is only supported for kind='full'"
    # Default empty lists for coordinates that may be skipped
    # for kind != "full".
    bends = list()
    dihedrals = list()
    if kind == "full":
        internal = RedundantCoords(
            atoms,
            geom.coords,
        )
        bonds = internal.bond_atom_indices
        bends = internal.bend_atom_indices
        dihedrals = internal.dihedral_atom_indices
    elif kind == "full_fast":
        bonds, bends, dihedrals = find_bonds_bends_dihedrals(geom)
    elif kind == "bonds_bends":
        bonds, bends = find_bonds_bends(geom)
    elif kind == "bonds":
        bonds = find_bonds_for_geom(geom)

    msg = (
        f"Constructing preconditioner from {len(bonds)} bonds, {len(bends)} bends "
        f"and {len(dihedrals)} dihedrals (kind='{kind}')."
    )
    log(logger, msg)

    def wrapper(coords):
        P = get_lindh_precon(
            atoms,
            coords,
            bonds,
            bends,
            dihedrals,
            c_stab=c_stab,
            logger=logger,
        )
        return P

    return wrapper

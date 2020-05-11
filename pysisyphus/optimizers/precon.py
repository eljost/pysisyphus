import itertools as it

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csc_matrix

from pysisyphus.intcoords.findbonds import get_pair_covalent_radii
from pysisyphus.intcoords.derivatives import dq_b, dq_a, dq_d
from pysisyphus.InternalCoordinates import RedundantCoords


# [1] https://www.nature.com/articles/s41598-018-32105-x
#     Mones, 2018


def get_lindh_alpha(atom1, atom2):
    first_period = ("h", "he")

    if (atom1 in first_period) and (atom2 in first_period):
        return 1.
    elif (atom1 in first_period) or (atom2 in first_period):
        return 0.3949
    else:
        return 0.28


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
    rhos = squareform(np.exp(alphas*(pair_cov_radii**2-cdm**2)))

    k_dict = {
        2: 0.45,  # Stretches/bonds
        3: 0.15,  # Bends/angles
        4: 0.005, # Torsions/dihedrals
    }
    ks = list()
    for inds in it.chain(bonds, angles, torsions):
        rho_product = 1
        for i in range(inds.size-1):
            i1, i2 = inds[i:i+2]
            rho_product *= rhos[i1, i2]
        ks.append(k_dict[len(inds)] * rho_product)
    return ks


def get_lindh_precon(atoms, coords, bonds=None, bends=None, dihedrals=None,
                     c_stab=0.0103):
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
        # Bond
        2: lambda i, j: dq_b(*c3d[i], *c3d[j]),
        # Bend
        3: lambda i, j, k: dq_a(*c3d[i], *c3d[j], *c3d[k]),
        # Dihedral
        4: lambda i, j, k, l: dq_d(*c3d[i], *c3d[j], *c3d[k]),
    }

    row = np.zeros(dim)
    P = np.zeros((dim, dim))
    for inds, k in zip(it.chain(bonds, bends, dihedrals), ks):
        # Construct full row
        cart_inds = list(it.chain(*[range(3*i,3*i+3) for i in inds]))

        # First derivatives of internal coordinates w.r.t cartesian coordinates
        int_grad = grad_funcs[len(inds)](*inds)
        # Assign to the correct cartesian indices
        full_row = row.copy()
        full_row[cart_inds] = int_grad
        P += np.outer(full_row, full_row*abs(k))

    # Add stabilization to diagonal
    P += c_stab*np.eye(dim)
    # Convert to sparse matrix
    P = csc_matrix(P)

    return P


def precon_getter(geometry, c_stab=0.0103):
    atoms = geometry.atoms
    internal = RedundantCoords(atoms, geometry.cart_coords)
    bonds  = internal.bond_indices
    bends = internal.bending_indices
    dihedrals = internal.dihedrals

    def wrapper(coords):
        P = get_lindh_precon(
                atoms, coords,
                bonds, bends, dihedrals,
                c_stab=c_stab,
        )
        return P
    return wrapper

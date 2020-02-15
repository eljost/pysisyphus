import itertools as it

import jax
import jax.numpy as jnp
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csc_matrix

from pysisyphus.intcoords.findbonds import get_pair_covalent_radii
from pysisyphus.intcoords.derivatives import dq_b, dq_a, dq_d, d2q_b, d2q_a


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


class Bond:

    def __init__(self, i, j, k, r_eq):
        self.i = i
        self.j = j
        self.k = k
        self.r_eq = r_eq


def bond_energy(coords3d, bond):
    r = jnp.linalg.norm(coords3d[bond.i] - coords3d[bond.j])

    energy = 0.5 * bond.k * (r - bond.r_eq)**2

    return energy


_bond_gradient = jax.grad(bond_energy)
def bond_gradient(coords3d, bond):
    return _bond_gradient(coords3d, bond).flatten()


_bond_hessian = jax.hessian(bond_energy)
def bond_hessian(coords3d, bond):
    return _bond_hessian(coords3d, bond).reshape(-1, coords3d.size)


def get_precon(atoms, coords, bonds=None, angles=None, torsions=None, c_stab=0.00103):
    """c_stab = 0.00103 hartree/bohr² corresponds to 0.1 eV/Å²"""

    if bonds is None:
        bonds = list()
    if angles is None:
        angles = list()
    if torsions is None:
        torsions = list()

    dim = coords.size
    c3d = coords.reshape(-1, 3)

    ks = get_lindh_k(atoms, c3d, bonds, angles)

    hess_funcs = {
        2: lambda i, j: d2q_b(*c3d[i], *c3d[j]),
        3: lambda i, j, k: d2q_a(*c3d[i], *c3d[j], *c3d[k]),
    }

    H_flat = np.zeros(dim*dim)
    P = np.zeros((dim, dim))
    for inds, k in zip(it.chain(bonds, angles), ks):
        int_hess = hess_funcs[len(inds)](*inds)
        # Construct full row
        H = H_flat.copy()
        cart_inds = list(it.chain(*[range(3*i,3*i+3) for i in inds]))
        flat_inds = [row*dim + col for row, col in it.product(cart_inds, cart_inds)]
        H[flat_inds] = int_hess
        # P += np.outer(full_row, full_row*abs(k))
        # P += H.reshape(dim, dim) * abs(k)
        i, j = inds
        req = np.linalg.norm(c3d[i]-c3d[j])
        b = Bond(i, j, abs(k), req)
        bh = bond_hessian(c3d, b)
        P += bh

    P += c_stab*np.eye(dim)
    P = csc_matrix(P)
    return P

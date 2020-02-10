#!/usr/bin/env python3

# [1] https://www.sciencedirect.com/science/article/pii/000926149500646L
#     Lindh, 1995
# [2] https://pubs.acs.org/doi/pdf/10.1021/j100203a036
#     Fischer, Almlöf, 1992
# [3] https://onlinelibrary.wiley.com/doi/full/10.1002/qua.21049
#     Swart, Bickelhaupt, 2006

from math import exp
import itertools as it

import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.calculators.XTB import XTB
from pysisyphus.intcoords.findbonds import get_pair_covalent_radii, get_bond_mat


def fischer_guess(geom):
    cdm = pdist(geom.coords3d)
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)

    bonds = geom.internal.bonds
    dihedrals = geom.internal.dihedrals
    # For the dihedral force constants we also have to count the number
    # of bonds formed with the centrals atoms of the dihedral.
    central_atoms = [dh.inds[1:3] for dh in dihedrals]
    bond_factor = geom.internal.bond_factor
    bond_mat = squareform(cdm <= (pair_cov_radii * geom.internal.bond_factor))
    bm = get_bond_mat(geom, geom.internal.bond_factor)
    np.testing.assert_allclose(bm, bond_mat)
    tors_atom_bonds = dict()
    for a, b in central_atoms:
        # Substract 2, as we don't want the bond between a and b,
        # but this bond will be in both rows of the bond_mat.
        bond_sum = bond_mat[a].sum() + bond_mat[b].sum() - 2
        tors_atom_bonds[(a, b)] = bond_sum

    dist_mat = squareform(cdm)
    pair_cov_radii_mat = squareform(pair_cov_radii)

    def h_bond(bond):
        a, b = bond.inds
        r_ab = dist_mat[a, b]
        r_ab_cov = pair_cov_radii_mat[a, b]
        return 0.3601 * exp(-1.944*(r_ab - r_ab_cov))

    def h_bend(bend):
        b, a, c = bend.inds
        r_ab = dist_mat[a, b]
        r_ac = dist_mat[a, c]
        r_ab_cov = pair_cov_radii_mat[a, b]
        r_ac_cov = pair_cov_radii_mat[a, c]
        return (0.089 + 0.11/(r_ab_cov*r_ac_cov)**(-0.42)
                * exp(-0.44*(r_ab + r_ac - r_ab_cov - r_ac_cov))
        )

    def h_dihedral(dihedral):
        c, a, b, d = dihedral.inds
        r_ab = dist_mat[a, b]
        r_ab_cov = pair_cov_radii_mat[a, b]
        bond_sum = tors_atom_bonds[(a, b)]
        return (0.0015 + 14.0*bond_sum**0.57 / (r_ab*r_ab_cov)**4.0
                * exp(-2.85*(r_ab - r_ab_cov))
        )
    h_funcs = {
        2: h_bond,
        3: h_bend,
        4: h_dihedral,
    }

    h_diag = list()
    for primitive in geom.internal._prim_internals:
        f = h_funcs[len(primitive.inds)]
        h_diag.append(f(primitive))
    H = np.array(np.diagflat(h_diag))
    return H


def lindh_guess(geom):
    """Slightly modified Lindh model hessian as described in [1].

    Instead of using the tabulated r_ref,ij values from [1] we will use the
    'true' covalent radii as pyberny. The tabulated r_ref,ij value for two
    carbons (2nd period) is 2.87 Bohr. Carbons covalent radius is ~ 1.44 Bohr,
    so 2*1.44 Bohr = 2.88 Bohr which fits nicely with the tabulate value.
    If values for elements > 3rd are requested the alpha values for the 3rd
    period will be (re)used.
    """
    first_period = "h he".split()
    def get_alpha(atom1, atom2):
        if (atom1 in first_period) and (atom2 in first_period):
            return 1.
        elif (atom1 in first_period) or (atom2 in first_period):
            return 0.3949
        else:
            return 0.28
    atoms = [a.lower() for a in geom.atoms]
    alphas = [get_alpha(a1, a2)
              for a1, a2 in it.combinations(atoms, 2)]
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)
    cdm = pdist(geom.coords3d)
    rhos = squareform(np.exp(alphas*(pair_cov_radii**2-cdm**2)))

    k_dict = {
        2: 0.45,  # Stretches/bonds
        3: 0.15,  # Bends/angles
        4: 0.005,  # Torsions/dihedrals
    }
    k_diag = list()
    for primitive in geom.internal._prim_internals:
        rho_product = 1
        for i in range(primitive.inds.size-1):
            i1, i2 = primitive.inds[i:i+2]
            rho_product *= rhos[i1, i2]
        k_diag.append(k_dict[len(primitive.inds)] * rho_product)
    H = np.diagflat(k_diag)
    return H


def simple_guess(geom):
    h_dict = {
        2: 0.5,  # Stretches/bonds
        3: 0.2,  # Bends/angles
        4: 0.1,  # Torsions/dihedrals
    }
    h_diag = [h_dict[len(prim.inds)] for prim in geom.internal._prim_internals]
    return np.diagflat(h_diag)


def swart_guess(geom):
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)
    cdm = pdist(geom.coords3d)
    rhos = squareform(np.exp(-cdm/pair_cov_radii + 1))
    k_dict = {
        2: 0.35,
        3: 0.15,
        4: 0.005,
    }
    k_diag = list()
    for primitive in geom.internal._prim_internals:
        rho_product = 1
        for i in range(primitive.inds.size-1):
            i1, i2 = primitive.inds[i:i+2]
            rho_product *= rhos[i1, i2]
        k_diag.append(k_dict[len(primitive.inds)] * rho_product)
    return np.diagflat(k_diag)


def xtb_hessian(geom):
    calc = geom.calculator
    xtb_kwargs = {
        "charge": calc.charge,
        "mult": calc.mult,
        "pal": calc.pal
    }
    xtb_calc = XTB(**xtb_kwargs)
    geom_ = geom.copy()
    geom_.set_calculator(xtb_calc)
    return geom_.hessian


def ts_hessian(hessian, coord_inds, damp=0.25):
    """According to [3]"""

    inds = list(coord_inds)

    # Use a copy as diag returns only a read-only view
    diag = np.diag(hessian).copy()

    # Reverse sign of reaction coordinates and damp them
    diag[inds] = -1 * damp * diag[inds]
    ts_hess = np.diag(diag)

    # Set off-diagonal elements
    for i, j in it.combinations(inds, 2):
        # fi = force_constants[i]
        # fj = force_constants[j]
        fi = diag[i]
        fj = diag[j]
        f = -(2*fi*fj)**0.5
        ts_hess[i,j] = f
        ts_hess[j,i] = f
    return ts_hess

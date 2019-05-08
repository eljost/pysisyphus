#!/usr/bin/env python3

# [1] https://www.sciencedirect.com/science/article/pii/000926149500646L
#     Lindh, 1995

import itertools as it

import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.elem_data import COVALENT_RADII


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
    cov_radii = np.array([COVALENT_RADII[a.lower()] for a in atoms])
    rref = np.array([r1+r2 for r1, r2 in it.combinations(cov_radii, 2)])
    cdm = pdist(geom.coords3d)
    diff = rref**2 - cdm**2
    rhos = squareform(np.exp(alphas*diff))

    k_dict = {
        2: 0.45,  # Stretches/bonds
        3: 0.15,  # Bends/angles
        4: 0.005, # Torsions/dihedrals
    }
    k_diag = list()
    for primitive in geom.internal._prim_coords:
        rho_product = 1
        for i in range(primitive.inds.size-1):
            i1, i2 = primitive.inds[i:i+2]
            rho_product *= rhos[i1, i2]
        k_diag.append(k_dict[len(primitive.inds)] * rho_product)
    H = np.diagflat(k_diag)
    return H

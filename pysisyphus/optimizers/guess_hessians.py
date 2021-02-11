# [1] https://www.sciencedirect.com/science/article/pii/000926149500646L
#     Lindh, 1995
# [2] https://pubs.acs.org/doi/pdf/10.1021/j100203a036
#     Fischer, Almlöf, 1992
# [3] https://onlinelibrary.wiley.com/doi/full/10.1002/qua.21049
#     Swart, Bickelhaupt, 2006
# [4] http://dx.doi.org/10.1063/1.4952956
#     Lee-Ping Wang 2016

from math import exp
import itertools as it

import h5py
import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.calculators.XTB import XTB
from pysisyphus.intcoords.setup import get_pair_covalent_radii
from pysisyphus.io.hessian import save_hessian


CART_F = 0.05


def fischer_guess(geom):
    cdm = pdist(geom.coords3d)
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)

    dihedrals = geom.internal.dihedrals
    # For the dihedral force constants we also have to count the number
    # of bonds formed with the centrals atoms of the dihedral.
    central_atoms = [dh.inds[1:3] for dh in dihedrals]
    bond_factor = geom.internal.bond_factor
    bond_mat = squareform(cdm <= (pair_cov_radii * bond_factor))
    tors_atom_bonds = dict()
    for a, b in central_atoms:
        # Substract 2, as we don't want the bond between a and b,
        # but this bond will be in both rows of the bond_mat.
        bond_sum = bond_mat[a].sum() + bond_mat[b].sum() - 2
        tors_atom_bonds[(a, b)] = bond_sum

    dist_mat = squareform(cdm)
    pair_cov_radii_mat = squareform(pair_cov_radii)

    def h_bond(bond):
        a, b = bond.indices
        r_ab = dist_mat[a, b]
        r_ab_cov = pair_cov_radii_mat[a, b]
        return 0.3601 * exp(-1.944*(r_ab - r_ab_cov))

    def h_bend(bend):
        b, a, c = bend.indices
        r_ab = dist_mat[a, b]
        r_ac = dist_mat[a, c]
        r_ab_cov = pair_cov_radii_mat[a, b]
        r_ac_cov = pair_cov_radii_mat[a, c]
        return (0.089 + 0.11/(r_ab_cov*r_ac_cov)**(-0.42)
                * exp(-0.44*(r_ab + r_ac - r_ab_cov - r_ac_cov))
        )

    def h_dihedral(dihedral):
        # import pdb; pdb.set_trace()
        c, a, b, d = dihedral.indices
        r_ab = dist_mat[a, b]
        r_ab_cov = pair_cov_radii_mat[a, b]
        bond_sum = max(tors_atom_bonds[(a, b)], 0)
        return (0.0015 + 14.0*bond_sum**0.57 / (r_ab*r_ab_cov)**4.0
                * exp(-2.85*(r_ab - r_ab_cov))
        )
    h_funcs = {
        1: lambda _: CART_F,  # See [4], last sentences of III.
        2: h_bond,
        3: h_bend,
        4: h_dihedral,
    }

    h_diag = list()
    for primitive in geom.internal.primitives:
        f = h_funcs[len(primitive.indices)]
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
    for prim_int in geom.internal.prim_internals:
        if len(prim_int.inds) == 1:
            k_diag.append(CART_F)  # See [4], last sentences of III.
            continue

        rho_product = 1
        for i in range(len(prim_int.inds)-1):
            i1, i2 = prim_int.inds[i:i+2]
            rho_product *= rhos[i1, i2]
        k_diag.append(k_dict[len(prim_int.inds)] * rho_product)
    H = np.diagflat(k_diag)
    return H


def simple_guess(geom):
    h_dict = {
        1: CART_F,  # See [4], last sentences of III.
        2: 0.5,  # Stretches/bonds
        3: 0.2,  # Bends/angles
        4: 0.1,  # Torsions/dihedrals
    }
    h_diag = [h_dict[len(prim.indices)] for prim in geom.internal.primitives]
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
    for primitive in geom.internal.primitives:
        if len(primitive.indices) == 1:
            k_diag.append(CART_F)  # See [4], last sentence of III.
            continue

        rho_product = 1
        for i in range(len(primitive.indices)-1):
            i1, i2 = primitive.indices[i:i+2]
            rho_product *= rhos[i1, i2]
        k_diag.append(k_dict[len(primitive.indices)] * rho_product)
    return np.diagflat(k_diag)


def xtb_hessian(geom, gfn=None):
    calc = geom.calculator
    xtb_kwargs = {
        "charge": calc.charge,
        "mult": calc.mult,
        "pal": calc.pal
    }
    if gfn is not None:
        xtb_kwargs["gfn"] = gfn
    xtb_calc = XTB(**xtb_kwargs)
    geom_ = geom.copy()
    geom_.set_calculator(xtb_calc)
    return geom_.hessian


def get_guess_hessian(geometry, hessian_init, int_gradient=None,
                      cart_gradient=None, h5_fn=None):
    """Obtain/calculate (model) Hessian.

    For hessian_init="calc" the Hessian will be in the coord_type
    of the geometry, otherwise a Hessian in primitive internals will
    be returned.
    """
    model_hessian = hessian_init in ("fischer", "lindh", "simple", "swart")
    target_coord_type = geometry.coord_type

    # Recreate geometry with internal coordinates if needed
    if model_hessian and (geometry.coord_type == "cart"):
        geometry = geometry.copy(coord_type="redund")

    hess_funcs = {
        # Calculate true hessian
        "calc": lambda: (geometry.hessian, "calculated exact"),
        # Unit hessian
        "unit": lambda: (np.eye(geometry.coords.size), "unit"),
        # Fischer model hessian
        "fischer": lambda: (fischer_guess(geometry), "Fischer"),
        # Lindh model hessian
        "lindh": lambda: (lindh_guess(geometry), "Lindh"),
        # Simple (0.5, 0.2, 0.1) model hessian
        "simple": lambda: (simple_guess(geometry), "simple"),
        # Swart model hessian
        "swart": lambda: (swart_guess(geometry), "Swart"),
        # XTB hessian using GFN-2
        "xtb": lambda: (xtb_hessian(geometry, gfn=2), "GFN2-XTB"),
        # XTB hessian using GFN-1
        "xtb1": lambda: (xtb_hessian(geometry, gfn=1), "GFN1-XTB"),
    }
    try:
        H, hess_str = hess_funcs[hessian_init]()
        # self.log(f"Using {hess_str} Hessian.")
    except KeyError:
        # Only cartesian hessians can be loaded
        # self.log(f"Trying to load saved Hessian from '{self.hessian_init}'.")
        if str(hessian_init).endswith(".h5"):
            with h5py.File(hessian_init, "r") as handle:
                cart_hessian = handle["hessian"][:]
        else:
            cart_hessian = np.loadtxt(hessian_init)
        geometry.cart_hessian = cart_hessian
        # Use the previously set hessian in whatever coordinate system we
        # actually employ.
        H = geometry.hessian
        hess_str = "saved"

    if (h5_fn is not None) and (hessian_init == "calc"):
        save_hessian(h5_fn, geometry)

    if model_hessian and target_coord_type == "cart":
        if cart_gradient is not None:
            int_gradient = geometry.internal.transform_forces(cart_gradient)
        H = geometry.internal.backtransform_hessian(H, int_gradient=int_gradient)

    return H, hess_str


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

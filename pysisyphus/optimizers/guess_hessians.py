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
from typing import Literal, Optional

import h5py
import numpy as np
from numpy.typing import ArrayLike
from scipy.spatial.distance import pdist, squareform

from pysisyphus.calculators.XTB import XTB
from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords.PrimTypes import PrimTypes as PT, Bonds, Bends, Dihedrals
from pysisyphus.intcoords.setup import get_pair_covalent_radii
from pysisyphus.io.hessian import save_hessian


HessInit = Literal[
    "calc", "unit", "fischer", "lindh", "simple", "swart", "xtb", "xtb1", "xtbff"
]


# See [4], last sentences of III.
CART_F = 0.05
# Default force constants
DEFAULT_F = {
    PT.BOND: 0.5,
    PT.AUX_BOND: 0.1,
    PT.HYDROGEN_BOND: 0.1,
    PT.INTERFRAG_BOND: 0.1,
    PT.AUX_INTERFRAG_BOND: 0.05,
    PT.BEND: 0.2,
    PT.LINEAR_BEND: 0.1,
    PT.LINEAR_BEND_COMPLEMENT: 0.1,
    PT.PROPER_DIHEDRAL: 0.1,
    PT.IMPROPER_DIHEDRAL: 0.1,
    PT.OUT_OF_PLANE: 0.1,
    PT.LINEAR_DISPLACEMENT: 0.1,
    PT.LINEAR_DISPLACEMENT_COMPLEMENT: 0.1,
    PT.TRANSLATION_X: CART_F,
    PT.TRANSLATION_Y: CART_F,
    PT.TRANSLATION_Z: CART_F,
    PT.ROTATION_A: CART_F,
    PT.ROTATION_B: CART_F,
    PT.ROTATION_C: CART_F,
    PT.CARTESIAN_X: CART_F,
    PT.CARTESIAN_Y: CART_F,
    PT.CARTESIAN_Z: CART_F,
    PT.BONDED_FRAGMENT: CART_F,
    PT.DUMMY_TORSION: 0.1,
    PT.DISTANCE_FUNCTION: 0.1,
}


def simple_guess(geom):
    """Default force constants."""
    h_diag = [DEFAULT_F[type_] for type_, *_ in geom.internal.typed_prims]
    return np.diag(h_diag)


def improved_guess(geom, bond_func, bend_func, dihedral_func):
    H_guess = simple_guess(geom)
    for i, (pt, *indices) in enumerate(geom.internal.typed_prims):
        if pt in Bonds:
            f_func = bond_func
        elif pt in Bends:
            f_func = bend_func
        elif pt in Dihedrals:
            f_func = dihedral_func
        else:
            continue
        new_f = f_func(indices)
        H_guess[i, i] = new_f
    return H_guess


def fischer_guess(geom):
    cdm = pdist(geom.coords3d)
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)

    # For the dihedral force constants we also have to count the number
    # of bonds formed with the centrals atoms of the dihedral.
    central_atoms = [inds[1:3] for inds in geom.internal.dihedral_atom_indices]
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

    def h_bond(indices):
        a, b = indices[:2]
        r_ab = dist_mat[a, b]
        r_ab_cov = pair_cov_radii_mat[a, b]
        return 0.3601 * exp(-1.944 * (r_ab - r_ab_cov))

    def h_bend(indices):
        b, a, c = indices
        r_ab = dist_mat[a, b]
        r_ac = dist_mat[a, c]
        r_ab_cov = pair_cov_radii_mat[a, b]
        r_ac_cov = pair_cov_radii_mat[a, c]
        return 0.089 + 0.11 / (r_ab_cov * r_ac_cov) ** (-0.42) * exp(
            -0.44 * (r_ab + r_ac - r_ab_cov - r_ac_cov)
        )

    def h_dihedral(indices):
        c, a, b, d = indices
        r_ab = dist_mat[a, b]
        r_ab_cov = pair_cov_radii_mat[a, b]
        bond_sum = max(tors_atom_bonds[(a, b)], 0)
        return 0.0015 + 14.0 * bond_sum ** 0.57 / (r_ab * r_ab_cov) ** 4.0 * exp(
            -2.85 * (r_ab - r_ab_cov)
        )

    H = improved_guess(
        geom, bond_func=h_bond, bend_func=h_bend, dihedral_func=h_dihedral
    )
    return H


def lindh_style_guess(geom, ks, rhos):
    """Approximate force constants according to Lindh.[1]

    Bonds:       k_ij = k_r * rho_ij
    Bends:      k_ijk = k_b * rho_ij * rho_jk
    Dihedrals: k_ijkl = k_d * rho_ij * rho_jk * rho_kl
    """

    def k_func(indices):
        rho_product = 1
        inds_len = len(indices)
        for i, ind in enumerate(indices[:-1], 1):
            i1, i2 = ind, indices[i]
            rho_product *= rhos[i1, i2]
        k = ks[inds_len] * rho_product
        return k

    H = improved_guess(geom, bond_func=k_func, bend_func=k_func, dihedral_func=k_func)
    return H


def get_lindh_alpha(atom1, atom2):
    first_period = "h", "he"
    if (atom1 in first_period) and (atom2 in first_period):
        return 1.0
    elif (atom1 in first_period) or (atom2 in first_period):
        return 0.3949
    else:
        return 0.28


def lindh_guess(geom):
    """Slightly modified Lindh model hessian as described in [1].

    Instead of using the tabulated r_ref,ij values from [1] we will use the
    'true' covalent radii as pyberny. The tabulated r_ref,ij value for two
    carbons (2nd period) is 2.87 Bohr. Carbons covalent radius is ~ 1.44 Bohr,
    so 2*1.44 Bohr = 2.88 Bohr which fits nicely with the tabulate value.
    If values for elements > 3rd are requested the alpha values for the 3rd
    period will be (re)used.
    """

    atoms = [a.lower() for a in geom.atoms]
    alphas = [get_lindh_alpha(a1, a2) for a1, a2 in it.combinations(atoms, 2)]
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)
    cdm = pdist(geom.coords3d)
    rhos = squareform(np.exp(alphas * (pair_cov_radii ** 2 - cdm ** 2)))

    ks = {
        2: 0.45,  # Stretches/bonds
        3: 0.15,  # Bends/angles
        4: 0.005,  # Torsions/dihedrals
    }
    return lindh_style_guess(geom, ks, rhos)


def swart_guess(geom):
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)
    cdm = pdist(geom.coords3d)
    rhos = squareform(np.exp(-cdm / pair_cov_radii + 1))
    ks = {
        2: 0.35,
        3: 0.15,
        4: 0.005,
    }
    return lindh_style_guess(geom, ks, rhos)


def xtb_hessian(geom, gfn=None):
    calc = geom.calculator
    xtb_kwargs = {"charge": calc.charge, "mult": calc.mult, "pal": calc.pal}
    if gfn is not None:
        xtb_kwargs["gfn"] = gfn
    xtb_calc = XTB(**xtb_kwargs)
    geom_ = geom.copy()
    geom_.set_calculator(xtb_calc)
    return geom_.hessian


def get_guess_hessian(
    geometry: Geometry,
    hessian_init: HessInit,
    int_gradient: Optional[ArrayLike] = None,
    cart_gradient: Optional[ArrayLike] = None,
    h5_fn: Optional[str] = None,
):
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
        # XTB hessian using GFN-FF
        "xtbff": lambda: (xtb_hessian(geometry, gfn="ff"), "GFN-FF"),
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
        f = -((2 * fi * fj) ** 0.5)
        ts_hess[i, j] = f
        ts_hess[j, i] = f
    return ts_hess

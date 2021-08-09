import logging

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import log
from pysisyphus.intcoords.setup import get_bond_sets
from pysisyphus.intcoords import RedundantCoords
from pysisyphus.intcoords.PrimTypes import PrimTypes


logger = logging.getLogger("internal_coords")


def augment_bonds(geom, root=0, proj=False):
    assert geom.coord_type != "cart"
    log(logger, "Trying to augment bonds.")

    hessian = geom.cart_hessian
    try:
        energy = geom.energy
    except AttributeError:
        energy = None

    func = find_missing_bonds_by_projection if proj else find_missing_strong_bonds

    missing_bonds = func(geom, hessian, root=root)

    if missing_bonds:
        aux_bond_pt = PrimTypes.AUX_BOND
        missing_aux_bonds = [(aux_bond_pt, *mbond) for mbond in missing_bonds]
        print("\t@Missing bonds:", missing_bonds)
        new_geom = Geometry(geom.atoms, geom.cart_coords,
                            coord_type=geom.coord_type,
                            coord_kwargs={"define_prims": missing_aux_bonds,},
        )
        new_geom.set_calculator(geom.calculator)
        new_geom.energy = energy
        new_geom.cart_hessian = hessian
        return new_geom
    else:
        return geom


def find_missing_strong_bonds(geom, hessian, bond_factor=1.7, thresh=0.3,
                              root=0):
    # Define only bonds
    red = RedundantCoords(geom.atoms, geom.cart_coords,
                          bond_factor=bond_factor, bonds_only=True)
    cur_bonds = set([frozenset(b) for b in geom.internal.bond_atom_indices])

    # Transform cartesian hessian to bond hessian
    bond_hess = red.transform_hessian(hessian)
    # Determine transisiton vector
    eigvals, eigvecs = np.linalg.eigh(bond_hess)
    # There are probably no bonds missing if there are no negative eigenvalues
    if sum(eigvals < 0) == 0:
        return list()

    trans_vec = eigvecs[:, root]
    # Find bonds that strongly contribute to the selected transition vector
    strong = np.abs(trans_vec) > thresh
    strong_bonds = np.array(red.bond_atom_indices)[strong]
    strong_bonds = set([frozenset(b) for b in strong_bonds])

    # Check which strong bonds are missing from the currently defiend bonds
    missing_bonds = strong_bonds - cur_bonds
    missing_bonds = [tuple(_) for _ in missing_bonds]
    return missing_bonds


def find_missing_bonds_by_projection(geom, hessian, bond_factor=2.0, bond_thresh=0.35,
                                     concerted_thresh=0.35, root=0):

    def array2set(arr):
        return set([tuple(_) for _ in arr])

    bonds_present = array2set(geom.internal.bond_atom_indices)
    eigvals, eigvecs = np.linalg.eigh(hessian)

    # There are probably no bonds missing if there are no negative eigenvalues
    if sum(eigvals < 0) == 0:
        return list()

    trans_vec = eigvecs[:, root]

    c3d = geom.coords3d
    bond_vec_empty = np.zeros_like(c3d)
    unique_bonds = array2set(get_bond_sets(geom.atoms, c3d, bond_factor=bond_factor))
    unique_bonds -= bonds_present
    unique_bonds = np.array(list(unique_bonds))

    bond_vecs = list()
    concerted_vecs = list()
    for m, k in unique_bonds:
        displ = c3d[k] - c3d[m]
        displ /= np.linalg.norm(displ)

        # Bond
        bond = bond_vec_empty.copy()
        bond[k] = displ
        bond[m] = -displ
        bond_vecs.append(bond)

        # Concerted movement
        conc = bond_vec_empty.copy()
        conc[k] = displ
        conc[m] = displ
        concerted_vecs.append(conc)

    def reshape(arr):
        return np.array(arr).reshape(-1, trans_vec.size)
    bond_vecs = reshape(bond_vecs)
    concerted_vecs = reshape(concerted_vecs)

    def overlaps(arr):
        return np.abs(arr.dot(trans_vec))
    bond_ovlps = overlaps(bond_vecs)
    concerted_ovlps = overlaps(concerted_vecs)

    unique_bonds = np.array(unique_bonds)
    missing_bonds = unique_bonds[bond_ovlps > bond_thresh]
    missing_concerted = unique_bonds[concerted_ovlps > concerted_thresh]

    missing_inds = array2set(missing_bonds) | array2set(missing_concerted)
    return missing_inds

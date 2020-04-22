import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.InternalCoordinates import RedundantCoords


def augment_bonds(geom, root=0):
    assert geom.coord_type != "cart"

    hessian = geom.cart_hessian
    energy = geom.energy

    missing_bonds = find_missing_strong_bonds(geom, hessian, root=root)
    if missing_bonds:
        print("\t@Missing bonds:", missing_bonds)
        new_geom = Geometry(geom.atoms, geom.cart_coords,
                            coord_type=geom.coord_type,
                            coord_kwargs={"define_prims": missing_bonds,},
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
    cur_bonds = set([frozenset(b) for b in geom.internal.bond_indices])

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
    strong_bonds = red.bond_indices[strong]
    strong_bonds = set([frozenset(b) for b in strong_bonds])

    # Check which strong bonds are missing from the currently defiend bonds
    missing_bonds = strong_bonds - cur_bonds
    missing_bonds = [tuple(_) for _ in missing_bonds]
    return missing_bonds

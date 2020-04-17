#!/usr/bin/env python3

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.Geometry import Geometry
from pysisyphus.xyzloader import write_geoms_to_trj
from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.interpolate.LST import LST
from pysisyphus.interpolate.IDPP import IDPP


np.set_printoptions(suppress=True, precision=4)


def test_lst():
    initial = geom_from_library("dipeptide_init.xyz")
    final = geom_from_library("dipeptide_fin.xyz")

    geoms = (initial, final)
    lst = LST(geoms, 28, align=True)
    geoms = lst.interpolate_all()
    lst.all_geoms_to_trj("lst_opt.trj")


def test_idpp():
    # initial = geom_from_library("dipeptide_init.xyz")
    # final = geom_from_library("dipeptide_fin.xyz")

    initial = geom_from_xyz_file("09_htransfer_product.xyz")
    final = geom_from_xyz_file("10_po_diss_product_xtbopt.xyz")

    geoms = (initial, final)
    idpp = IDPP(geoms, 18, align=True)
    geoms = idpp.interpolate_all()
    idpp.all_geoms_to_trj("idpp_opt.trj")


def dlc_interpolate(initial, final, between=18):
    print("initial primitives", initial.coords.size)
    print("final primitives", final.coords.size)
    def to_set(iterable):
        return set([tuple(_) for _ in iterable])

    def get_ind_sets(geom):
        bonds, bends, dihedrals = geom.internal.prim_indices
        return to_set(bonds), to_set(bends), to_set(dihedrals)

    def merge_coordinate_definitions(geom1, geom2):
        bonds1, bends1, dihedrals1 = get_ind_sets(initial)
        bonds2, bends2, dihedrals2 = get_ind_sets(final)
        # Form new superset of coordinate definitions that contain
        # all definitions from geom1 and geom2.
        all_bonds = bonds1 | bonds2
        all_bends = bends1 | bends2
        all_dihedrals = dihedrals1 | dihedrals2
        all_prim_indices = (all_bonds, all_bends, all_dihedrals)
        # Check if internal coordinates that are only present in one
        # of the two geometries are valid in the other. If not we omit
        # this coordinate definition in the end.
        redundant = RedundantCoords(geom1.atoms, geom1.cart_coords,
                                    prim_indices=all_prim_indices)
        bonds, bends, dihedrals = redundant.prim_indices
        return to_set(bonds), to_set(bends), to_set(dihedrals)

    bonds1, bends1, dihedrals1 = merge_coordinate_definitions(initial, final)
    bonds2, bends2, dihedrals2 = merge_coordinate_definitions(final, initial)
    # Only use primitive coordinate definitions that are valid for both
    valid_bonds = bonds1 & bonds2
    valid_bends = bends1 & bends2
    valid_dihedrals = dihedrals1 & dihedrals2
    prim_indices = (valid_bonds, valid_bends, valid_dihedrals)
    print("union of primitives", len(valid_bonds) + len(valid_bends) + len(valid_dihedrals))

    geom1 = Geometry(initial.atoms, initial.cart_coords,
                     coord_type="redund", coord_kwargs={"prim_indices": prim_indices,}
    )
    geom2 = Geometry(final.atoms, final.cart_coords,
                     coord_type="redund", coord_kwargs={"prim_indices": prim_indices,}
    )

    def update_internals(prev_internals, new_internals, bonds_bends, d):
        internal_diffs = np.array(new_internals - prev_internals)
        dihedral_diffs = internal_diffs[bonds_bends:]
        # Find differences that are shifted by 2*pi
        shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2*np.pi) < np.pi/2
        new_dihedrals = new_internals[bonds_bends:]
        new_dihedrals[shifted_by_2pi] -= 2*np.pi * np.sign(dihedral_diffs[shifted_by_2pi])

        new_internals[bonds_bends:] = new_dihedrals
        return new_internals

    def get_tangent(prims1, prims2, bonds_bends):
        diff = prims2 - prims1
        diheds = diff[bonds_bends:].copy()
        diheds_plus = diheds.copy() + 2*np.pi
        diheds_minus = diheds.copy() - 2*np.pi
        bigger = np.abs(diheds) > np.abs(diheds_plus)
        diheds[bigger] = diheds_plus[bigger]
        bigger = np.abs(diheds) > np.abs(diheds_minus)
        diheds[bigger] = diheds_minus[bigger]
        diff[bonds_bends:] = diheds
        return diff

    bonds_bends = len(valid_bonds) + len(valid_bends)
    initial_tangent = get_tangent(geom1.coords, geom2.coords, bonds_bends)
    initial_diff = np.linalg.norm(initial_tangent)
    approx_stepsize = initial_diff / (between+1)
    final_prims = geom2.internal.prim_coords

    geoms = [geom1, ]
    for i in range(between+1):
        new_geom = geoms[-1].copy()
        prim_tangent = get_tangent(new_geom.coords, final_prims, bonds_bends)
        # Form active set
        B = new_geom.internal.B_prim
        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)
        U = eigvectors[:, np.abs(eigvals) > 1e-6]
        reduced_tangent = (np.einsum("i,ij->j", prim_tangent, U) * U).sum(axis=1)
        reduced_tangent /= np.linalg.norm(reduced_tangent)
        step = approx_stepsize * reduced_tangent
        new_coords = new_geom.coords + step
        new_geom.coords = new_coords
        geoms.append(new_geom)
    return geoms


def test_redund():
    initial = geom_from_xyz_file("bare_split.image_000.xyz", coord_type="redund")
    final = geom_from_xyz_file("bare_split.image_056.xyz", coord_type="redund")

    from pysisyphus.interpolate.helpers import interpolate_all

    geoms = interpolate_all((initial, final), 18, kind="redund", align=True)
    out_fn = "dlc_interpolate.trj"
    write_geoms_to_trj(geoms, out_fn)
    print("Wrote ", out_fn)

    # geoms_ref = dlc_interpolate(initial, final)
    # out_fn = "dlc_interpolate_ref.trj"
    # write_geoms_to_trj(geoms_ref, out_fn)
    # print("Wrote ", out_fn)


if __name__ == "__main__":
    # test_lst()
    # test_idpp()
    test_redund()

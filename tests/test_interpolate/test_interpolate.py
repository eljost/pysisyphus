#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.interpolate.LST import LST
from pysisyphus.interpolate.IDPP import IDPP
from pysisyphus.interpolate import interpolate_all


def test_lst():
    initial = geom_from_library("dipeptide_init.xyz")
    final = geom_from_library("dipeptide_fin.xyz")

    geoms = (initial, final)
    lst = LST(geoms, 50, align=True)
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


def test_redund():
    import numpy as np

    from pysisyphus.InternalCoordinates import RedundantCoords
    from pysisyphus.Geometry import Geometry
    from pysisyphus.xyzloader import write_geoms_to_trj

    np.set_printoptions(suppress=True, precision=4)

    # initial = geom_from_xyz_file("bare_split.image_000.xyz", coord_type="redund")
    # final = geom_from_xyz_file("bare_split.image_056.xyz", coord_type="redund")

    # initial = geom_from_xyz_file("01_ed.xyz", coord_type="redund")
    # final = geom_from_xyz_file("01_prod.xyz", coord_type="redund")
    # initial = geom_from_xyz_file("h2o2_hf_321g_opt.xyz", coord_type="redund")
    # final = geom_from_xyz_file("h2o2_rot.xyz", coord_type="redund")
    # initial = geom_from_xyz_file("min115.xyz", coord_type="redund")
    # final = geom_from_xyz_file("plu115.xyz", coord_type="redund")

    # initial = geom_from_xyz_file("hcn.xyz", coord_type="redund")
    # final = geom_from_xyz_file("nhc.xyz", coord_type="redund")

    # initial = geom_from_xyz_file("h2_sih2_start.xyz", coord_type="redund")
    # final = geom_from_xyz_file("h2_sih2_end.xyz", coord_type="redund")

    # initial = geom_from_library("dipeptide_init.xyz", coord_type="redund")
    # final = geom_from_library("dipeptide_fin.xyz", coord_type="redund")

    initial = geom_from_xyz_file("butadiene_ethene.xyz", coord_type="redund")
    final = geom_from_xyz_file("cyclohexene.xyz", coord_type="redund")
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
                     coord_type="redund", prim_indices=prim_indices
    )
    geom2 = Geometry(final.atoms, final.cart_coords,
                     coord_type="redund", prim_indices=prim_indices
    )
    d = geom2.coords - geom1.coords
    c1 = geom1.coords.copy()
    c2 = geom2.coords.copy()

    def update_internals(prev_internals, new_internals, bonds_bends, d):
        internal_diffs = np.array(new_internals - prev_internals)
        dihedral_diffs = internal_diffs[bonds_bends:]
        # Find differences that are shifted by 2*pi
        shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2*np.pi) < np.pi/2
        new_dihedrals = new_internals[bonds_bends:]
        new_dihedrals[shifted_by_2pi] -= 2*np.pi * np.sign(dihedral_diffs[shifted_by_2pi])

        new_internals[bonds_bends:] = new_dihedrals
        return new_internals

    bonds_bends = len(valid_bonds) + len(valid_bends)
    d_diheds = d[bonds_bends:]
    coords2_ = update_internals(geom1.coords, geom2.coords, bonds_bends, d)
    d_ = coords2_ - geom1.coords
    num = 20
    bl = len(valid_bonds)
    ba = len(valid_bends)
    bd = len(valid_dihedrals)
    bo1 = geom1.coords[:bl]
    bo2 = coords2_[:bl]
    be1 = geom1.coords[bl:bl+ba]
    be2 = coords2_[bl:bl+ba]
    d1 = geom1.coords[bl+ba:]
    d2 = coords2_[bl+ba:]
    import pdb; pdb.set_trace()
    base_step = d_ / (num-1)
    # import pdb; pdb.set_trace()
    geoms = [geom1, ]
    # import pdb; pdb.set_trace()
    print("base_step", base_step)
    for i in range(num):
        print(i)
        # step = base_step
        # if i == 11:
            # import pdb; pdb.set_trace()
        new_geom = geoms[-1].copy()
        try:
            new_coords = new_geom.coords + base_step
        except ValueError:
            import pdb; pdb.set_trace()
        new_geom.coords = new_coords
        # print(i, new_coords)
        geoms.append(new_geom)
        write_geoms_to_trj(geoms, f"redund_{i:02d}.trj")
    write_geoms_to_trj(geoms, f"liiic_{i:02d}.trj")


if __name__ == "__main__":
    # test_lst()
    # test_idpp()
    test_redund()

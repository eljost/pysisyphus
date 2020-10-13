import argparse
import itertools as it
from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
from natsort import natsorted
import numpy as np
import rmsd as rmsd
import yaml

from pysisyphus.constants import BOHR2ANG, AU2KJPERMOL
from pysisyphus.cos import *
from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords.setup import get_fragments
from pysisyphus.helpers import geom_loader, procrustes, get_coords_diffs, shake_coords
from pysisyphus.interpolate import *
from pysisyphus.intcoords.helpers import form_coordinate_union
from pysisyphus.io.pdb import geom_to_pdb_str
from pysisyphus.stocastic.align import match_geom_atoms
from pysisyphus.xyzloader import split_xyz_str


INTERPOLATE = {
    "idpp": IDPP.IDPP,
    "lst": LST.LST,
    "linear": Interpolator.Interpolator,
    "redund": Redund.Redund,
}


def parse_args(args):
    parser = argparse.ArgumentParser("Utility to transform .xyz and .trj files.")

    parser.add_argument("fns", nargs="+",
        help="Filenames of .xyz and/or .trj files (xyz and trj can be mixed)."
    )

    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--between", type=int, default=0,
        help="Interpolate additional images."
    )
    action_group.add_argument("--align", action="store_true",
        help="Align geometries onto the first geometry."
    )
    action_group.add_argument("--split", action="store_true",
        help="Split a supplied geometries in multiple .xyz files."
    )
    action_group.add_argument("--reverse", action="store_true",
        help="Reverse a .trj file."
    )
    action_group.add_argument("--cleantrj", action="store_true",
        help="Keep only the first four columns of xyz/trj files."
    )
    action_group.add_argument("--spline", action="store_true",
        help="Evenly redistribute geometries along a splined path."
    )
    action_group.add_argument("--first", type=int,
        help="Copy the first N geometries to a new .trj file."
    )
    action_group.add_argument("--every", type=int,
        help="Create new .trj with every N-th geometry. "
             "Always includes the first and last point."
    )
    action_group.add_argument("--center", action="store_true",
        help="Move the molecules centroid into the origin."
    )
    action_group.add_argument("--centerm", action="store_true",
        help="Move the molecules center of mass into the origin."
    )
    action_group.add_argument("--translate", nargs=3, type=float,
        help="Translate the molecule by the given vector given " \
             "in Ångström."
    )
    action_group.add_argument("--append", action="store_true",
        help="Combine the given .xyz files into one .xyz file."
    )
    action_group.add_argument("--join", action="store_true",
        help="Combine the given .xyz/.trj files into one .trj file."
    )
    action_group.add_argument("--match", action="store_true",
        help="Resort the second .xyz file so the atom order matches the "
             "first .xyz file. Uses the hungarian method."
    )
    action_group.add_argument("--std", action="store_true",
        help="Move supplied geometry to its standard orientation."
    )
    action_group.add_argument("--shake", action="store_true",
        help="Shake (randomly displace) coordiantes."
    )
    action_group.add_argument("--internals", action="store_true",
        help="Print automatically generated internals."
    )
    action_group.add_argument("--get", type=int,
        help="Get n-th geometry. Expects 0-based index input."
    )
    action_group.add_argument("--geti", action="store_true",
        help="Decide on geometry interactively."
    )
    action_group.add_argument("--origin", action="store_true",
        help="Translate geometry, so that min(X/Y/Z) == 0."
    )
    action_group.add_argument("--fragsort", action="store_true",
        help="Resort atoms by fragments."
    )
    action_group.add_argument("--topdb", action="store_true",
        help="Convert given geometry to PDB with automatic fragment detection."
    )

    shake_group = parser.add_argument_group()
    shake_group.add_argument("--scale", type=float, default=0.1,
        help="Scales the displacement in --shake."
    )
    shake_group.add_argument("--seed", type=int, default=None,
        help="Initialize the RNG for reproducible results."
    )

    interpolate_group = parser.add_mutually_exclusive_group()
    interpolate_group.add_argument("--idpp", action="store_true",
        help="Interpolate using Image Dependent Pair Potential."
    )
    interpolate_group.add_argument("--lst", action="store_true",
        help="Interpolate by linear synchronous transit."
    )
    interpolate_group.add_argument("--redund", action="store_true",
        help="Interpolate in internal coordinates."
    )

    parser.add_argument("--noipalign", action="store_false",
        help="Don't align geometries when interpolating."
    )
    parser.add_argument("--bohr", action="store_true",
        help="Input geometries are in Bohr instead of Angstrom."
    )
    parser.add_argument("--noxyz", action="store_false",
        help="Disable dumping of single .xyz files."
    )
    parser.add_argument("--atoms", nargs="+", type=int, default=list(),
        help="Used with --internals. Only print primitives including the given atoms."
    )
    parser.add_argument("--add_prims", type=str, default="",
        help="Used with --internals. Define additional primitives. Expects a "
             "string representation of a nested list that can be parsed as YAML "
             "e.g. [[10,30],[1,2,3],[4,5,6,7]]."
    )

    return parser.parse_args()


def read_geoms(xyz_fns, in_bohr=False, coord_type="cart",
               define_prims=None, typed_prims=None):
    if isinstance(xyz_fns, str):
        xyz_fns = [xyz_fns, ]

    geoms = list()
    geom_kwargs = {
        "coord_type": coord_type,
    }
    # Dont supply coord_kwargs with coord_type == "cart"
    if coord_type != "cart":
        geom_kwargs["coord_kwargs"] = {
            "define_prims": define_prims,
            "typed_prims": typed_prims,
        }

    for fn in xyz_fns:
        if "*" in fn:
            cwd = Path(".")
            geom = [geom_loader(xyz_fn, **geom_kwargs)
                    for xyz_fn in natsorted(cwd.glob(fn))]
        # Simplify this to use geom_loader...
        elif fn.endswith(".xyz"):
            geom = [geom_loader(fn, **geom_kwargs), ]
        elif fn.endswith(".trj"):
            geom = geom_loader(fn, **geom_kwargs)
        elif fn.endswith(".pdb"):
            geom = [geom_loader(fn, **geom_kwargs), ]
        elif fn.endswith(".cjson"):
            geom = [geom_loader(fn, **geom_kwargs), ]
        else:
            continue
        geoms.extend(geom)

    # Try to parse as inline xyz formatted string
    if len(geoms) == 0:
        try:
            atoms_coords = split_xyz_str(fn)
            # We excpect the coordinates to be given in Angstrom
            geoms = [Geometry(atoms, coords/BOHR2ANG, **geom_kwargs)
                     for atoms, coords in atoms_coords]
        except AssertionError:
            raise Exception("Could not parse supplied 'xyz' values as either "
                            ".xyz, .trj or xyz-formatted string.!"
            )

    # Original coordinates are in bohr, but pysisyphus expects them
    # to be in Angstrom, so right now they are already multiplied
    # by ANG2BOHR. We revert this.
    if in_bohr:
        for geom in geoms:
            geom.coords *= BOHR2ANG
    return geoms


def get_geoms(xyz_fns, interpolate=None, between=0, coord_type="cart",
              comments=False, in_bohr=False, define_prims=None, union=None,
              interpolate_align=True, same_prims=True, quiet=False):
    """Returns a list of Geometry objects in the given coordinate system
    and interpolates if necessary."""

    if union is not None:
        assert coord_type != "cart", \
            "union must not be used with coord_type=cart!"
        union_geoms = read_geoms(union, coord_type=coord_type)
        assert len(union_geoms) == 2, \
            f"Got {len(union_geoms)} geometries for 'union'! Please give only two!"
        typed_prims_union = form_coordinate_union(*union_geoms)
    else:
        typed_prims_union = None

    assert interpolate in list(INTERPOLATE.keys()) + [None], \
        f"Unsupported type: '{interpolate}' given. Valid arguments are " \
        f"{list(INTERPOLATE.keys())}'"

    geoms = read_geoms(xyz_fns, in_bohr, coord_type=coord_type,
                       define_prims=define_prims, typed_prims=typed_prims_union)
    if not quiet:
        print(f"Read {len(geoms)} geometries.")

    atoms_0 = geoms[0].atoms
    atoms_strs = [" ".join(geom.atoms).lower() for geom in geoms]
    atoms_0_str = atoms_strs[0]
    assert all([atoms_str == atoms_0_str for atoms_str in atoms_strs]), \
        "Atom ordering/numbering in the geometries is inconsistent!"

    # Dont try to align 1-atom species
    interpolate_align = interpolate_align and len(geoms[0].atoms) > 1

    # TODO:
    # Multistep interpolation (when more than two geometries are specified)
    # in internal coordinates may lead to a different number of defined coordinates.
    # Maybe the union between geom0 and geom1 contains 6 internals and the union
    # betweeen geom1 and geom2 contains 8 primtives. Then the number of coordinates
    # at all images in the final list may be non-constant.
    if interpolate:
        interpolate_class = INTERPOLATE[interpolate]
        interpolator = interpolate_class(geoms, between, align=interpolate_align)
        geoms = interpolator.interpolate_all()

    # Recreate Geometries so they have the correct coord_type. There may
    # be a difference between the coord_type used for interpolation and
    # the desired coord_type as specified in the function arguments.
    # prim_indices = [None for geom in geoms]
    if coord_type != geoms[0].coord_type:
        recreated_geoms = list()
        for geom in geoms:
            coord_kwargs = None
            if coord_type != "cart":
                coord_kwargs = {
                    "typed_prims": geom.internal.typed_prims,
                }
            geom = Geometry(geom.atoms, geom.cart_coords, coord_type=coord_type,
                            comment=geom.comment, coord_kwargs=coord_kwargs
            )
            recreated_geoms.append(geom)
        geoms = recreated_geoms

    if not same_prims:
        return geoms

    same_prim_inds = True
    if coord_type != "cart":
        geom_prim_inds = [geom.internal.typed_prims for geom in geoms]
        first_set = geom_prim_inds[0]
        same_prim_inds = all([ith_set == first_set for ith_set in geom_prim_inds[1:]])
    # Recreate geometries with the same primitive internal coordinates
    if not same_prim_inds:
        typed_prims = form_coordinate_union(geoms[0], geoms[-1])
        geoms = [Geometry(atoms_0, geom.cart_coords,
                         coord_type=coord_type,
                         coord_kwargs={"typed_prims": typed_prims},
                 )
                 for geom in geoms]
        # assert all([geom.coords.size == geoms[0].coords.size for geom in geoms])
        assert all([len(geom.internal.typed_prims) == len(typed_prims) for geom in geoms])
    return geoms


def dump_geoms(geoms, fn_base, trj_infix="", dump_trj=True, dump_xyz=True, dump_pdb=False,
               ang=False):
    xyz_per_geom = [geom.as_xyz() for geom in geoms]
    if dump_trj:
        trj_str = "\n".join(xyz_per_geom)
        trj_fn = f"{fn_base}{trj_infix}.trj"
        with open(trj_fn, "w") as handle:
            handle.write(trj_str)
        print(f"Wrote all geometries to {trj_fn}.")
    if dump_xyz:
        for i, xyz in enumerate(xyz_per_geom):
            geom_fn = f"{fn_base}.geom_{i:03d}.xyz"
            with open(geom_fn, "w") as handle:
                handle.write(xyz)
            print(f"Wrote geom {i:03d} to {geom_fn}.")
    elif dump_pdb:
        for i, geom in enumerate(geoms):
            geom_fn = f"{fn_base}.geom_{i:03d}.pdb"
            pdb_str = geom_to_pdb_str(geom, detect_fragments=True)
            with open(geom_fn, "w") as handle:
                handle.write(pdb_str)
            print(f"Wrote geom {i:03d} to {geom_fn}.")
    print()


def align(geoms):
    """Align all geometries onto the first using partical procrustes."""
    cos = ChainOfStates.ChainOfStates(geoms)
    procrustes(cos)
    return [geom for geom in cos.images]


def spline_redistribute(geoms):
    szts = SimpleZTS.SimpleZTS(geoms)
    pre_diffs = get_coords_diffs([image.coords for image in szts.images])
    szts.reparametrize()
    post_diffs = get_coords_diffs([image.coords for image in szts.images])
    cds_str = lambda cds: " ".join([f"{cd:.2f}" for cd in cds])
    print("Normalized path segments before splining:")
    print(cds_str(pre_diffs))
    print("Normalized path segments after redistribution along spline:")
    print(cds_str(post_diffs))
    return szts.images


def every(geoms, every_nth):
    # every_nth_geom = geoms[::every_nth]
    # The first geometry is always present, but the last geometry
    # may be missing.
    sampled_indices = list(range(0, len(geoms), every_nth))
    if sampled_indices[-1] != len(geoms)-1:
        sampled_indices.append(len(geoms)-1)
    sampled_inds_str = ", ".join([str(i) for i in sampled_indices])
    print(f"Sampled indices {sampled_inds_str}")
    # if every_nth_geom[-1] != geoms[-1]:
        # every_nth_geom.append(geoms[-1])
    every_nth_geom = [geoms[i] for i in sampled_indices]
    return every_nth_geom


def center(geoms):
    for geom in geoms:
        geom.coords3d = geom.coords3d - geom.centroid
    return geoms


def centerm(geoms):
    for geom in geoms:
        geom.coords3d = geom.coords3d - geom.center_of_mass
    return geoms


def translate(geoms, trans):
    for geom in geoms:
        geom.coords3d += trans
    return geoms


def append(geoms):
    atoms = geoms[0].atoms * len(geoms)
    coords = list(it.chain([geom.coords for geom in geoms]))
    return [Geometry(atoms, coords), ]


def match(ref_geom, geom_to_match):
    rmsd_before = rmsd.kabsch_rmsd(ref_geom.coords3d, geom_to_match.coords3d)
    print(f"Kabsch RMSD before: {rmsd_before:.4f}")
    matched_geom = match_geom_atoms(ref_geom, geom_to_match, hydrogen=True)

    # Right now the atoms are not in the right order as we only sorted the
    # individual coord blocks by atom.
    # This dictionary will hold the counter indices for the individual atom
    atom_type_inds = {atom: 0 for atom in ref_geom.atom_types}
    matched_coord_blocks, _ = matched_geom.coords_by_type
    new_coords = list()
    for atom in ref_geom.atoms:
        # Get the current counter/index from the dicitonary for the given atom
        cur_atom_ind = atom_type_inds[atom]
        # Select the appropriate atom from the coords block
        atom_coords = matched_coord_blocks[atom][cur_atom_ind]
        new_coords.append(atom_coords)
        # Increment the counter so the next time the same atom type comes up
        # we fetch the next entry of the coord block.
        atom_type_inds[atom] += 1
    # Assign the updated atom order and corresponding coordinates
    matched_geom.atoms = ref_geom.atoms
    matched_geom.coords = np.array(new_coords).flatten()
    rmsd_after = rmsd.kabsch_rmsd(ref_geom.coords3d, matched_geom.coords3d)
    print(f"Kabsch RMSD after: {rmsd_after:.4f}")
    return [matched_geom, ]


def standard_orientation(geoms):
    [geom.standard_orientation() for geom in geoms]
    return geoms


def shake(geoms, scale=0.1, seed=None):
    for geom in geoms:
        geom.coords = shake_coords(geom.coords, scale=scale, seed=seed)
    return geoms


def print_internals(geoms, filter_atoms=None, add_prims=""):
    if filter_atoms is None:
        filter_atoms = set()

    filter_set = set(filter_atoms)
    pi_types = {2: "B", 3: "A", 4: "D"}
    add_prims = yaml.safe_load(add_prims)

    for i, geom in enumerate(geoms):
        print(geom)
        atom_num = len(geom.atoms)
        atom_inds = set(range(atom_num))
        # Atom indices must superset of filter_atoms
        invalid = filter_set - atom_inds
        assert not invalid, \
            f"Filter indices {invalid} are outside of the " \
            f"valid range for the {i}-th geometry '{geom}' with {atom_num} " \
            f"atoms (valid indices: range(0,{atom_num}))."

        int_geom = Geometry(geom.atoms, geom.coords, coord_type="redund",
                            coord_kwargs={"define_prims": add_prims,},
        )

        prim_counter = 0
        prev_len = 2
        for j, pi in enumerate(int_geom.internal._prim_internals):
            pi_type = pi_types[len(pi.inds)]
            val = pi.val
            len_ = len(pi.inds)
            if len_ > 2:
                unit = "°"
                val = np.rad2deg(val)
            else:
                val *= BOHR2ANG
                unit = " Å"

            if len_ > prev_len:
                prim_counter = 0
                print()

            # Continue if we want to filter and there is no intersection
            # between the atoms of the current primitive and the atoms we
            # want to filter for.
            if filter_set and not (set(pi.inds) & filter_set):
                continue
            print(f"{j:04d}: {pi_type}{prim_counter:03d}{str(pi.inds): >20} {val: >10.4f}"
                  f"{unit}")
            prim_counter += 1
            prev_len = len_

        print(f"Printed {j+1} primitive internals.")
        print()


def get(geoms, index):
    # Convert to positive index. Right now this doesn't do anything useful.
    # Could be used to generated more meaningful filenames if we could somehow
    # also return a index number to generate the filename.
    if index < 0:
        index += len(geoms)
    return [geoms[index], ]


class GotNoGeometryException(Exception):
    pass


def get_interactively(geoms):
    # Try to parse energies from geoms
    energy_re = re.compile("[-\.\d]+")
    energies = list()
    for geom in geoms:
        mobj = energy_re.search(geom.comment)
        if mobj:
            energy = float(mobj[0])
        else:
            energy = np.nan
        energies.append(energy)
    energies = np.array(energies)
    energies -= np.nanmin(energies)
    energies *= AU2KJPERMOL
    min_ind = np.nanargmin(energies)
    print(f"Minimum energy at index {min_ind}")
    print("(q) to quit\n(p) to plot energies")
    msg = f"Input index (0-{len(geoms)-1})/p/q: "
    while True:
        try:
            selection = input(msg)
            if selection == "q":
                raise GotNoGeometryException()
            elif selection == "p":
                fig, ax = plt.subplots()
                ax.plot(energies, "o-")
                ax.set_xlabel("Index")
                ax.set_ylabel("ΔE / kJ mol⁻¹")
                plt.show()
                continue
            selection = int(selection)
        except ValueError:
            print("Invalid input!")
            continue
        en = energies[selection]
        print(f"ΔE at geometry {selection} is {en:+.2f} kJ mol⁻¹.")
        yn = input("Get this geometry (y/n)? ").lower()
        if yn == "y":
            return get(geoms, selection)
        else:
            continue


def origin(geoms):
    for i, geom in enumerate(geoms):
        print(f"{i:02d}: {geom}")
        geom.coords3d -= geom.coords3d.min(axis=0)
        print(f"\tmax(coords3d): {geom.coords3d.max(axis=0)}")
    return geoms


def frag_sort(geoms):
    sorted_geoms = list()
    for i, geom in enumerate(geoms):
        print(f"{i:02d}: {geom}")
        frags = get_fragments(geom.atoms, geom.coords)
        print(f"\tFound {len(frags)} fragments.\n"
               "\tResorting atoms and coordinates.")
        new_indices = list(it.chain(*frags))
        new_atoms = [geom.atoms[i] for i in new_indices]
        new_coords = geom.coords3d[new_indices]
        sorted_geoms.append(
                Geometry(new_atoms, new_coords)
        )
    return sorted_geoms


def run():
    args = parse_args(sys.argv[1:])

    if args.idpp:
        interpolate = "idpp"
    elif args.lst:
        interpolate = "lst"
    elif args.redund:
        interpolate = "redund"
    elif args.between:
        interpolate = "linear"
    else:
        interpolate = None

    # Read supplied files and create Geometry objects
    geoms = get_geoms(args.fns, interpolate, args.between, in_bohr=args.bohr,
                      interpolate_align=args.noipalign)

    to_dump = geoms
    dump_trj = True
    dump_xyz = args.noxyz
    dump_pdb = False
    trj_infix = ""
    if args.between:
        fn_base = "interpolated"
    elif args.align:
        to_dump = align(geoms)
        fn_base = "aligned"
    elif args.split:
        fn_base = "split"
        dump_trj = False
    elif args.reverse:
        to_dump = geoms[::-1]
        fn_base = "reversed"
    elif args.cleantrj:
        fn_base = "cleaned"
    elif args.first:
        to_dump = geoms[:args.first]
        fn_base = "first"
        trj_infix = f"_{args.first}"
    elif args.spline:
        to_dump = spline_redistribute(geoms)
        fn_base = "splined"
    elif args.every:
        to_dump = every(geoms, args.every)
        fn_base = "every"
        trj_infix = f"_{args.every}th"
    elif args.center:
        to_dump = center(geoms)
        fn_base = "centered"
    elif args.centerm:
        to_dump = centerm(geoms)
        fn_base = "centeredm"
    elif args.translate:
        trans = np.array(args.translate) / BOHR2ANG
        to_dump = translate(geoms, trans)
        fn_base = "translated"
    elif args.append:
        to_dump = append(geoms)
        fn_base = "appended"
    elif args.join:
        to_dump = geoms
        fn_base = "joined"
    elif args.match:
        to_dump = match(*geoms)
        fn_base = "matched"
    elif args.std:
        to_dump = standard_orientation(geoms)
        fn_base = "standard"
    elif args.shake:
        to_dump = shake(geoms, args.scale, args.seed)
        fn_base = "shaked"
    elif args.get or (args.get == 0):
        to_dump = get(geoms, args.get)
        fn_base = "got"
    elif args.geti:
        try:
            to_dump = get_interactively(geoms)
            fn_base = "got"
        except GotNoGeometryException:
            return
    elif args.internals:
        print_internals(geoms, args.atoms, args.add_prims)
        return
    elif args.fragsort:
        to_dump = frag_sort(geoms)
        fn_base = "frag_sorted"
    elif args.origin:
        origin(geoms)
        fn_base = "origin"
    elif args.topdb:
        fn_base = "as_pdb"
        dump_pdb = True
        dump_xyz = False

    # Write transformed geometries
    dump_trj = dump_trj and (len(to_dump) > 1)

    dump_geoms(to_dump, fn_base, trj_infix=trj_infix, dump_trj=dump_trj,
               dump_xyz=dump_xyz, dump_pdb=dump_pdb)

if __name__ == "__main__":
    run()

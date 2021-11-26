import argparse
import sys

import numpy as np

from pysisyphus.calculators import Composite, HardSphere, TransTorque
from pysisyphus.calculators.OBabel import OBabel
from pysisyphus.drivers.precon_pos_rot import (
    center_fragments,
    get_steps_to_active_atom_mean,
    form_A,
    get_which_frag,
    SteepestDescent,
)
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import align_coords, geom_loader
from pysisyphus.io import geom_to_crd_str
from pysisyphus.linalg import get_rot_mat_for_coords
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.xyzloader import coords_to_trj


def merge_geoms(geom1, geom2, geom1_del=None, geom2_del=None, make_bonds=None):
    """Merge geom1 and geom2 while keeping the original coordinates.

    Supports deleting certain atoms.
    """

    if geom1_del is None:
        geom1_del = list()
    if geom2_del is None:
        geom2_del = list()

    atom_num1 = len(geom1.atoms)
    geom2_del = np.array(geom2_del)
    # Offset by number of atoms in geom1
    geom2_del += atom_num1

    if make_bonds is not None:
        make_bonds = np.array(make_bonds, dtype=int)
        geom1_bond_inds, geom2_bond_inds = make_bonds.T
        geom2_bond_inds += atom_num1

        # Correct original bond indices given in make_bonds
        to_del = np.concatenate((geom1_del, geom2_del))
        # If geom1_del/geom2_del are empty then there is nothing to do
        for to_del in geom1_del:
            mask = geom1_bond_inds > to_del
            geom1_bond_inds[mask] -= 1
            geom2_bond_inds -= 1

        for to_del in geom2_del:
            mask = geom2_bond_inds > to_del
            geom2_bond_inds[mask] -= 1
        make_bonds_cor = np.stack((geom1_bond_inds, geom2_bond_inds), axis=1)
    else:
        make_bonds_cor = None

    union = geom1 + geom2
    union = union.del_atoms(list(geom1_del) + list(geom2_del))

    # Set appropriate fragments
    atom_num1 -= len(geom1_del)
    frag1 = np.arange(atom_num1)
    lig_atoms = geom2.atoms
    frag2 = np.arange(atom_num1, atom_num1 + len(lig_atoms) - len(geom2_del))
    union.fragments = {"geom1": frag1.tolist(), "geom2": frag2.tolist()}

    return union, make_bonds_cor


def hardsphere_merge(geom1, geom2):
    union = geom1 + geom2

    atom_num = len(geom1.atoms)
    # Set up lists containing the atom indices for the two fragments
    frag_lists = [
        [i for i, _ in enumerate(geom1.atoms)],
        [atom_num + i for i, _ in enumerate(geom2.atoms)],
    ]

    def get_hs(kappa=1.0):
        return HardSphere(
            union,
            frag_lists,
            permutations=True,
            kappa=kappa,
        )

    union.set_calculator(get_hs(1.0))
    opt_kwargs = {
        "max_cycles": 1000,
        "max_step": 0.5,
        "rms_force": 0.0005,
    }
    opt = SteepestDescent(union, **opt_kwargs)
    opt.run()
    return union


def prepare_merge(geom1, bond_diff, geom2=None, del1=None, del2=None, dump=False):
    if del1:
        geom1 = geom1.del_atoms(del1)
    if del2:
        geom2 = geom2.del_atoms(del2)

    if geom2:
        union = geom1 + geom2
        atom_num = len(geom1.atoms)
        # Set up lists containing the atom indices for the two fragments
        frag_lists = [
            [i for i, _ in enumerate(geom1.atoms)],
            [atom_num + i for i, _ in enumerate(geom2.atoms)],
        ]
        fragments = {"geom1": frag_lists[0], "geom2": frag_lists[1]}
        union.fragments = fragments
    else:
        union = geom1
    frag_lists = [union.fragments["geom1"], union.fragments["geom2"]]

    backup = list()

    def keep(comment):
        backup.append((union.cart_coords.copy(), comment))

    keep("Initial union")

    which_frag = get_which_frag(frag_lists)
    AR = form_A(frag_lists, which_frag, bond_diff)

    # Center fragments at their geometric average
    center_fragments(frag_lists, union)
    keep("Centered fragments")

    # Translate centroid of active atoms to origin
    alphas = get_steps_to_active_atom_mean(frag_lists, frag_lists, AR, union.coords3d)
    for frag, alpha in zip(frag_lists, alphas):
        union.coords3d[frag] += alpha
    keep("Shifted centroids to active atoms")

    def get_hs(kappa=1.0):
        return HardSphere(
            union,
            frag_lists,
            permutations=True,
            kappa=kappa,
        )

    union.set_calculator(get_hs(1.0))
    opt_kwargs = {
        "max_cycles": 1000,
        "max_step": 0.5,
        "rms_force": 0.0005,
    }
    opt = SteepestDescent(union, **opt_kwargs)
    opt.run()
    keep("Initial hardsphere optimization")

    keys_calcs = {
        "hs": get_hs(10.0),
        "tt": TransTorque(frag_lists, frag_lists, AR, AR, kappa=2, do_trans=True),
    }
    comp = Composite("hs + tt", keys_calcs, remove_translation=True)

    union.set_calculator(comp)
    max_cycles = 15_000
    max_dist = 50
    max_step = max_dist / max_cycles
    opt_kwargs = {
        "max_cycles": max_cycles,
        "max_step": max_step,
        "rms_force": 0.1,
        # "dump": True,
    }
    opt = SteepestDescent(union, **opt_kwargs)
    opt.run()
    keep("Hardsphere + TransTorque optimization")

    if dump:
        kept_coords, kept_comments = zip(*backup)
        align_coords(kept_coords)
        fn = "merge_dump.trj"
        coords_to_trj(fn, union.atoms, kept_coords, comments=kept_comments)

    return union


def merge_opt(union, bond_diff, ff="uff"):
    """Fragment merging along given bond by forcefield optimization."""
    from openbabel import pybel

    geom1 = union.get_fragments("geom1")
    freeze = list(range(len(geom1.atoms)))

    # Create pybel.Molecule/OBMol to set the missing bonds
    mol = pybel.readstring("xyz", union.as_xyz())
    obmol = mol.OBMol
    for from_, to_ in bond_diff:
        obmol.AddBond(int(from_ + 1), int(to_ + 1), 1)

    # Use modified pybel.Molecule
    calc = OBabel(mol=mol, ff=ff)
    funion = union.copy()
    # Only releax second fragment
    funion.freeze_atoms = freeze
    funion.set_calculator(calc)

    opt = LBFGS(funion, max_cycles=1000, max_step=0.5, dump=False)
    opt.run()

    return funion


def align_on_subset(geom1, union, del1=None):
    """Align 'union' onto subset of 'geom1'"""

    # Delete some coordinates (atoms) in geom1, that are not present in union
    coords3d_1 = geom1.coords3d.copy()
    atoms1 = geom1.atoms
    if del1 is not None:
        atoms1 = [atom for i, atom in enumerate(atoms1) if i not in del1]
        coords3d_1 = np.delete(coords3d_1, del1, axis=0)
    atoms1 = tuple(atoms1)
    num1 = coords3d_1.shape[0]
    # Restrict length of union to length of coords3d_1
    coords3d_2 = union.coords3d.copy()
    coords3d_2_subset = coords3d_2[:num1]
    atoms2 = union.atoms
    assert atoms1 == tuple(atoms2[:num1])

    *_, rot_mat = get_rot_mat_for_coords(coords3d_1, coords3d_2_subset)

    # Align merged system
    coords3d_2_aligned = (coords3d_2 - coords3d_2.mean(axis=0)[None, :]).dot(rot_mat)

    # Translate aligned system so that centroids of subsets match
    centroid_1 = coords3d_1.mean(axis=0)
    centroid_2 = coords3d_2_aligned[:num1].mean(axis=0)
    coords3d_2_aligned += centroid_1 - centroid_2

    aligned = Geometry(atoms2, coords3d_2_aligned)
    subset = Geometry(atoms2[:num1], coords3d_2_aligned[:num1])
    rest = Geometry(atoms2[num1:], coords3d_2_aligned[num1:])
    return aligned, subset, rest


def merge_with_frozen_geom(frozen_geom, lig_geom, make_bonds, frozen_del, lig_del):
    union, make_bonds_cor = merge_geoms(
        frozen_geom, lig_geom, frozen_del, lig_del, make_bonds
    )
    atoms = union.atoms
    print("Docking to form bonds:")
    for i, (from_, to_) in enumerate(make_bonds):
        from_atom = atoms[from_]
        to_atom = atoms[to_]
        print(f"\t{i:02d} {from_atom}{from_}-{to_atom}{to_}")

    union = prepare_merge(union, make_bonds_cor, dump=True)
    opt_union = merge_opt(union, make_bonds_cor, ff="uff")
    # aligned: whole system
    # subset: frozen_geom - deleted atoms
    # rest: aligned ligand
    aligned, subset, rest = align_on_subset(frozen_geom, opt_union, frozen_del)
    return aligned, subset, rest


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("frozen_fn", help="Filename of the frozen geometry.")
    parser.add_argument("lig_fn", help="Filename of the ligand to be merged.")
    parser.add_argument(
        "--frozen-del",
        nargs="+",
        type=int,
        help="0-based atom indices to be deleted from frozen_fn.",
    )
    parser.add_argument(
        "--lig-del",
        nargs="+",
        type=int,
        help="0-based atom indices to be deleted from lig_fn.",
    )
    parser.add_argument(
        "--make-bonds",
        nargs="+",
        type=int,
        help="0-based indices of atom pairs (frozen, lig), forming bonds "
        "in the merged geometry. lig indices should start at 0.",
    )
    parser.add_argument("--res", default="LIG")
    parser.add_argument("--resno", type=int)

    return parser.parse_args(args)


def run_merge():
    args = parse_args(sys.argv[1:])

    frozen_fn = args.frozen_fn
    lig_fn = args.lig_fn
    protein_pdb = frozen_fn
    lig_pdb = lig_fn

    frozen_geom = geom_loader(protein_pdb)
    lig_geom = geom_loader(lig_pdb)

    prot_del = args.frozen_del
    lig_del = args.lig_del
    make_bonds = args.make_bonds
    make_bonds = np.array(make_bonds, dtype=int).reshape(-1, 2)
    aligned, subset, rest = merge_with_frozen_geom(
        frozen_geom, lig_geom, make_bonds, prot_del, lig_del
    )
    aligned.jmol()

    trj_fn = "merged.trj"
    trj = "\n".join([geom.as_xyz() for geom in (aligned, rest)])
    with open(trj_fn, "w") as handle:
        handle.write(trj)
    print(f"Dumped geometries to '{trj_fn}'.")

    res = args.res
    resno = args.resno
    if (res is not None) and (resno is not None):
        crd_str = geom_to_crd_str(
            rest, res=res, resno=resno, ref_atoms=lig_geom.atoms, del_atoms=lig_del
        )
        crd_fn = "lig_aligned.crd"
        with open(crd_fn, "w") as handle:
            handle.write(crd_str)
        print(f"Dumped ligand coordinates to  to '{crd_fn}'.")


if __name__ == "__main__":
    run_merge()

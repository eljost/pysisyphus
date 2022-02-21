import argparse
import sys

from pysisyphus import logger
from pysisyphus.calculators import XTB
from pysisyphus.calculators.OBabel import OBabel
from pysisyphus.config import LIB_DIR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import chunks
from pysisyphus.intcoords.setup import get_bond_sets
from pysisyphus.linalg import get_rot_mat_for_coords
from pysisyphus.optimizers.LBFGS import LBFGS


def get_bond_subgeom(geom, ind, invert=False):
    bond_inds = get_bond_sets(geom.atoms, geom.coords3d)
    # Determine bonds that are connected to atom 'ind' in 'geom'
    conn_bond_inds = [bond for bond in bond_inds if ind in bond]
    assert len(conn_bond_inds) == 1
    ind0, ind1 = conn_bond_inds[0]
    # Index of atom to which 'ind' is connected
    anchor_ind = ind0 if (ind == ind1) else ind1
    # conn_bond_ind = (ind0, ind1)
    conn_bond_ind = (ind, anchor_ind) if invert else (anchor_ind, ind)
    return geom.get_subgeom(conn_bond_ind), anchor_ind


def run_opt(geom, freeze_inds, ff="uff", use_xtb=False, charge=0, mult=1):
    def get_xtb_calc():
        return XTB(charge=charge, mult=mult, pal=2)

    if not use_xtb:
        try:
            from openbabel import pybel

            # Create pybel.Molecule/OBMol to set the missing bonds
            mol = pybel.readstring("xyz", geom.as_xyz())
            # Use modified pybel.Molecule
            calc = OBabel(mol=mol, ff=ff)
        except ImportError:
            logger.warning("Could not import openbabel. Trying fallback to GFN2-xTB.")
            calc = get_xtb_calc()
    else:
        calc = get_xtb_calc()

    frozen_geom = geom.copy()
    # Only some atoms
    frozen_geom.freeze_atoms = freeze_inds
    frozen_geom.set_calculator(calc)

    opt = LBFGS(frozen_geom, max_cycles=1000, max_step=0.5, dump=False)
    opt.run()

    return frozen_geom


def replace_atom(
    geom,
    ind,
    repl_geom,
    repl_ind,
    return_full=True,
    opt=False,
    use_xtb=False,
    charge=0,
    mult=1,
):
    """Replace atom with fragment."""

    geom = geom.copy(coord_type="cart")
    if len(repl_geom.atoms) == 1:
        atoms = list(geom.atoms)
        atoms[ind] = repl_geom.atoms[0]
        coords3d = geom.coords3d
        return Geometry(atoms, coords3d)

    repl_geom = repl_geom.copy(coord_type="cart")

    # Determine bond-subgeometries, which will be used to calculate a rotation
    # matrix.
    geom_sub, _ = get_bond_subgeom(geom, ind)
    # Call with invert=True, as we want the bonds to point in opposite directions.
    # Without invert=True, 'get_bond_subgeom' will always return with the atom order
    # (ind, anchor_ind). This would lead to overlapping fragments.
    repl_sub, repl_anchor = get_bond_subgeom(repl_geom, repl_ind, invert=True)

    # Determine step that translates 'repl_geom' so that the atom 'ind' in 'geom'
    # and 'anchor_ind' in 'repl_geom' coincide.
    # step = geom.coords3d[ind] - repl_geom.coords3d[repl_anchor]
    # repl_sub.coords3d += step[None, :]

    # Determine rotation matrix that aligns 'repl_geom' on on 'geom_bond'
    *_, rot_mat = get_rot_mat_for_coords(geom_sub.coords3d, repl_sub.coords3d)
    # Rotate whole 'repl_geom' with 'rot_mat'
    repl_c3d = repl_geom.coords3d
    repl_centroid = repl_c3d.mean(axis=0)
    repl_c3d_rot = (repl_c3d - repl_centroid[None, :]).dot(rot_mat) + repl_centroid[
        None, :
    ]
    repl_geom.coords = repl_c3d_rot
    # Translate repl_geom so that repl_anchor and ind coincide
    step = geom.coords3d[ind] - repl_c3d_rot[repl_anchor]
    repl_c3d_rot += step[None, :]
    repl_geom.coords = repl_c3d_rot

    # Create geometries without atoms that will be deleted/replaces
    geom_ = geom.get_subgeom_without([ind])
    repl_geom_ = repl_geom.get_subgeom_without([repl_ind])
    # Union
    union = geom_ + repl_geom_
    if opt:
        freeze_inds = [ind for ind, _ in enumerate(geom_.atoms)]
        union = run_opt(union, freeze_inds, use_xtb=use_xtb, charge=charge, mult=mult)
        # Update coords in 'repl_geom' with optimized coordinates
        repl_geom_.coords3d = union.coords3d[len(freeze_inds) :]

    # This is useful when multiple replacements are to be done. If multiple
    # replacements would be done sequentially on the same geometry, while
    # returning the union, would require taking into account altered atom ordering/
    # numbering. By disabling return_full we can only return the replacement fragment,
    # that can be added to the original geometry, after all replacement fragments
    # have been calculated.
    if return_full:
        return_geom = union
    else:
        return_geom = repl_geom_
    return return_geom


# Located in ../geom_library/replacements/
REPLACEMENTS = {
    "Ph": ("benzene.xyz", 8),
    "OMe": ("methanol.xyz", 4),
    "OEt": ("ethanol.xyz", 8),
    "OH": ("water.xyz", 1),
    "F": ("hf.xyz", 1),
    "Cl": ("hcl.xyz", 1),
    "Br": ("hbr.xyz", 1),
}


def normalize_replacements(replacements):
    replacements = list(replacements)
    for i, repl in enumerate(replacements):
        if isinstance(repl[1], str):
            ind, repl_key = repl
            repl_fn, repl_ind = REPLACEMENTS[repl_key]
            repl_geom = geom_loader(LIB_DIR / "replacements" / repl_fn)
            replacements[i] = (ind, repl_ind, repl_geom)
        elif len(repl) == 3:
            continue
        else:
            raise Exception("Invalid replacements!")
    return replacements


def replace_atoms(geom, replacements, opt=False, use_xtb=False, charge=0, mult=1):
    replacements = normalize_replacements(replacements)
    repl_frags = list()
    del_inds = list()
    for ind, repl_ind, repl_geom in replacements:
        repl_frag = replace_atom(geom, ind, repl_geom, repl_ind, False, opt=opt, use_xtb=use_xtb)
        repl_frags.append(repl_frag)
        del_inds.append(ind)
    org_atom_num = len(geom.atoms)
    union = geom.get_subgeom_without(del_inds)
    for frag in repl_frags:
        union += frag

    if opt:
        freeze_inds = [ind for ind in range(org_atom_num - len(del_inds))]
        union = run_opt(union, freeze_inds, use_xtb=use_xtb)

    return union


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("geom_fn")
    parser.add_argument("replacements", nargs="+")
    parser.add_argument("--opt", action="store_true")
    parser.add_argument("--jmol", action="store_true")
    parser.add_argument("--out_fn", type=str, default="union.xyz")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--mult", type=int, default=1)
    parser.add_argument("--xtb", action="store_true")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    geom_fn = args.geom_fn
    replacements_ = args.replacements
    opt = args.opt
    jmol = args.jmol
    out_fn = args.out_fn
    charge = args.charge
    mult = args.mult
    xtb = args.xtb

    # Activate opt when xtb is given
    opt = opt or xtb

    assert len(replacements_) % 2 == 0
    replacements = [(int(ind), key) for ind, key in chunks(replacements_, 2)]

    geom = geom_loader(geom_fn)
    union = replace_atoms(geom, replacements, opt=opt, charge=charge, mult=mult, use_xtb=xtb)
    print(union.as_xyz())

    if union:
        with open(out_fn, "w") as handle:
            handle.write(union.as_xyz())
        print(f"Saved new geometry to '{out_fn}'.")
    if jmol:
        union.jmol()

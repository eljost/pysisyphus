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
from pysisyphus.helpers import align_coords
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.xyzloader import coords_to_trj


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


def get_rot_mat(coords3d_1, coords3d_2):
    coords3d_1 = coords3d_1.copy()
    coords3d_2 = coords3d_2.copy()
    centroid_1 = coords3d_1.mean(axis=0)
    centroid_2 = coords3d_2.mean(axis=0)
    coords3d_1 -= centroid_1[None, :]
    coords3d_2 -= centroid_2[None, :]

    tmp = coords3d_2.T.dot(coords3d_1)
    U, W, Vt = np.linalg.svd(tmp)
    rot_mat = U.dot(Vt)
    if np.linalg.det(rot_mat) < 0:
        U[:, -1] *= -1
        rot_mat = U.dot(Vt)
    return coords3d_1, coords3d_2, rot_mat


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

    *_, rot_mat = get_rot_mat(coords3d_1, coords3d_2_subset)

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

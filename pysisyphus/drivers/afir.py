# [1] https://doi.org/10.1002/jcc.23481
#     Exploring transition state structures for intramolecular pathways
#     by the artificial force induced reaction method
#     Maeda, Morokuma et al, 2013

import itertools as it
from functools import reduce
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist
from scipy.optimize import least_squares

from pysisyphus.calculators.AFIR import AFIR, AFIRPath
from pysisyphus.config import OUT_DIR_DEFAULT
from pysisyphus.cos.NEB import NEB
from pysisyphus.drivers.opt import run_opt
from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import pick_image_inds
from pysisyphus.intcoords.helpers import get_bond_difference
from pysisyphus.intcoords.setup_fast import find_bonds


##########################
#                        #
#  Multi-component AFIR  #
#       MC-AFIR          #
#                        #
##########################


##########################
#                        #
#  Multi-component AFIR  #
#       MC-AFIR          #
#                        #
##########################


def generate_random_union(geoms, offset=1.0, copy=True):
    assert 2 <= len(geoms) <= 6
    # Center, rotate and displace from origin acoording to approximate radius
    # and an offset.
    # Displace along +x, -x, +y, -y, +z, -z. Affords at max 6 fragments.
    axis_inds = (0, 0, 1, 1, 2, 2)
    randomized = list()
    for i, geom in enumerate(geoms):
        if copy:
            geom = geom.copy()
        geom.center()
        geom.rotate()
        step_size = geom.approximate_radius() + offset
        step = np.zeros(3)
        axis = axis_inds[i]
        # Alternate between negative and positive direction along x/y/z
        step[axis] = (-1) ** i * step_size
        geom.coords3d += step[None, :]
        randomized.append(geom)
    union = reduce(lambda geom1, geom2: geom1 + geom2, randomized)
    return union


def relax_afir_path(atoms, cart_coords, calc_getter, images=15, out_dir=None):
    image_inds = pick_image_inds(cart_coords, images=images)
    images = [Geometry(atoms, cart_coords[i]) for i in image_inds]

    # Relax last image
    opt_kwargs = {
        "dump": True,
        "prefix": "last",
        "out_dir": out_dir,
    }
    last_image = images[-1]
    last_image_backup = last_image.copy()
    run_opt(last_image, calc_getter, opt_key="rfo", opt_kwargs=opt_kwargs)
    _, broken = get_bond_difference(last_image, last_image_backup)
    if broken:
        return

    cos_kwargs = {}
    cos = NEB(images, **cos_kwargs)
    cos_opt_kwargs = {
        "align": True,
        "dump": True,
        "max_cycles": 30,
        "out_dir": out_dir,
    }
    run_opt(cos, calc_getter, opt_key="lbfgs", opt_kwargs=cos_opt_kwargs)


def run_mc_afir_paths(geoms, calc_getter, num=5, gamma=None, t=None, T=298.15):
    assert gamma or t, "Either parameter gamme of time t must be given!"

    # Set up list of fragments
    i = 0
    fragments = list()
    for geom in geoms:
        atom_num = len(geom.atoms)
        fragments.append(np.arange(atom_num) + i)
        i += atom_num

    unions = [generate_random_union(geoms, copy=True) for _ in range(num)]
    # with open("unions.trj", "w") as handle:
    # handle.write("\n".join([union.as_xyz() for union in unions]))

    afir_paths = list()
    for i, union in enumerate(unions):
        out_dir = Path(f"out_{i:03d}")

        union_backup = union.copy()
        # Preoptimize union, so we start from an actual minimum
        union_opt_kwargs = {
            "prefix": "union",
            "out_dir": out_dir,
        }
        run_opt(union, calc_getter, opt_key="rfo", opt_kwargs=union_opt_kwargs)

        afir_path = multicomponent_afir(
            union, calc_getter, fragments, gamma=gamma, out_dir=out_dir
        )
        afir_paths.append(afir_path)
        formed, broken = get_bond_difference(union, union_backup)
        if formed or broken:
            relax_afir_path(
                union.atoms, afir_path.cart_coords, calc_getter, out_dir=out_dir
            )


def multicomponent_afir(geom, calc_getter, fragments, gamma, out_dir=None):
    actual_calc = calc_getter(out_dir=out_dir / OUT_DIR_DEFAULT)

    def afir_calc_getter():
        afir_calc = AFIR(
            actual_calc, fragment_indices=fragments, gamma=gamma, out_dir=out_dir
        )
        return afir_calc

    opt_kwargs = {
        "dump": True,
        "out_dir": out_dir,
        "prefix": "afir",
        "max_cycles": 200,
    }
    opt_result = run_opt(geom, afir_calc_getter, opt_key="rfo", opt_kwargs=opt_kwargs)
    opt = opt_result.opt

    afir_path = AFIRPath(
        atoms=geom.atoms,
        cart_coords=np.array(opt.cart_coords),
        energies=np.array(opt.true_energies),
        forces=np.array(opt.true_forces),
        opt_is_converged=opt.is_converged,
    )

    return afir_path


def analyze_afir_path(energies):
    energies = np.array(energies)
    # min_ind = energies.argmin()
    max_ind = energies.argmax()

    local_minima = list()
    local_maxima = list()
    stationary_points = list()
    for i, en in enumerate(energies[1:-1], 1):
        prev_en = energies[i - 1]
        next_en = energies[i + 1]
        if is_minimum := prev_en > en < next_en:
            local_minima.append(i)
        if is_ts := prev_en < en > next_en:
            local_maxima.append(i)
        if is_minimum or is_ts:
            stationary_points.append(i)

    if max_ind in stationary_points:
        ts_ind = np.where(stationary_points == max_ind)[0][0]
        prev_min_ind = stationary_points[ts_ind - 1]
        next_min_ind = stationary_points[ts_ind + 1]
        sp_inds = [prev_min_ind, max_ind, next_min_ind]
    else:
        sp_inds = list()
    return sp_inds


###########################
#                         #
#  Single-component AFIR  #
#        SC-AFIR          #
#                         #
###########################


def decrease_distance(coords3d, m, n, frac=0.9):
    c3d_m = coords3d[m]
    c3d_n = coords3d[n]
    dist_vec = c3d_n - c3d_m
    step = (1 - frac) / 2 * dist_vec
    c3d_new = coords3d.copy()
    c3d_new[m] += step
    c3d_new[n] -= step
    return c3d_new


def lstsqs_with_reference(coords3d, ref_coords3d, freeze_atoms=None):
    """Least-squares w.r.t. reference coordinates while keeping some
    atoms frozen."""

    if freeze_atoms is None:
        freeze_atoms = []
    else:
        freeze_atoms = list(freeze_atoms)
    ref_dists = pdist(ref_coords3d)

    mask = np.ones_like(coords3d[:, 0], dtype=bool)
    mask[freeze_atoms] = False
    # All atoms w/o the frozen atoms
    coords = coords3d[mask].flatten()
    x0 = coords

    coords_full = coords3d.copy()

    def fun(x):
        # Consider all distances, including distances to the fixed atoms 'm' and 'n'.
        coords_full[mask] = x.reshape(-1, 3)
        dists = pdist(coords_full)
        return dists - ref_dists

    res = least_squares(fun, x0)
    opt_coords = res.x
    coords_full[mask] = opt_coords.reshape(-1, 3)
    return res, coords_full


def weight_function(atoms, coords3d, i, j, p=6):
    cr_sum = sum([CR[atoms[k].lower()] for k in (i, j)])
    r_ij = np.linalg.norm(coords3d[i] - coords3d[j])
    omega = (cr_sum / r_ij) ** 6
    return omega


def find_candidates(center, bond_sets):
    center_set = {
        center,
    }
    bonded_to_center = list()
    for bond in bond_sets:
        if center not in bond:
            continue
        bonded_to_center.append(*set(bond) - center_set)
    return bonded_to_center


def automatic_fragmentation(
    atoms, coords3d, frag1, frag2, cycles=2, p=6, bond_factor=1.25
):
    """Automatic fragmentation scheme as described in SC-AFIR paper [1]."""

    frag1 = set(frag1)
    frag2 = set(frag2)

    def w(m, n):
        """Shortcut for weight function"""
        return weight_function(atoms, coords3d, m, n, p=p)

    pairs = list(it.product(frag1, frag2))
    weights = [w(m, n) for m, n in pairs]
    max_weight = max(weights)

    bonds = find_bonds(atoms, coords3d, bond_factor=bond_factor)
    bond_sets = [set(bond) for bond in bonds.tolist()]

    def filter_candidates(candidates, partners, max_weight):
        to_keep = set()
        for candidate in candidates:
            for partner in partners:
                weight = w(candidate, partner)
                if weight > max_weight:
                    break
            else:
                to_keep.add(candidate)
        return to_keep

    def grow_fragment(frag1, frag2):
        # Find candidates that are bonded to atoms in frag1. Step 2 in [1].
        candidates = [find_candidates(m, bond_sets) for m in frag1]
        # Filter out candidates that are already contained in frag1
        candidates = [c for c in it.chain(*candidates) if c not in frag1]
        # Filter out candidates with weights that are too big. Step 3/4 in [1].
        candidates = set(filter_candidates(candidates, frag2, max_weight))
        return candidates

    for _ in range(cycles):
        f1_candidates = grow_fragment(frag1, frag2)
        f2_candidates = grow_fragment(frag2, frag1)

        # Step 5 in [1].
        f1_candidates = filter_candidates(f1_candidates, f2_candidates, max_weight)
        f2_candidates = filter_candidates(f2_candidates, f1_candidates, max_weight)

        # Step 6 in [1].
        frag1.update(f1_candidates)
        frag2.update(f2_candidates)
    return frag1, frag2


def prepare_single_component_afir(geom, m, n, calc_getter, afir_kwargs):
    """Create perturbed geometry, determine fragments and set AFIR calculator."""
    # Move target atoms closer together along distance vector (decrease distance)
    tmp_coords3d = decrease_distance(geom.coords3d, m, n)
    # Optimize remaining coordinates using least-squares, while keeping target
    # atom pair fixed.
    _, opt_coords3d = lstsqs_with_reference(tmp_coords3d, geom.coords3d, (m, n))
    # Determine fragments, using the automated fragmentation
    frag1, frag2 = automatic_fragmentation(geom.atoms, opt_coords3d, [m], [n])
    fragment_indices = [list(frag) for frag in (frag1, frag2)]
    # Set lstsq-optimized coordinates and created wrapped calculator
    geom.coords3d = opt_coords3d
    calc = calc_getter()
    afir_defaults = {
        "calculator": calc,
        "fragment_indices": fragment_indices,
        "complete_fragments": False,
        "ignore_hydrogen": False,
    }
    afir_kwargs = afir_kwargs.copy()
    # Force the use of the determined fragment_indices, by override any
    # potential user input.
    afir_kwargs.update(afir_defaults)
    afir_calc = AFIR(**afir_kwargs)
    geom.set_calculator(afir_calc)

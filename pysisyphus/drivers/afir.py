# [1] https://doi.org/10.1002/jcc.23481
#     Exploring transition state structures for intramolecular pathways
#     by the artificial force induced reaction method
#     Maeda, Morokuma et al, 2013
# [2] https://pubs.acs.org/doi/10.1021/ct200290m
#     Finding Reaction Pathways of Type A + B → X: Toward Systematic
#     Prediction of Reaction Mechanisms
#     Maeda, Morokuma, 2011

from dataclasses import dataclass
import itertools as it
import logging
from functools import reduce
import os
from pathlib import Path
from pprint import pformat
import shutil
import traceback
from typing import Callable, Dict, List, Tuple, Optional

import numpy as np
from numpy.typing import NDArray
from rmsd import kabsch_rmsd
from scipy.spatial.distance import pdist
from scipy.optimize import least_squares

from pysisyphus import logger as pysis_logger
from pysisyphus.calculators.AFIR import AFIR, CovRadiiSumZero
from pysisyphus.calculators import HardSphere
from pysisyphus.config import AFIR_RMSD_THRESH, OUT_DIR_DEFAULT
from pysisyphus.constants import BOHR2ANG
from pysisyphus.cos.NEB import NEB
from pysisyphus.drivers.opt import run_opt
from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import pick_image_inds, check_for_end_sign
from pysisyphus.helpers_pure import to_sets
from pysisyphus.intcoords.helpers import get_bond_difference
from pysisyphus.intcoords.setup import get_pair_covalent_radii
from pysisyphus.intcoords.setup_fast import find_bonds
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.xyzloader import make_xyz_str


logger = pysis_logger.getChild("afir")
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler("afir.log", mode="w", delay=True)
logger.addHandler(file_handler)


AFIR_BOND_FACTOR = 1.2


@dataclass
class AFIRPath:
    atoms: tuple
    cart_coords: np.ndarray
    energies: np.ndarray
    forces: np.ndarray
    charge: int
    mult: int
    opt_is_converged: Optional[bool] = None
    gamma: Optional[float] = None
    path_indices: Optional[List[int]] = None

    def compatible(self, other):
        check = ("atoms", "charge", "mult")
        return all([getattr(self, name) == getattr(other, name) for name in check])

    def __add__(self, other):
        """This assumes that the first item in other.cart_coords is the same as the
        last item of self.cart_coords."""
        assert self.compatible(other)

        def conc(name):
            return np.concatenate(
                (getattr(self, name), getattr(other, name)[1:]), axis=0
            )

        cart_coords = conc("cart_coords")
        energies = conc("energies")
        forces = conc("forces")

        if self.path_indices is None:
            path_indices = [0] * len(self.cart_coords)
        else:
            path_indices = self.path_indices
        path_indices += [path_indices[-1] + 1] * (len(other.cart_coords) - 1)

        return AFIRPath(
            atoms=self.atoms,
            cart_coords=cart_coords,
            energies=energies,
            forces=forces,
            charge=self.charge,
            mult=self.mult,
            path_indices=path_indices,
        )

    def dump_trj(self, fn):
        geom = Geometry(self.atoms, self.cart_coords[0])
        xyzs = [geom.as_xyz(cart_coords=cc) for cc in self.cart_coords]
        with open(fn, "w") as handle:
            handle.write("\n".join(xyzs))


def analyze_afir_path(energies):
    energies = np.array(energies)
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


##########################
#  Multi-component AFIR  #
#        Helper          #
##########################


def generate_random_union(geoms, offset=2.0, rng=None):
    """Unite fragments into one Geometry with random fragment orientations.

    Center, rotate and translate from origin acoording to approximate radius
    and an offset.
    Displace along +x, -x, +y, -y, +z, -z.

    Results for > 3 fragments don't look so pretty ;).
    """

    axis_inds = (0, 0, 1, 1, 2, 2)
    axis_inds_num = len(axis_inds)
    axis_translations = np.zeros(axis_inds_num)
    randomized = list()
    for i, geom in enumerate(geoms):
        i_mod = i % axis_inds_num
        geom = geom.copy()
        geom.center()
        geom.rotate(rng=rng)
        axis = axis_inds[i_mod]
        step_size = axis_translations[i_mod] + geom.approximate_radius() + offset
        axis_translations[i_mod] = step_size
        step = np.zeros(3)
        # Alternate between negative and positive direction along x/y/z
        step[axis] = (-1) ** i * step_size
        geom.coords3d += step[None, :]
        randomized.append(geom)
    union = reduce(lambda geom1, geom2: geom1 + geom2, randomized)
    return union


def generate_random_union_ref(geoms, rng=None, opt_kwargs=None):
    """Unite fragments into one Geometry with random fragment orientations."""

    geoms = [geom.copy() for geom in geoms]

    if rng is None:
        rng = np.random.default_rng()

    if opt_kwargs is None:
        opt_kwargs = {}

    # Random rotations
    for geom in geoms:
        geom.center()
        geom.rotate(rng=rng)

    # HardSphere optimization to fix overlapping fragments
    #
    # Set up fragment lists.
    fragments = list()
    for geom in geoms:
        geom_inds = np.arange(len(geom.atoms))
        try:
            geom_inds += fragments[-1][-1] + 1
        except IndexError:
            pass
        fragments.append(geom_inds.tolist())
    union = reduce(lambda geom1, geom2: geom1 + geom2, geoms)
    calc = HardSphere(union, frags=fragments, permutations=False, kappa=1.0)
    union.set_calculator(calc)
    _opt_kwargs = {
        "max_step": 0.2,
        "max_cycles": 500,
    }
    _opt_kwargs.update(opt_kwargs)
    opt = FIRE(union, **_opt_kwargs)
    opt.run()
    if not opt.is_converged:
        union = None
    else:
        # Remove HardSphere calculator
        union.clear()
        del union.calculator

    return union


def prepare_mc_afir(geoms, rng=None, **kwargs):
    """Wrapper for generate_random_union(_ref)."""
    union = generate_random_union(geoms, rng=rng, **kwargs)

    # Set up list of fragments
    i = 0
    fragments = list()
    for geom in geoms:
        atom_num = len(geom.atoms)
        fragments.append(np.arange(atom_num) + i)
        i += atom_num

    afir_kwargs = {
        "fragment_indices": fragments,
    }
    broken_bonds = []
    return union, afir_kwargs, broken_bonds


###########################
#  Single-component AFIR  #
#        Helper           #
###########################


def decrease_distance(coords3d, m, n, frac=0.8):
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
    omega = (cr_sum / r_ij) ** p
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
                if candidate == partner:
                    break
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
        assert frag1.isdisjoint(
            frag2
        ), "Overlapping fragments detected!"  # Sanity check
    return frag1, frag2


def prepare_sc_afir(geom, m, n, bond_factor=AFIR_BOND_FACTOR):
    """Create perturbed geometry, determine fragments and set AFIR calculator."""
    geom = geom.copy()
    atoms = geom.atoms
    org_coords3d = geom.coords3d.copy()

    def bond_sets(coords3d):
        bonds = find_bonds(atoms, coords3d, bond_factor=bond_factor)
        return set([frozenset(bond) for bond in bonds.tolist()])

    org_bond_sets = bond_sets(org_coords3d)

    # Move target atoms closer together along distance vector (decrease distance)
    decr_coords3d = decrease_distance(geom.coords3d, m, n)
    # Optimize remaining coordinates using least-squares, while keeping target
    # atom pair fixed.
    _, opt_coords3d = lstsqs_with_reference(decr_coords3d.copy(), geom.coords3d, (m, n))
    # Determine fragments, using the automated fragmentation
    frag1, frag2 = automatic_fragmentation(atoms, opt_coords3d, [m], [n])
    logger.debug(f"Fragments for target pair [{m}, {n}]: ({frag1}, {frag2})")
    fragment_indices = [list(frag) for frag in (frag1, frag2)]
    # Set lstsq-optimized coordinates and created wrapped calculator
    geom.coords3d = opt_coords3d

    opt_bond_sets = bond_sets(opt_coords3d)
    broken_bonds = org_bond_sets - opt_bond_sets

    afir_kwargs = {
        "fragment_indices": fragment_indices,
        # If 'complete_fragments' would be True, all remaining atom indices not
        # present in 'fragment_indices' would be assigned to a third fragment.
        "complete_fragments": False,
    }

    def set_atoms(inds, atom_type="X", mod_atoms=None):
        if mod_atoms is None:
            mod_atoms = list(atoms)
        for i in inds:
            mod_atoms[i] = atom_type
        return mod_atoms

    atoms_target = set_atoms((m, n))
    atoms_fragments = set_atoms(frag1)
    atoms_fragments = set_atoms(frag2, atom_type="Q", mod_atoms=atoms_fragments)

    atoms_coords3d = {
        "original": (atoms, org_coords3d),
        "original w/ target atoms": (atoms_target, org_coords3d),
        "decreased distance": (atoms, decr_coords3d),
        "decreased distance w/ target atoms": (atoms_target, decr_coords3d),
        "lstsq optimized": (atoms, opt_coords3d),
        "lstsq optimized w/ target atoms": (atoms_target, opt_coords3d),
        "lstsq optimized w/ fragments": (atoms_fragments, opt_coords3d),
    }
    trj = "\n".join(
        [
            make_xyz_str(atoms, BOHR2ANG * coords3d, comment=key)
            for key, (atoms, coords3d) in atoms_coords3d.items()
        ]
    )
    return geom, afir_kwargs, broken_bonds, trj


def determine_target_pairs(
    atoms: Tuple[str],
    coords3d: NDArray,
    min_: float = 1.25,
    max_: float = 5.0,
    active_atoms=None,
) -> List[Tuple[int]]:
    """Determine possible target m, n atom pairs for SC-AFIR calculations."""
    if active_atoms is None:
        active_atoms = range(len(atoms))
    active_atoms = set(active_atoms)

    pair_cov_radii = get_pair_covalent_radii(atoms)
    pair_dists = pdist(coords3d)
    quots = pair_dists / pair_cov_radii
    pair_inds = it.combinations(range(len(atoms)), 2)
    target_pairs = list()
    for pair_ind, quot in zip(pair_inds, quots):
        if (min_ <= quot <= max_) and (set(pair_ind) & active_atoms):
            target_pairs.append(pair_ind)
    return target_pairs


def determine_target_pairs_for_geom(geom: Geometry, **kwargs) -> List[Tuple[int]]:
    """Determine possible target m, n atom pairs for SC-AFIR calculations
    from geom."""
    target_pairs = determine_target_pairs(geom.atoms, geom.coords3d, **kwargs)
    return target_pairs


def coordinates_similar(
    test_coords3d: NDArray, ref_coords3d: List[NDArray], rmsd_thresh: float = 1e-2
) -> Tuple[bool, int]:
    # When the reference coordinates are an empty list.
    if len(ref_coords3d) == 0:
        return False, -1
    test_centered3d = test_coords3d - test_coords3d.mean(axis=0)[None, :]
    for i, rcoords3d in enumerate(ref_coords3d):
        ref_centered3d = rcoords3d - rcoords3d.mean(axis=0)[None, :]
        rmsd_ = kabsch_rmsd(test_centered3d, ref_centered3d)
        if rmsd_ <= rmsd_thresh:
            break
    else:
        return False, -1
    return True, i


def geom_similar(test_geom: Geometry, ref_geoms: List[Geometry], **kwargs) -> bool:
    return coordinates_similar(
        test_geom.coords3d, [geom.coords3d for geom in ref_geoms], **kwargs
    )


########################
#  Actual AFIR drivers #
########################


def opt_afir_path(geom, calc_getter, afir_kwargs, opt_kwargs=None, out_dir=None):
    """Minimize geometry with AFIR calculator."""
    if opt_kwargs is None:
        opt_kwargs = dict()
    if out_dir is None:
        out_dir = "."
    out_dir = Path(out_dir)

    actual_calc = calc_getter(out_dir=out_dir / OUT_DIR_DEFAULT)

    def afir_calc_getter():
        afir_calc = AFIR(actual_calc, out_dir=out_dir, **afir_kwargs)
        return afir_calc

    _opt_kwargs = {
        "dump": True,
        "out_dir": out_dir,
        "prefix": "afir",
        "max_cycles": 125,
        "overachieve_factor": 3,
        "hessian_update": "flowchart",
    }
    logger.debug(
        "\n".join(
            (
                "afir_kwargs:",
                "\t" + pformat(afir_kwargs),
                "opt_kwargs:",
                "\t" + pformat(_opt_kwargs),
            )
        )
    )
    _opt_kwargs.update(opt_kwargs)
    opt_result = run_opt(geom, afir_calc_getter, opt_key="rfo", opt_kwargs=_opt_kwargs)
    opt = opt_result.opt

    afir_path = AFIRPath(
        atoms=geom.atoms,
        cart_coords=np.array(opt.cart_coords),
        energies=np.array(opt.true_energies),
        forces=np.array(opt.true_forces),
        opt_is_converged=opt.is_converged,
        charge=actual_calc.charge,
        mult=actual_calc.mult,
        gamma=afir_kwargs["gamma"],
    )

    return afir_path


def run_afir_path(
    geom,
    calc_getter,
    out_dir,
    gamma_max,
    gamma_interval: Tuple[float, float],
    rng,
    ignore_bonds=None,
    bond_factor=AFIR_BOND_FACTOR,
    afir_kwargs=None,
    opt_kwargs=None,
):
    """Driver for AFIR minimizations with increasing gamma values."""
    if ignore_bonds is None:
        ignore_bonds = list()
    if afir_kwargs is None:
        afir_kwargs = dict()
    if opt_kwargs is None:
        opt_kwargs = {}

    if out_dir.exists():
        dir_contents = os.listdir(out_dir)
        for fn in dir_contents:
            fn = out_dir / fn
            try:
                os.remove(fn)
            except IsADirectoryError:
                shutil.rmtree(fn)
    else:
        os.mkdir(out_dir)

    # Decreasing the distance between two atoms in SC-AFIR may lead to broken
    # bonds for these two bonds, e.g., hydrogen atoms "that are left behind".
    # These bonds will be formed again when the AFIR function is minimized, but
    # we are not interested in these changes. So they can be ignored here.
    ignore_bonds = set([frozenset(bond) for bond in ignore_bonds])

    ref_calc = calc_getter()
    ref_energy = ref_calc.get_energy(geom.atoms, geom.cart_coords)["energy"]
    geom_backup = geom.copy()

    # By using (1.0, 1.0) as interval we can directly start at gamma_max, e.g.,
    # in SC-AFIR.
    gamma_low, gamma_high = gamma_interval
    assert gamma_high >= gamma_low
    gamma_spread = gamma_high - gamma_low
    gamma_0 = (gamma_low + (gamma_spread * rng.random(1)[0])) * gamma_max
    gamma_inrc = 0.1 * gamma_max

    afir_paths = list()
    best_afir_path = None
    lowest_barrier = None
    gamma = gamma_0

    # Minimize AFIR functions until gamma exceeds gamma_max.
    logger.info(f"New AFIR run with γ_max={gamma_max:.6f} au")
    while gamma <= gamma_max:
        gamma_ratio = gamma / gamma_max
        logger.info(f"AFIR run with γ={gamma:.6f} au, γ/γ_max={gamma_ratio: >6.2%}")
        _afir_kwargs = afir_kwargs.copy()
        _afir_kwargs["gamma"] = gamma
        try:
            afir_path = opt_afir_path(
                geom,
                calc_getter,
                afir_kwargs=_afir_kwargs,
                opt_kwargs=opt_kwargs,
                out_dir=out_dir,
            )
        # Can happen in SC-AFIR runs when fragments comprise only hydrogens.
        except CovRadiiSumZero:
            logger.warning("Sum of covalent radii is 0.0!")
            best_afir_path = None
            break
        except Exception:
            logger.error(f"Optimization crashed!\n{traceback.format_exc()}")
            best_afir_path = None
            break

        afir_paths.append(afir_path)
        true_energies = afir_path.energies

        # Check for changes in bond topology by comparing the result of the current
        # minimization to the initial geometry.
        formed, broken = get_bond_difference(geom_backup, geom, bond_factor=bond_factor)
        formed = to_sets(formed) - ignore_bonds
        broken = to_sets(broken) - ignore_bonds
        if formed or broken:
            max_energy = true_energies.max()
            barrier = max_energy - ref_energy

            # Always store the first path that leads to a change in bond topology.
            if lowest_barrier is None:
                best_afir_path = afir_path
                lowest_barrier = barrier
            # Break when a path with a lower or higher barrier is detected. This
            # should happen in the cycle after the first change in bond toplogy was
            # detected.
            if barrier < lowest_barrier:
                best_afir_path = afir_path
                lowest_barrier = barrier
                break
            elif barrier > lowest_barrier:
                break
        else:
            barrier = None

        # Update values for next cycle
        gamma += gamma_inrc

    # reduce(lambda ap1, ap2: ap1 + ap2, afir_paths).dump_trj("dumped.trj")

    # Construct TS guess from highest energy point along the AFIR path.
    if best_afir_path:
        guess_ind = best_afir_path.energies.argmax()
        guess_coords = best_afir_path.cart_coords[guess_ind]
        ts_guess = Geometry(geom.atoms, guess_coords)
        ts_guess.dump_xyz(out_dir / "ts_guess.xyz")
        afir_path_merged = reduce(lambda ap1, ap2: ap1 + ap2, afir_paths)
    else:
        ts_guess = None
        afir_path_merged = None
    return ts_guess, afir_path_merged


def relax_afir_path(atoms, cart_coords, calc_getter, images=15, out_dir=None):
    """Sample imagef from AFIR path and do COS relaxation."""
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


def run_mc_afir_paths(
    geoms: List,
    calc_getter: Callable,
    gamma_max: float,
    rng,
    N_max: int = 5,
    gamma_interval: Tuple[float, float] = (0.0, 1.0),
    afir_kwargs: Optional[Dict] = None,
    opt_kwargs: Optional[Dict] = None,
):
    if afir_kwargs is None:
        afir_kwargs = dict()

    N = 0
    N_0 = 0
    fmt = " >4d"
    while True:
        geom, _afir_kwargs, *_ = prepare_mc_afir(geoms, rng=rng)
        afir_kwargs = afir_kwargs.copy()
        afir_kwargs.update(_afir_kwargs)
        out_dir = Path(f"out_{N:03d}")

        ts_guess, afir_path = run_afir_path(
            geom,
            calc_getter,
            out_dir,
            gamma_max,
            gamma_interval,
            rng=rng,
            afir_kwargs=afir_kwargs,
            opt_kwargs=opt_kwargs,
        )
        ts_guess_is_new = yield N, ts_guess, afir_path
        if ts_guess_is_new:
            N_0 = N
        # Allow up to consecutive N_max failures
        if (N - N_0) > N_max:
            break
        logger.debug(f"{N_0=:{fmt}}, {N=:{fmt}}, {N-N_0=:{fmt}}, {N_max=:{fmt}}, ")
        N += 1


def run_sc_afir_paths(
    geom: Geometry,
    calc_getter: Callable,
    gamma_max: float,
    rng,
    N_max: int = 5,
    N_sample: int = 0,
    gamma_interval: Tuple[float, float] = (1.0, 1.0),  # Start directly with gamma_max
    afir_kwargs: Optional[Dict] = None,
    opt_kwargs: Optional[Dict] = None,
    target_pairs: Optional[List] = None,
):
    if target_pairs is None:
        target_pairs = determine_target_pairs_for_geom(geom)

    if afir_kwargs is None:
        afir_kwargs = dict()
    if opt_kwargs is None:
        opt_kwargs = dict()
    _opt_kwargs = opt_kwargs.copy()

    i = 0
    while len(target_pairs) > 0:
        m, n = target_pairs.pop(0)
        logger.info(f"Running SC-AFIR with target_pair ({m}, {n}).")
        # _afir_kwargs will contain the automatically determined fragments
        geom_mod, _afir_kwargs, broken_bonds, trj = prepare_sc_afir(geom, m, n)
        afir_kwargs = afir_kwargs.copy()
        afir_kwargs.update(_afir_kwargs)

        # _opt_kwargs.update({
        # "fragments": [[m, ], [n, ]],
        # "monitor_frag_dists": 5,
        # })

        out_dir = Path(f"out_{i:03d}")
        ts_guess, afir_path = run_afir_path(
            geom_mod,
            calc_getter,
            out_dir,
            gamma_max,
            gamma_interval,
            rng=rng,
            ignore_bonds=broken_bonds,
            afir_kwargs=afir_kwargs,
            opt_kwargs=_opt_kwargs,
        )
        ts_guess_is_new = yield i, ts_guess, afir_path
        i += 1
        # Here, TS optimization & IRC integration could take place, when the TS is new.


def run_afir_paths(
    afir_key,
    geoms,
    calc_getter,
    afir_kwargs=None,
    opt_kwargs=None,
    seed=None,
    N_sample=None,
    rmsd_thresh: float = AFIR_RMSD_THRESH,
    **kwargs,
):
    if seed is None:
        rng = np.random.default_rng()
        seed = rng.integers(1_000_000_000_000)
    logger.info(f"{seed=}")
    rng = np.random.default_rng(seed)

    if afir_key == "sc":
        assert (
            len(geoms) == 1
        ), f"Expected only 1 geometry for SC-AFIR, but got {len(geoms)}!."
        geoms = geoms[0]

    afir_funcs = {
        "mc": run_mc_afir_paths,
        "sc": run_sc_afir_paths,
    }
    afir_func = afir_funcs[afir_key]
    afir_coroutine = afir_func(
        geoms,
        calc_getter,
        rng=rng,
        afir_kwargs=afir_kwargs,
        opt_kwargs=opt_kwargs,
        **kwargs,
    )
    ts_guesses = list()
    afir_paths = list()
    stop_sign = "afir_stop"

    i, ts_guess, afir_path = next(afir_coroutine)
    while True:
        prefix = f"{i:03d}"
        logger.info(f"AFIR run {prefix}")
        # Check if we found enough TS guesses
        if N_sample and (len(ts_guesses) >= N_sample):
            break

        # Check similarity of TS guess to previously found TS guesses
        try:
            rmsds = [ts_guess.rmsd(guess) for guess in ts_guesses]
            min_rmsd = min(rmsds)
        except AttributeError:  # Raised when ts_guess is None
            min_rmsd = None
        except ValueError:
            min_rmsd = rmsd_thresh if (ts_guess is not None) else None

        if ts_guess_is_new := ((min_rmsd is not None) and (min_rmsd >= rmsd_thresh)):
            logger.debug(f"rmsds different enough! {min_rmsd=:.6f} au")
            ts_guesses.append(ts_guess)
            afir_paths.append(afir_path)
            afir_path.dump_trj(f"{prefix}_afir_path.trj")
            ts_guess.dump_xyz(f"{prefix}_ts_guess.xyz")
        elif min_rmsd:
            logger.debug(f"rmsds too similar! {min_rmsd=:.6f} au)")

        sign = check_for_end_sign(add_signs=(stop_sign,))
        if sign == stop_sign:
            break
        logger.info(f"{len(ts_guesses)=: >5d}")
        try:
            i, ts_guess, afir_path = afir_coroutine.send(ts_guess_is_new)
        except StopIteration:
            break

    return ts_guesses, afir_paths
